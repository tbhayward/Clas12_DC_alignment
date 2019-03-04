/*
 * author Timothy B. Hayward
 * 2018
 * CLAS (12) Collaboration Service Work 
 * (supervised by Mac Mestayer)
 */

import java.io.File;

import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.clas.physics.*;
import org.jlab.clas12.physics.*;

import javax.swing.JFrame;
import org.jlab.groot.data.*;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.io.hipo.HipoDataBank;

import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.fitter.DataFitter;

public class theta_studies {

	public static boolean banks_test(HipoDataEvent event){
		boolean banks_result = true; // check to see if the event has all of the banks present
		// TBHits required for residuals, TBTracks required for angular limit
		if (!(event.hasBank("TimeBasedTrkg::TBHits"))) {
        	banks_result = false;
    	} else if (!(event.hasBank("TimeBasedTrkg::TBTracks"))) {
    		banks_result = false;
    	} else if (!(event.hasBank("TimeBasedTrkg::Trajectory"))) {
    		banks_result = false;
    	} else if (!(event.hasBank("REC::Particle"))) {
    		banks_result = false;
    	} else if (!(event.hasBank("REC::Cherenkov"))) {
    		banks_result = false;
    	} else if (!(event.hasBank("REC::Calorimeter"))) {
    		banks_result = false;
    	}
    	return banks_result;
	}

	public static boolean track_check(HipoDataBank recBank, HipoDataBank ccBank, 
		HipoDataBank calBank) {
		boolean track_check = true;
		if (!(recBank.getFloat("beta", 0)>0)) { // require ftof record a beta 
			track_check=false; 
		}
		for (int ccBankRow = 0; ccBankRow<ccBank.rows(); ccBankRow++) {
			if (ccBank.getInt("pindex", ccBankRow)==0) {
				if (ccBank.getFloat("nphe", ccBankRow)<2) {
					track_check=false;
				}
			}
		}
		double cal_energy = 0;
		for (int calBankRow = 0; calBankRow<calBank.rows(); calBankRow++) {
			if (calBank.getInt("pindex", calBankRow)==0) {
				cal_energy+=calBank.getFloat("energy", calBankRow);
			}
		}
		if (cal_energy<0.06) {
			track_check=false;
		}
		return track_check;
	}

	private static double phi_calculation (float x, float y, float z) {
		return (double) Math.atan2(x,y);	
	}

	private static double theta_calculation (float x, float y, float z) {
		// convert cartesian coordinates to polar angle
		double r = Math.pow(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2),0.5);
		return (double) Math.acos(z/r);
	}

	public static F1D fit_function(int n_bins, double min_bin, double max_bin, H1F histogram) {
		// first fit to a Gaussian and check chi^2/dof 
		double[] func_limits = [0.6, 1.5, -500, 500, 0, 15000];
		String func_string = "[A0]+[A1]*x+[A2]*x*x+[A3]*x*x*x+[A4]*x*x*x*x+[A5]*x*x*x*x*x+[A6]*x*x*x*x*x*x";
		F1D func = new F1D("func",func_string, min_bin, max_bin);

		func.setParameter(0, 1); 
		func.setParameter(1, 1); 
		func.setParameter(2, 1);
		func.setParameter(3, 1); 
		func.setParameter(4, 1);  
		func.setParameter(5, 1); 
		func.setParameter(6, 1);  
	
		DataFitter.fit(func, histogram, "Q"); //No options uses error for sigma
		// double chi2dof = func.getChiSquare()/(n_bins-1);
		return func;
	}

	public static boolean keep_track(double theta) {
		def bins = [1.375, 4.125, 6.875, 9.625, 12.375, 15.125, 17.875, 20.625, 23.375, 
			26.125, 28.875, 31.625, 34.375, 37.125, 39.875, 42.625, 45,375, 48.125, 50.875, 53.625];
		def values = [57, 2927, 84437, 55560, 29554, 21354, 17850, 15520, 13278, 11782, 11098, 9403,
			5700, 3322, 2177, 1498, 727, 301, 190, 129];
		double factor = 2177; 

		if ((theta < 4.125) || (theta > 42.625)) {
			return true;
		}

		boolean theta_exceeds = false;
		int index = 0;
		while(!theta_exceeds) {
			if (theta > bins[index]) {
				index++;
			} else {
				theta_exceeds = true;
			}
		}

		Random random = new Random();
		double rand_val = (Math.abs(new Random().nextInt() % 10000)+1)/10000;
		// println(rand_val + " " + factor/values[index]);
		if (rand_val < factor/values[index]) {
			return true;
		} 

		return false;

	}

	public static void main(String[] args) {
		File[] hipo_list;
		if (args.length == 0) {
			// exits program if input directory not specified 
    	   	println("ERROR: Please enter a hipo file directory as the first argument");
    	  	System.exit(0);
    	} else {
    		File directory = new File(args[0]);
    		hipo_list = directory.listFiles();
    	}

    	int n_files;
		if ((args.length < 2)||(Integer.parseInt(args[1])>hipo_list.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			n_files = hipo_list.size();
		} else{
			// if specified, convert to int
			n_files = Integer.parseInt(args[1]);
		}
		HipoDataSource reader = new HipoDataSource();

		int n_bins = 300;
		double min_bin = 0;
		double max_bin = 40;

		H1F theta_hist = new H1F("theta",n_bins,min_bin,max_bin);

		for (int current_file; current_file<n_files; current_file++) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files);
			// limit to a certain number of files defined by n_files
			
			reader.open(hipo_list[current_file]); // open next hipo file

			while(reader.hasEvent()==true){ // cycle through events
				HipoDataEvent event = reader.getNextEvent(); // load next event in the hipo file
				if (banks_test(event)) { // check that the necessary banks are present
					HipoDataBank hitBank= (HipoDataBank) event.getBank("TimeBasedTrkg::TBHits");
					HipoDataBank trkBank= (HipoDataBank) event.getBank("TimeBasedTrkg::TBTracks");
					HipoDataBank trajBank= (HipoDataBank) event.getBank("TimeBasedTrkg::Trajectory");
					HipoDataBank recBank= (HipoDataBank) event.getBank("REC::Particle");
					HipoDataBank ccBank= (HipoDataBank) event.getBank("REC::Cherenkov");
					HipoDataBank calBank= (HipoDataBank) event.getBank("REC::Calorimeter");

					int trkID, sector;
					float theta;

					if (track_check(recBank, ccBank, calBank)&&(trkBank.rows()==1)) {
						sector = trkBank.getInt("sector", 0);
						float px = trkBank.getFloat("p0_x",0);
						float py = trkBank.getFloat("p0_y",0);
						float pz = trkBank.getFloat("p0_z",0);
						float phi = Math.toDegrees(phi_calculation(px,py,pz));
						theta = Math.toDegrees(theta_calculation(px,py,pz));
						// if (keep_track(theta)) {
						// if (true) {
						if (sector == 6) {
							theta_hist.fill(theta);
						}
						// println(Math.abs(new Random().nextInt() % 10000) + 1);
					}
				}
			}
		}

		JFrame frame_1 = new JFrame("");
		frame_1.setSize(900,470);
		EmbeddedCanvas canvas_1 = new EmbeddedCanvas();
		theta_hist.setTitleX("#theta (Deg)");
		theta_hist.setTitleY("Counts");
		theta_hist.unit();
		canvas_1.draw(theta_hist);
		frame_1.setLocationRelativeTo(null);
		frame_1.add(canvas_1); frame_1.setVisible(true);
		// F1D fits = fit_function(n_bins, min_bin, max_bin, theta_hist);
		// fits.show();
		// fits.setLineColor(2);
		// fits.setLineWidth(5); 
		// fits.setLineStyle(0);
		// canvas_1.draw(fits,"same");
		// fits.setOptStat(11111110);
		// theta_hist.setOptStat(1110);
		// canvas_1.getPad(0).getAxisX().setRange(-1999,1999);
		canvas_1.getPad(0).getAxisY().setRange(0,1);
		canvas_1.setAxisTitleSize(28);
		canvas_1.setAxisLabelSize(22);
		// canvas_1.getPad(0).setTitle("Sector 1"); 
		// canvas_1.setTitleSize(34);
	}
}