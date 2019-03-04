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

public class RMSE{

	public static boolean banks_test(HipoDataEvent event){
		boolean banks_result = true; // check to see if the event has all of the banks present
		// TBHits required for residuals, TBTracks required for angular limit
		if (!(event.hasBank("TimeBasedTrkg::TBHits"))) {
        	banks_result = false;
    	} else if (!(event.hasBank("TimeBasedTrkg::TBTracks"))) {
    		banks_result = false;
    	} else if (!(event.hasBank("REC::Particle"))) {
    		banks_result = false;
    	} else if (!(event.hasBank("REC::Calorimeter"))) {
    		banks_result = false;
    	} else if (!(event.hasBank("REC::Cherenkov"))) {
    		banks_result = false;
    	} 

    	return banks_result;
	}

	private static double theta_calculation (float x, float y, float z) {
		// convert cartesian coordinates to polar angle
		double r = Math.pow(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2),0.5);
		return (double) Math.acos(z/r);
	}

	private static double round (double value, int precision) {
		// just round a number to a certain precision
   		int scale = (int) Math.pow(10, precision);
    	return (double) Math.round(value * scale) / scale;
	}

	// !! MAIN program !!
	public static void main(String[] args) {
		// create hipo file list
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

		// declare histograms
		int n_bins = 150;
		double min_bin = -0.5*10000;
		double max_bin = 0.5*10000;
		H1F residuals_histogram01 = new H1F("",n_bins,min_bin,max_bin);
		H1F residuals_histogram02 = new H1F("",n_bins,min_bin,max_bin);
		H1F residuals_histogram03 = new H1F("",n_bins,min_bin,max_bin);

		for (int current_file; current_file<n_files; current_file++) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files);
			// limit to a certain number of files defined by n_files
			HipoDataSource reader = new HipoDataSource();
			reader.open(hipo_list[current_file]); // open next hipo file

			while(reader.hasEvent()==true){ // cycle through events
				HipoDataEvent event = reader.getNextEvent(); // load next event in the hipo file

				if (banks_test(event)) {
					HipoDataBank hitBank= (HipoDataBank) event.getBank("TimeBasedTrkg::TBHits");
					HipoDataBank trkBank= (HipoDataBank) event.getBank("TimeBasedTrkg::TBTracks");
					HipoDataBank recBank= (HipoDataBank) event.getBank("REC::Particle");
					HipoDataBank ccBank= (HipoDataBank) event.getBank("REC::Cherenkov");
					HipoDataBank calBank= (HipoDataBank) event.getBank("REC::Calorimeter");

					boolean track_check = true;
					if (track_check) {
						for(int hitBankRow=0; hitBankRow<hitBank.rows(); hitBankRow++){
							// cycle over the entires in the hitbank (separate layers and regions)

							int trkID = hitBank.getInt("trkID", hitBankRow); // ID of reconstructed
							if (trkID >0) {
								float residual = 10000*hitBank.getFloat("fitResidual", hitBankRow);
								int superlayer = hitBank.getInt("superlayer", hitBankRow);
								// int layer = hitBank.getInt("layer", hitBankRow);

								if (superlayer==1||superlayer==2) {
									residuals_histogram01.fill(residual);
								} else if (superlayer==3||superlayer==4) {
									residuals_histogram02.fill(residual);
								} else {
									residuals_histogram03.fill(residual);
								}
							}
						}
					}
				}
			}
		}

		// Region 1
		double region01mean = round(residuals_histogram01.getMean(),2); 
		double region01rms = round(residuals_histogram01.getRMS(),2);
		double region01rmse = round(Math.pow(Math.pow(residuals_histogram01.getMean(),2)+
			Math.pow(residuals_histogram01.getRMS(),2),0.5),2);
		double region01entries = round(residuals_histogram01.getEntries(),2);
		println("Region 1 entries:"+region01entries+", mean: "+region01mean+", rms: "+
			region01rms+", RMSE: "+region01rmse)
		// Region 2
		double region02mean = round(residuals_histogram02.getMean(),2); 
		double region02rms = round(residuals_histogram02.getRMS(),2);
		double region02rmse = round(Math.pow(Math.pow(residuals_histogram02.getMean(),2)+
			Math.pow(residuals_histogram02.getRMS(),2),0.5),2);
		double region02entries = round(residuals_histogram02.getEntries(),2);
		println("Region 2 entries:"+region02entries+", mean: "+region02mean+", rms: "+
			region02rms+", RMSE: "+region02rmse)
		// Region 1
		double region03mean = round(residuals_histogram03.getMean(),2); 
		double region03rms = round(residuals_histogram03.getRMS(),2);
		double region03rmse = round(Math.pow(Math.pow(residuals_histogram03.getMean(),2)+
			Math.pow(residuals_histogram03.getRMS(),2),0.5),2);
		double region03entries = round(residuals_histogram03.getEntries(),2);
		println("Region 3 entries:"+region03entries+", mean: "+region03mean+", rms: "+
			region03rms+", RMSE: "+region03rmse);
		println();
		println("Mean RMSE: "+round(((region01rmse+region02rmse+region03rmse)/3),2));

		JFrame frame = new JFrame("some frame");
		frame.setSize(900,470);
		//Initialize EmbeddedCanvas
		EmbeddedCanvas canvas = new EmbeddedCanvas();
		residuals_histogram01.setLineColor(2);
		residuals_histogram02.setLineColor(4);
		residuals_histogram03.setLineColor(3);
		residuals_histogram01.setTitleX("Residuals (micrometers)");
		residuals_histogram01.setTitleY("Counts");
		canvas.draw(residuals_histogram01);
		canvas.draw(residuals_histogram02,"same")
		canvas.draw(residuals_histogram03,"same")
		frame.add(canvas);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);

	}

}
