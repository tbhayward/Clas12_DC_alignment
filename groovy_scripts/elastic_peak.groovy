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
import org.jlab.io.base.DataEvent;
import org.jlab.clas.physics.LorentzVector;

import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.fitter.DataFitter;


public class elastic_peak{

	public static boolean banks_test(HipoDataEvent event) {
		boolean banks_exist = true; // check to see if the event has all of the banks present
        if (!(event.hasBank("REC::Particle"))) {
            banks_exist = false;
        } 
        return banks_exist;
	}

	public static double calculate_energy(double px, double py, double pz, double mass) {
		return Math.sqrt(px*px + py*py + pz*pz + mass*mass);
	}

	private static double theta_calculation (float x, float y, float z) {
		// convert cartesian coordinates to polar angle
		double r = Math.pow(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2),0.5);
		return (double) Math.acos(z/r);
	}

	private static double w_calculation (double proton_mass, double beam_energy, double energy,
		double theta) {
		// calculates W^2
		double first_term = proton_mass*proton_mass;
		double second_term = 2*proton_mass*(beam_energy-energy);
		double third_term = - 4*beam_energy*energy*Math.pow(Math.sin(theta/2),2);
		double W2 = first_term + second_term + third_term;
		return W2;
	}

	public static double[] bins(int n_bins, double min_bin, double max_bin) {
		double[] bin_centers = new double[n_bins-1];
		double running_center = min_bin+(0.5)*(max_bin-min_bin)/n_bins;
		for (int i=0; i<n_bins-1; i++) {
			bin_centers[i] = running_center;
			running_center+= (max_bin-min_bin)/n_bins;
		}

		return bin_centers;
	}

	public static F1D fit_function(int n_bins, double min_bin, double max_bin, H1F histogram,
		double mu, double sigma) {

		// first fit to a Gaussian and check chi^2/dof 
		String func_string = "[amp]*gaus(x,[mean],[sigma])";
		F1D func = new F1D("gauss",func_string,min_bin,max_bin);
		func.setParameter(0, 1.0);
		func.setParameter(1, mu);
		func.setParameter(2, sigma);
		DataFitter.fit(func, histogram, "Q"); //No options uses error for sigma
		double chi2dof = func.getChiSquare()/(n_bins-1);

		// add polynomials until chi2dof stops improving (<1 or >chi2dof)
		boolean chi2dof_check = true; 
		F1D test_func = new F1D("gauss+?",func_string,min_bin,max_bin);
		test_func.setParameter(0, 1.0);
		test_func.setParameter(1, mu);
		test_func.setParameter(2, sigma);
		int poly_order = 0; // current order of polynomial added to func
		while(chi2dof_check) {
		// for (int current_order=0; current_order < 3; current_order++) {
			func_string+="+[p"+poly_order+"]"
			for (int i=0; i<poly_order; i++) {
				func_string+="*x";
				test_func = new F1D("gauss+?",func_string,min_bin,max_bin);
				test_func.setParameter(i+3, 1.0);
			}
			test_func.setParameter(0, 1.0);
			test_func.setParameter(1, mu);
			test_func.setParameter(2, sigma);
			DataFitter.fit(test_func, histogram, "Q"); //No options uses error for sigma

			if ((test_func.getChiSquare()/(n_bins-1))<1) {
				println("Chi2/dof dropped below 1.0. Ending fitting routine.");
				chi2dof_check = false;
			} else if ((test_func.getChiSquare()/(n_bins-1))>chi2dof) {
				println("Chi2/dof of gaus+p(i+1) > gaus+p(i). Ending fitting routine.");
				chi2dof_check = false;
			} else {
				chi2dof = test_func.getChiSquare()/(n_bins-1);
				func = test_func;
			}
			poly_order++;	
		}
		println(chi2dof);
		return func;
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

		double min_bin;
		double max_bin;
		if (args.length > 3) {
			min_bin = Double.parseDouble(args[2]);
			max_bin = Double.parseDouble(args[3]);
		} else {
			// arguments 3 and 4 are the min and max values to bin, should be numbers
			println("Warning: Bin limits not specified. Setting to [0.65,1.2].");
			println("Setting limits strongly advised.");
			min_bin = 0.65;
			max_bin = 1.2;
		}

		int n_bins;
		if (args.length < 5) {
			// if number of bins not specified set to 250
			println("WARNING: Number of bins not specified. Set to 250 bins.");
			n_bins = 250;
		} else {
			n_bins = Integer.parseInt(args[4]);
		}

		JFrame frame = new JFrame("some frame here");
		frame.setSize(900,470);

		// initialize embedded canvas 

		double[] bin_centers = new double[n_bins-1];
		bin_centers = bins(n_bins, min_bin, max_bin);
		EmbeddedCanvas canvas = new EmbeddedCanvas();	
		H1F histogram = new H1F("",n_bins,min_bin,max_bin);
		histogram.setTitleX("W^2 (Gev^2)");
		histogram.setTitleY("Counts (Normalized)");

		for (int current_file; current_file<n_files; current_file++) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files);
			// limit to a certain number of files defined by n_files
			HipoDataSource reader = new HipoDataSource();
			reader.open(hipo_list[current_file]); // open next hipo file

			while(reader.hasEvent()==true){ // cycle through events
				HipoDataEvent event = reader.getNextEvent();
				boolean banks_exist = banks_test(event); // check to see if REC::Particle exists

				double electron_mass = 0.00051099894; // electron mass (GeV)
				double proton_mass = 0.93827208; // proton mass (GeV)
				LorentzVector beam = new LorentzVector(0, 0, 2.2, 2.2); // 2.2 GeV beam (hard coded)
				LorentzVector target = new LorentzVector(0, 0, 0, proton_mass) // LH2 target
				LorentzVector electron = new LorentzVector(0, 0, 0, 0);
				LorentzVector proton;

				if (banks_exist) {
					// load particle bank
					HipoDataBank eventBank = (HipoDataBank) event.getBank("REC::Particle"); 

					for(int particle = 0; particle < eventBank.rows(); particle++) {
						int pid = eventBank.getInt("pid", particle);
						if (pid==11) {
							float px = eventBank.getFloat("px", particle);
							float py = eventBank.getFloat("py", particle);
							float pz = eventBank.getFloat("pz", particle);
							double energy = calculate_energy(px, py, pz, electron_mass);
							if (energy>electron.e()) {
								// if new electron has higher energy, assume it is outgoing beam
								electron = new LorentzVector(px, py, pz, energy);
								double theta = theta_calculation(px,py,pz);
								double W2 = w_calculation(proton_mass, beam.e(), energy, theta);
								histogram.fill(W2);

							}
						}
					}
				}
			}
		}

		histogram.unit();
		double mu = bin_centers[histogram.getMaximumBin()];
		int half_height_index = 0;
		double half_height_test = 0;
		while(half_height_test < 0.5) {
			half_height_test = histogram.getDataY(half_height_index);
			half_height_index++;
		}
		double sigma = mu-bin_centers[half_height_index];
		println(mu+" +/- "+sigma);

		F1D func = fit_function(n_bins, min_bin, max_bin, histogram, mu, sigma);
		// draw the histogram
		canvas.draw(histogram);
		frame.add(canvas);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		canvas.setStatBoxFontSize(18);
		histogram.setOptStat(10);
		func.show();
		func.setLineColor(2); func.setLineWidth(5); func.setLineStyle(0);
		canvas.draw(func,"same");
		func.setOptStat(1110);
		println("Fit has chi2dof = "+(func.getChiSquare()/(n_bins-1)));
    }
}