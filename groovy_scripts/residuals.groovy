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

public class residuals{
	// groovy program to create plots / calculate stats from CLAS12 B=0 data that are of the form 
	// of those from clas_notes02/02-010.pdf by S. A. Morrow and M. D. Mestayer

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

	private static double rootMeanSquare(double... nums) {
		// calculate the root mean square of an array
        double sum = 0.0;
        for (double num : nums)
            sum += num * num;
        return Math.sqrt(sum / nums.length);
    }

    private static double round (double value, int precision) {
		// just round a number to a certain precision
    int scale = (int) Math.pow(10, precision);
    return (double) Math.round(value * scale) / scale;
	}

	private static double theta_calculation (float x, float y, float z) {
		// convert cartesian coordinates to polar angle
		double r = Math.pow(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2),0.5);
		return (double) Math.acos(z/r);
	}

	public static F1D fit_function(int n_bins, double min_bin, double max_bin, H1F histogram) {

		// first fit to a Gaussian and check chi^2/dof 
		double[] func_limits = [0.6, 1.5, -500, 500, 200, 1500];
		String func_string = "[amp]*gaus(x,[mean],[sigma])";
		F1D func = new F1D("func",func_string,min_bin,max_bin);
		func.setParameter(0, 1.0); func.setParLimits(0, func_limits[0], func_limits[1]);
		func.setParameter(1, histogram.getMean()); 
		func.setParLimits(1, func_limits[2], func_limits[3]);
		func.setParameter(2, 250); func.setParLimits(2, func_limits[4], func_limits[5]);
		DataFitter.fit(func, histogram, "Q"); //No options uses error for sigma
		double chi2dof = func.getChiSquare()/(n_bins-1);

		// add polynomials until test_chi2dof stops improving (<1 or >chi2dof)
		boolean chi2dof_check = true; 
		F1D test_func = new F1D("gauss+?",func_string,min_bin,max_bin);
		test_func.setParameter(0, 1.0); test_func.setParLimits(0, func_limits[0], func_limits[1]);
		test_func.setParameter(1, histogram.getMean()); 
		test_func.setParLimits(1, func_limits[2], func_limits[3]);
		test_func.setParameter(2, 250); test_func.setParLimits(2, func_limits[4], func_limits[5]);
		int poly_order = 0; // current order of polynomial added to func
		while(chi2dof_check) {
			func_string+="+[p"+poly_order+"]"
			for (int i=0; i<poly_order; i++) {
				func_string+="*x";
				test_func = new F1D("gauss+?",func_string,min_bin,max_bin);
				test_func.setParameter(i+3, 0.0);
			}
			test_func.setParameter(0, 1.0); 
			test_func.setParLimits(0, func_limits[0], func_limits[1]);
			test_func.setParameter(1, histogram.getMean());
			test_func.setParLimits(1, func_limits[2], func_limits[3]);
			test_func.setParameter(2, 250);
			test_func.setParLimits(2, func_limits[4], func_limits[5]);
			DataFitter.fit(test_func, histogram, "Q"); //No options uses error for sigma

			double temp_chi2dof = (test_func.getChiSquare()/(n_bins-1))

			if (temp_chi2dof<1) {
				println("Chi2/dof dropped below 1.0. Ending fitting routine.");
				chi2dof_check = false;
			} else if (temp_chi2dof>chi2dof) {
				println("Chi2/dof of gaus+p(i+1) > gaus+p(i). Ending fitting routine.");
				println(chi2dof+"    "+temp_chi2dof);
				chi2dof_check = false;
			} else {
			chi2dof = temp_chi2dof;
			func = test_func;
			}
		poly_order++;
		println(poly_order+" "+chi2dof);	
		}
		return func;
	}

	// !! MAIN program !!
	public static void main(String[] args) {
		// accepts up to 5 input arguments: 
		// args[0] = location of hipo file directory 
		// args[1] = number of files to analyze, if not specified takes all files in directory
		// args[2] = minimum polar angle of track to analyze, if not specified = 0 degrees
		// args[3] = maximum polar angle of track to analyze, if not specified = 90 degrees
		// args[4] = plots requested, 1 = sector vs superlayer, 2 = sector vs layer, else no plots

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

		// check for theta restrictions
		double theta_min;
		double theta_max;
		if (args.length < 3) {
			theta_min = 0.0;
			theta_max = 90.0; // if arguments not entered, assumes you want the whole range
		} else if (args.length == 3) {
			println("ERROR: You entered a minimum theta but not a maximum. Enter a 4th arg.");
			System.exit(1);
		} else if (Integer.parseInt(args[2])>Integer.parseInt(args[3])) { // check min < max
			println("ERROR: You entered a minimum theta that was greater than the maximum.");
			println("Please have arg[2]<arg[3].");
			System.exit(2);
		}
		else if (args.length >= 4) {
			theta_min = Integer.parseInt(args[2]);
			theta_max = Integer.parseInt(args[3]);
		}

		int plot_test = -1; // 1 = sector vs layer, 2 = sector vs superlayer, else = no plots
		if (args.length < 5) {
			println("Plots not requested."); 
			println("5th argument should be 1 (sector vs superlayer) or 2 (sector vs layer).");
			println();
		}
		if (args.length >= 5 ) {
			if (Integer.parseInt(args[4])==1) {
				println("Requested sector vs superlayer plots."); println();
				plot_test = 1;
			} else if (Integer.parseInt(args[4])==2) {
				println("Requested sector vs layer plots."); println();
				plot_test = 2;
			} else {
				println("Plots not requested."); 
				println("Please give an int (1 or 2) as the fifth argument."); println();
				System.exit(3);
			}
		}
		

		//Create JFrame and set to a nice size and aspect ratio
		JFrame[] frame = new JFrame[6];
		//Initialize EmbeddedCanvas
		EmbeddedCanvas[] sVsl_canvas = new EmbeddedCanvas[6];
		EmbeddedCanvas[] sVl_canvas = new EmbeddedCanvas[6];
		// declare histograms
		H1F[][] sVsl_array = new H1F[6][6]; // (i.e., sector vs superlayer array)
		H1F[][] sVl_array = new H1F[6][36]; // (i.e., sector vs layer array, will just scale 1-32)
		// declare fits for those histograms
		F1D[][] sVsl_fits = new F1D[6][6]; // 
		F1D[][] sVl_fits = new F1D[6][36]; 
		// declare number of bins and boundaries
		int n_bins = 150;
		double min_bin = -0.5*10000;
		double max_bin = 0.5*10000;
		int current_canvas = 0;
		// set up
		for (int sector; sector<6; sector++){
			frame[sector] = new JFrame("Spatial Residual Distributions, Hall B Run 002467");
			frame[sector].setSize(1200,900);
			sVsl_canvas[sector] = new EmbeddedCanvas();
			sVsl_canvas[sector].divide(2,3);
			sVl_canvas[sector] = new EmbeddedCanvas();	
			sVl_canvas[sector].divide(6,6);		
			for (int superlayer; superlayer<6; superlayer++){
				sVsl_array[sector][superlayer] = 
					new H1F("S"+Integer.toString(sector+1)+" SL"+
					Integer.toString(superlayer+1),n_bins,min_bin,max_bin); // micrometers
				sVsl_array[sector][superlayer].setTitleX("Spatial Residual (micrometers)");
				sVsl_array[sector][superlayer].setTitleY("Counts");
				sVsl_array[sector][superlayer].setOptStat(10); //Each IDSet has own OptStat
				sVsl_array[sector][superlayer].setLineWidth(2);
				sVsl_array[sector][superlayer].setLineColor(21);
				sVsl_array[sector][superlayer].setFillColor(30 + (current_canvas % 6) + 2);
				// sVsl_fits[sector][superlayer].setOptStat(1110);
				for (int layer; layer<6; layer++){
					sVl_array[sector][(superlayer)*6+(layer)] = 
					new H1F("S"+Integer.toString(sector+1)+" L"+
					Integer.toString((superlayer)*6+(layer)+1),n_bins,min_bin,max_bin); // mum
					sVl_array[sector][(superlayer)*6+(layer)].setTitleX(
						"Residuals (mum)");
					sVl_array[sector][(superlayer)*6+(layer)].setTitleY("Counts");
				}
			}
		current_canvas++;
		}

		for (int current_file; current_file<n_files; current_file++) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files);
			// limit to a certain number of files defined by n_files
			HipoDataSource reader = new HipoDataSource();
			reader.open(hipo_list[current_file]); // open next hipo file

			while(reader.hasEvent()==true){ // cycle through events
				HipoDataEvent event = reader.getNextEvent(); // load next event in the hipo file

				if (banks_test(event)) { // check that the necessary banks are present
					HipoDataBank hitBank= (HipoDataBank) event.getBank("TimeBasedTrkg::TBHits");
					HipoDataBank trkBank= (HipoDataBank) event.getBank("TimeBasedTrkg::TBTracks");
					HipoDataBank recBank= (HipoDataBank) event.getBank("REC::Particle");
					if (recBank.rows()>0) {
						
					for(int hitBankRow=0; hitBankRow<hitBank.rows(); hitBankRow++){
						// cycle over the entires in the hitbank (separate layers and regions)

						int trkID = hitBank.getInt("trkID", hitBankRow); // ID of reconstructed
						// track in the hitbank (many hits per track);

						float px = trkBank.getFloat("p0_x",trkID-1);
						float py = trkBank.getFloat("p0_y",trkID-1);
						float pz = trkBank.getFloat("p0_z",trkID-1);
						float theta = Math.toDegrees(theta_calculation(px,py,pz)); 

						boolean calculation_test = true;
						if (trkID==0) { // hit not assigned to track
							calculation_test = false;
						} else if (theta<theta_min) { // track theta outside requested range
							calculation_test = false;
						} else if (theta>theta_max) { // track theta outside requested range
							calculation_test = false;
						}

						if (calculation_test) {
							float doca = hitBank.getFloat("doca",hitBankRow);
							// distance of closest approach to wire
							float trkDoca = hitBank.getFloat("trkDoca",hitBankRow);
							// fitted distance of closest approach to wire
							int LR = hitBank.getInt("LR", hitBankRow);
							// trkDoca to the (L)eft or (R)ight of wire, determines sign of residual
							double residual = LR*(trkDoca-doca)*10000; //*10000 to convert to microm

							int sector = hitBank.getInt("sector", hitBankRow);
							int superlayer = hitBank.getInt("superlayer", hitBankRow);
							int layer = hitBank.getInt("layer", hitBankRow);

							sVsl_array[sector-1][superlayer-1].fill(residual);
							sVl_array[sector-1][(superlayer-1)*6+(layer)-1].fill(residual);
						}
					}
					}
				}
			}
		}

		// normalize all of the histograms created in the previous for loops
		for (int sector = 0; sector < 6; sector++) {
			for (int superlayer = 0; superlayer < 6; superlayer++) {
				sVsl_array[sector][superlayer].unit();
			}
		}
		for (int sector = 0; sector < 6; sector++) {
			for (int layer = 0; layer < 36; layer++) {
				sVl_array[sector][layer].unit();
			}
		}

		for (int sector=0; sector<6; sector++){
			println("For sector "+Integer.toString(sector+1)+" ...,");
			if (plot_test==2) {
				double[] means = new double[36];
				double[] sigmas = new double[36];
				double[] entries = new double[36];
				double[] meanError = new double[36];
				for (int layer = 0; layer<36; layer++){
					sVl_canvas[sector].cd(layer);
					sVl_canvas[sector].draw(sVl_array[sector][layer]);

					// calculate statistics about that sector
					means[layer] = sVl_array[sector][layer].getMean();
					sigmas[layer] = sVl_array[sector][layer].getRMS();
					entries[layer] = sVl_array[sector][layer].getEntries();
					meanError[layer] = sigmas[layer]/Math.pow(entries[layer],0.5);
					double roundedMean = round(means[layer],2);
					double roundedSigma = round( sigmas[layer],2);
					double roundedMeanError = round(meanError[layer],2);
					println("layer "+Integer.toString(layer+1)+", mean = "+roundedMean+" +/- "
						+roundedMeanError+" micrometers, sigma = "+roundedSigma+" micrometers");
				}
				double RMS = round(rootMeanSquare(means),2);
				double sigmaRMS = round(rootMeanSquare(sigmas),2);
				double errorMean = round(rootMeanSquare(meanError),2);
				println("The RMS of the means is "+RMS+" +/- "+errorMean+" (micrometers)");
				println();

				frame[sector].add(sVl_canvas[sector]);
				frame[sector].setLocationRelativeTo(null);
				frame[sector].setVisible(true);
			}
			if (plot_test==1) {
				double[] means = new double[6];
				double[] sigmas = new double[6];
				double[] entries = new double[6];
				double[] meanError = new double[6];
				for (int superlayer=0; superlayer<6; superlayer++){
					// cd to appropriate canvas (separate canvases for each sector)
					sVsl_canvas[sector].cd(superlayer);
					sVsl_canvas[sector].draw(sVsl_array[sector][superlayer]);

					// fit the histogram
					sVsl_fits[sector][superlayer] = fit_function(n_bins, min_bin, 
						max_bin, sVsl_array[sector][superlayer]); 
					sVsl_fits[sector][superlayer].show();
					sVsl_fits[sector][superlayer].setLineColor(2);
					sVsl_fits[sector][superlayer].setLineWidth(5); 
					sVsl_fits[sector][superlayer].setLineStyle(0);
					sVsl_canvas[sector].draw(sVsl_fits[sector][superlayer],"same");
					sVsl_fits[sector][superlayer].setOptStat(10001100);
					sVsl_canvas[sector].setAxisLabelSize(14);
					sVsl_canvas[sector].setStatBoxFontSize(14);
					sVsl_canvas[sector].setAxisTitleSize(18);
					println("Fit has chi2dof = "+
						(sVsl_fits[sector][superlayer].getChiSquare()/(n_bins-1)));

					// calculate statistics about that sector
					means[superlayer] = sVsl_array[sector][superlayer].getMean();
					sigmas[superlayer] = sVsl_array[sector][superlayer].getRMS();
					entries[superlayer] = sVsl_array[sector][superlayer].getEntries();
					meanError[superlayer] = sigmas[superlayer]/Math.pow(entries[superlayer],0.5);
					double roundedMean = round(means[superlayer],2);
					double roundedSigma = round( sigmas[superlayer],2);
					double roundedMeanError = round(meanError[superlayer],2);
					println("superlayer "+Integer.toString(superlayer+1)+", mean = "+roundedMean
						+" +/- "+roundedMeanError+" micrometers, sigma = "+roundedSigma
						+" micrometers");
				}
				double RMS = round(rootMeanSquare(means),2);
				double sigmaRMS = round(rootMeanSquare(sigmas),2);
				double errorMean = round(rootMeanSquare(meanError),2);
				println("The RMS of the means is "+RMS+" +/- "+errorMean+" (micrometers)");
				println();	

				frame[sector].add(sVsl_canvas[sector]);
				frame[sector].setLocationRelativeTo(null);
				frame[sector].setVisible(true);
			}
		}
	}
}