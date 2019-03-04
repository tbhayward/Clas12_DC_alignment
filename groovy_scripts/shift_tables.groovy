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

public class shift_tables {

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
		// convert cartesian coordinates to polar angle
		double r = Math.pow(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2),0.5);
		return (double) Math.atan(y/x);
	}

	private static double theta_calculation (float x, float y, float z) {
		// convert cartesian coordinates to polar angle
		double r = Math.pow(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2),0.5);
		return (double) Math.acos(z/r);
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
		File[] hipo_list_nominal;
		if (args.length == 0) {
			// exits program if input directory not specified 
    	   	println("ERROR: Please enter a hipo file directory as the first argument");
    	  	System.exit(0);
    	} else {
    		File directory = new File(args[0]);
    		hipo_list_nominal = directory.listFiles();
    	}

    	File[] hipo_list_shift;
		if (args.length < 2) {
			// exits program if input directory not specified 
    	   	println("ERROR: Please enter a hipo file directory as the second argument");
    	  	System.exit(1);
    	} else {
    		File directory = new File(args[1]);
    		hipo_list_shift = directory.listFiles();
    	}

		int lower_theta = Integer.parseInt(args[2]);
		int upper_theta = Integer.parseInt(args[3]);

		int n_files;
		if (args.length < 5) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting number of files to one.");
			n_files = hipo_list_shift.size();
		} else{
			// if specified, convert to int
			n_files = Integer.parseInt(args[4]);
		}

		HipoDataSource reader = new HipoDataSource();
		int n_bins = 400;
		double min_bin = -0.8*10000;
		double max_bin = 0.8*10000;

		def nominal_vertex_cuts_lower = [5.462-17.940, 19.245-21.314, -14.711-15.480, 
			-18.594-22.112, -32.732-30.953, -14.989-24.797];

		def r1x_vertex_cuts_lower = [-7.243-19.445, 7.661-26.814, -30.153-29.022, 
			-15.668-24.304, -37.766-25.516, -7.075-1.679];
		def r2x_vertex_cuts_lower = [-11.675-26.098, 1.807-13.24, -51.134-42.341, 
			-15.756-27.500, -51.160-3.845, -5.290-22.835];
		def r3x_vertex_cuts_lower = [26.254-25.167, 43.082-26.496, 26.862-21.171, 
			-15.986-24.312, -21.865-21.734, -18.589-25.926];

		def r1y_vertex_cuts_lower = [4.564-24.024, 21.114-22.427, -13.731-18.119, 
			-15.560-25.417, -33.890-27.804, -14.080-22.010];
		def r2y_vertex_cuts_lower = [5.431-20.690, 21.618-22.450, -14.286-19.127, 
			-15.758-26.248, -36.063-3.0830, -13.434-21.929];
		def r3y_vertex_cuts_lower = [6.009-22.820, 19.291-25.523, -13.456-15.935, 
			-16.381-25.817, -32.145-25.533, -10.717-24.484];

		def r1z_vertex_cuts_lower = [7.948-27.332, 24.320-21.917, -15.052-20.319, 
			-14.343-21.917, -32.849-32.897, -10.197-23.551];
		def r2z_vertex_cuts_lower = [9.350-25.257, 22.912-24.082, -11.026-16.619, 
			-15.240-18.504, -30.039-25.004, -8.0721-23.243];
		def r3z_vertex_cuts_lower = [19.000-21.124, 18.535-25.112, -17.500-14.826, 
			-20.718-19.328, -38.927-29.106, -12.438-23.471];

		def r1cy_vertex_cuts_lower = [11.782-17.433, 15.625-22.742, -18.844-21.661, 
			-21.450-25.883, -38.008-28.690, -16.764-25.982];
		def r2cy_vertex_cuts_lower = [5.431-20.690, 21.618-22.450, -14.286-19.127, 
			-15.758-26.248, -36.063-3.0830, -13.434-21.929];
		def r3cy_vertex_cuts_lower = [6.009-22.820, 19.291-25.523, -13.456-15.935, 
			-16.381-25.817, -32.145-25.533, -10.717-24.484];

		def current_shift = nominal_vertex_cuts_lower;


		H1F[][] nominal_residuals = new H1F[6][36]; // sector, layer
		H1F[][] shift_residuals = new H1F[6][36]; // sector, layer
		for (int sector = 0; sector < 6; sector ++) {
			for (int layer = 0; layer < 36; layer++) {
				nominal_residuals[sector][layer] = new H1F("S"+Integer.toString(sector+1)+" L"+
					Integer.toString(layer+1),n_bins,min_bin,max_bin);
				shift_residuals[sector][layer] = new H1F("S"+Integer.toString(sector+1)+" L"+
					Integer.toString(layer+1),n_bins,min_bin,max_bin); 
			}
		}

		for (int current_file; current_file<n_files; current_file++) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files+" of the nominal geometries.");
			// limit to a certain number of files defined by n_files
			reader.open(hipo_list_nominal[current_file]); // open next hipo file

			while(reader.hasEvent()==true){ // cycle through events
				HipoDataEvent event = reader.getNextEvent(); // load next event in the hipo file
				if (banks_test(event)) { // check that the necessary banks are present
					HipoDataBank hitBank= (HipoDataBank) event.getBank("TimeBasedTrkg::TBHits");
					HipoDataBank trkBank= (HipoDataBank) event.getBank("TimeBasedTrkg::TBTracks");
					HipoDataBank trajBank= (HipoDataBank) event.getBank("TimeBasedTrkg::Trajectory");
					HipoDataBank recBank= (HipoDataBank) event.getBank("REC::Particle");
					HipoDataBank ccBank= (HipoDataBank) event.getBank("REC::Cherenkov");
					HipoDataBank calBank= (HipoDataBank) event.getBank("REC::Calorimeter");

					int sector;
					float vtx0_z = 10*trkBank.getFloat("Vtx0_z", 0);

					if (track_check(recBank, ccBank, calBank)&&(trkBank.rows()==1)) {
						sector = trkBank.getInt("sector", 0);
						float px = trkBank.getFloat("p0_x",0);
						float py = trkBank.getFloat("p0_y",0);
						float pz = trkBank.getFloat("p0_z",0);
						// float phi = Math.toDegrees(phi_calculation(px,py,pz));
						double theta = Math.toDegrees(theta_calculation(px,py,pz));
						// if (vtx0_z > current_shift[sector-1]) {
						// if (keep_track(theta)) {
						if (theta>lower_theta&&theta<upper_theta) {
							for(int hitBankRow=0; hitBankRow<hitBank.rows(); hitBankRow++){
								if (hitBank.getInt("trkID", hitBankRow) == 1) {
									float residual = 
										10000*hitBank.getFloat("fitResidual", hitBankRow);
									int superlayer = hitBank.getInt("superlayer", hitBankRow);
									int layer = hitBank.getInt("layer", hitBankRow);

									nominal_residuals[sector-1][(superlayer-1)*6+(layer)-1].fill(residual)
								}
							}
						}
					}
				}
			}
		}

		int num_sec1 = 0;
		int num_sec2 = 0;
		int num_sec3 = 0;
		int num_sec4 = 0;
		int num_sec5 = 0;
		int num_sec6 = 0;

		for (int current_file; current_file<n_files; current_file++) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files+" of the shifted geometries.");
			// limit to a certain number of files defined by n_files
			reader.open(hipo_list_shift[current_file]); // open next hipo file

			while(reader.hasEvent()==true){ // cycle through events
				HipoDataEvent event = reader.getNextEvent(); // load next event in the hipo file
				if (banks_test(event)) { // check that the necessary banks are present
					HipoDataBank hitBank= (HipoDataBank) event.getBank("TimeBasedTrkg::TBHits");
					HipoDataBank trkBank= (HipoDataBank) event.getBank("TimeBasedTrkg::TBTracks");
					HipoDataBank trajBank= (HipoDataBank) event.getBank("TimeBasedTrkg::Trajectory");
					HipoDataBank recBank= (HipoDataBank) event.getBank("REC::Particle");
					HipoDataBank ccBank= (HipoDataBank) event.getBank("REC::Cherenkov");
					HipoDataBank calBank= (HipoDataBank) event.getBank("REC::Calorimeter");

					int sector;
					float vtx0_z = 10*trkBank.getFloat("Vtx0_z", 0);

					if (track_check(recBank, ccBank, calBank)&&(trkBank.rows()==1)) {
						sector = trkBank.getInt("sector", 0);
						float px = trkBank.getFloat("p0_x",0);
						float py = trkBank.getFloat("p0_y",0);
						float pz = trkBank.getFloat("p0_z",0);
						// float phi = Math.toDegrees(phi_calculation(px,py,pz));
						double theta = Math.toDegrees(theta_calculation(px,py,pz));
						// if (vtx0_z > current_shift[sector-1]) {
						// if (keep_track(theta)) {
						if (theta>lower_theta&&theta<upper_theta) {

							switch(sector) {
							case 1:
								num_sec1++;
								break
							case 2:
								num_sec2++;
								break
							case 3:
								num_sec3++;
								break
							case 4:
								num_sec4++;
								break
							case 5:
								num_sec5++;
								break
							case 6:
								num_sec6++;
								break	
							}

							for(int hitBankRow=0; hitBankRow<hitBank.rows(); hitBankRow++){
								if (hitBank.getInt("trkID", hitBankRow) == 1) {
									float residual = 
										10000*hitBank.getFloat("fitResidual", hitBankRow);
									int superlayer = hitBank.getInt("superlayer", hitBankRow);
									int layer = hitBank.getInt("layer", hitBankRow);

									shift_residuals[sector-1][(superlayer-1)*6+(layer)-1].fill(residual)
								}
							}
						}
					}
				}
			}
		}

		double[][][] values = new double[2][6][36];

		println(num_sec1);
		println(num_sec2);
		println(num_sec3);
		println(num_sec4);
		println(num_sec5);
		println(num_sec6);

		for (int sector = 0; sector < 6; sector++) {
			GraphErrors nominal_residuals_points = new GraphErrors();
			GraphErrors shift_residuals_points = new GraphErrors();
			GraphErrors difference_points = new GraphErrors();
			double[][] differences = new double[6][36];

			println(); println("Sector "+(sector+1));
			print("{");
			for (int layer = 0; layer < 36; layer++) {
				double shift_mean = shift_residuals[sector][layer].getMean();
				double shift_std = shift_residuals[sector][layer].getRMS();
				double shift_counts = shift_residuals[sector][layer].getEntries();
				double shift_error = shift_std/Math.sqrt(shift_counts);
				double nominal_mean = nominal_residuals[sector][layer].getMean();
				double nominal_std = nominal_residuals[sector][layer].getRMS();
				double nominal_counts = nominal_residuals[sector][layer].getEntries();
				double nominal_error = nominal_std/Math.sqrt(nominal_counts);
				double error = Math.sqrt(shift_error*shift_error+nominal_error*nominal_error);
				values[0][sector][layer] = nominal_mean-shift_mean;
				values[1][sector][layer] = nominal_std;

				if (layer<35) {
					print((nominal_mean-shift_mean).round(1)+", ");	
				} else {
					print((nominal_mean-shift_mean).round(1));
				}
				// if (layer<35) {
				// 	print((error).round(1)+", ");	
				// } else {
				// 	print((error).round(1));
				// }
			}
			print("}"); println(); println();
		}

		println(); print("{")
		for (int layer = 0; layer < 36; layer++) {
			def mean_values = [values[0][0][layer], values[0][1][layer], 
				values[0][2][layer], values[0][3][layer], values[0][4][layer], values[0][5][layer]];
			def std_values = [values[1][0][layer], values[1][1][layer], 
				values[1][2][layer], values[1][3][layer], values[1][4][layer], values[1][5][layer]];

			if (layer<35) {
				print((mean_values.sum()/mean_values.size()).round(1)+", ");
			} else	{
				print((mean_values.sum()/mean_values.size()).round(1));
			}
		}
		print("}"); println();

		println(); print("{")
		for (int layer = 0; layer < 36; layer++) {
			def std_values = [values[1][0][layer], values[1][1][layer], 
				values[1][2][layer], values[1][3][layer], values[1][4][layer], values[1][5][layer]];

			if (layer<35) {
				print((std_values.sum()/std_values.size()).round(3)+", ");
			} else	{
				print((std_values.sum()/std_values.size()).round(3));
			}
		}
		print("}"); println(); println();

	}
		
}