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

public class chi2_results {

	public static boolean banks_test(HipoDataEvent event){
		boolean banks_result = true; // check to see if the event has all of the banks present
		// TBHits required for residuals, TBTracks required for angular limit
		if (!(event.hasBank("TimeBasedTrkg::TBHits"))) {
        	banks_result = false;
    	} else if (!(event.hasBank("TimeBasedTrkg::TBTracks"))) {
    		banks_result = false;
    	} else if (!(event.hasBank("REC::Particle"))) {
    		banks_result = false;
    	}
    	return banks_result;
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

		double chi2 = 0;
		int n_hits = 0;
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
					if (recBank.rows()==1) {
						for(int hitBankRow=0; hitBankRow<hitBank.rows(); hitBankRow++){
							// cycle over the entires in the hitbank (separate layers and regions)

							int trkID = hitBank.getInt("trkID", hitBankRow); // ID of reconstructed
							// track in the hitbank (many hits per track); 

							int superlayer = hitBank.getInt("superlayer", hitBankRow);

							boolean calculation_test = true;
							if (trkID==0) { // hit not assigned to track
								calculation_test = false;
							} else if (superlayer<5) {
								calculation_test = false;
							}

							if (calculation_test) {

								float doca = hitBank.getFloat("doca",hitBankRow);
								// distance of closest approach to wire

								float docaError = hitBank.getFloat("docaError",hitBankRow);
								// uncertainty on doca of the hit calculated from TDC (in cm)

								float trkDoca = hitBank.getFloat("trkDoca",hitBankRow);
								// fitted distance of closest approach to wire

								float timeResidual = hitBank.getFloat("timeResidual",hitBankRow);
								// time residual of the hit (in cm)

								// println( Math.pow(trkDoca-doca,2) / 
								// 	(Math.pow(docaError,2) + Math.pow(timeResidual,2)) )

								chi2+= Math.pow(trkDoca-doca,2) / 
									(Math.pow(docaError,2) + Math.pow(timeResidual,2)); 
								n_hits++;
							}
						}
					}
				}
			}
		}
		println("chi^2 = "+chi2);
		println(chi2/n_hits);
    }
}