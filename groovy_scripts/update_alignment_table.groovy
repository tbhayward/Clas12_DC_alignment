/*
 * author Timothy B. Hayward
 * April 2018
 * CLAS (12) Collaboration Service Work 
 * (with Torri Roark and Tyler Viducic, supervised by Mac Mestayer)
 */

public class update_alignment_table {

	public static double ccdb(String filename) {
		// ccdb mysql command for updating alignment data table
		// note the spaces at the end of the strings
		// user must have ccdb write permission
		String ccdb_command = "ccdb -c mysql://clas12writer:geom3try@clasdb/clas12 add ";
		String table_location = "geometry/dc/alignment ";
		String alignment_table = filename;
		String change_table_command = ccdb_command+table_location+alignment_table;
		change_table_command.execute(); // executes command line input

		return 0;
	}

	public static double write_file(File file, int region, int sector, int var, double value) {
		// write the alignment table data file 
     	def lines = file.readLines(); // read in all the lines (rows) in the file
     	String new_table = ""; // ongoing string of new table to be written to the file
     	for (int current_line=0; current_line<lines.size(); current_line++) { // iterate over lines
     		String[] s = lines[current_line].split(" "); // splits line into array without spaces
     		if (current_line == ( (region-1)*6+(sector)-1 ) ) { 
     			// if current line is the input region and sector to be modified first 3 entries
     			// are the same as before
     			new_table+=(s[0]+" "+s[1]+" "+s[2]+" ");
     			for (int current_var=0; current_var<6; current_var++) {
     				// write values into we arrive at the variable to be edited
     				if (current_var==var-1) {
     					new_table+=(Double.toString(value)+" ");
     				} else {
     					new_table+=(s[3+current_var]+" ");
     				}
     			}
     		} else {
     			// otherwise the new line is the same as the old line
     			for (int i=0; i<s.size(); i++){
     			new_table+=(s[i]+" ");
     			}   
     		}
     		new_table+=("\n"); // add new line to string
     	}
     	file.write(new_table);

		return 0;
	}

	public static double clear_file(File file) {
		// clear the alignment table data file 
		// similar to write_file except we place zeroes everywhere
     	def lines = file.readLines();
     	String new_table = "";
     	for (int current_line=0; current_line<lines.size(); current_line++) {
     		String[] s = lines[current_line].split(" ");  
     		for (int i=0; i<3; i++) {
     			new_table+=(s[i]+" ");
     		}
     		for (int i=3; i<s.size(); i++) {
     			new_table+=("0 ");
     		}   
     		new_table+=("\n");
     	}
     	file.write(new_table);
     	
		return 0;
	}

	public static void main(String[] args) {
		println("Beginning update_alignment_table."); println();
		// Accepts input arguments: data_table_location.txt, region, sector and variable
		// to modify and the value to place in that location. args[1-5].
		// If only data_table_location specified, clears data table of all shifts.

		if (args.length == 0) {
			// exits program if input alignment_table not specified 
        	println("ERROR: Enter a text file for the alignment table.");
       		System.exit(0);

    	} else if (args.length == 1) {
    		// clears data table (sets all values to 0) if only input file given

    		// check that first argument is a file
    		String filename = args[0]; // input txt file
    		File file = new File(filename);
    		if (!file.isFile()) {
    			println("ERROR: Please enter a text file as the first argument.");
    			System.exit(1);
    		}

    		// create backup before writing over file incase user was unaware of input arguments
    		println("WARNING: No input region, sector, variable or value given.")
    		println("Creating back up file.");
    		String command = "cp "+filename+" "+filename+"_backup";  
    		command.execute();

    		println("Clearing geometry/dc/alignment and resetting all values to 0.");
    		clear_file(file);
    		println("Uploading reset data table to ccdb.");
    		ccdb(filename); // update the ccdb table with args[0] file

    	} else {
    		if (args.length < 5 ) {
    			// exits program if region, sector, variable and value not given
    			println("ERROR: Please enter a region, sector and variable to edit and a value.");
    			System.exit(3);
    		} else if (!args[1].isNumber()) {
    			println("ERROR: Please enter an int [1-3] for the region as the 2nd argument.");
    			System.exit(4);
    		} else if (!args[2].isNumber()) {
    			println("ERROR: Please enter an int [1-6] for the sector as the 3rd argument.");
    			System.exit(5);
    		} else if (!args[3].isNumber()) {
    			println("ERROR: Please enter an int [1-6] for the variable as the 4th argument.");
    			System.exit(6);
    		} else if (!args[4].isNumber()) {
    			println("ERROR: Please enter a double for the input variable as the 5th argument.");
    			System.exit(7);
    		}

   			int region_edit = Integer.parseInt(args[1]); // region number to be edited
    		int sector_edit = Integer.parseInt(args[2]); // sector number to be edited
    		int var_edit = Integer.parseInt(args[3]); // variable number to be edited
    											  	  // [dx dy dz dtheta_x dtheta_y dtheta_z]
    		double value = Double.parseDouble(args[4]); // value to be set

    		String filename = args[0]; // input txt file
    		File file = new File(filename);

    		write_file(file, region_edit, sector_edit, var_edit, value);
			ccdb(filename); // update the ccdb table with args[0] file
    	} 
    	println(); println("Ending update_alignment_table."); println();
	}
}	