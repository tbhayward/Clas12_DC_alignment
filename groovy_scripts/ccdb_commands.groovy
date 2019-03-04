public class ccdb_commands {

	public static void main(String[] args) {

		int alignment_run;
		if (args.length == 0) {
			// exits program if input directory not specified 
    	   	println("ERROR: Please enter a ccdb variation name as the first argument.");
    	  	System.exit(0);
    	} else if (args.length == 1) {
    		// exits program if alignment period not specified
    		println("ERROR: Please specify an alignment period as the second argument.");
    		println("Spring 2018 (enter 0), Fall 2018 (enter 1), Spring 2019 (2).");
    		println("(Fall 2019 not yet programmed.)");
    		System.exit(1);
    	} else {
    		alignment_run = Integer.parseInt(args[1]);
    	}

    	print("Specified ");
    	if (alignment_run==0) {
    		println("Spring 2018 alignment table.");
    	} else if (alignment_run==1) {
    		println("Fall 2018 alignment table.");
    	} else {
    		println("Spring 2019 alignment table.");
    	}
    	println(); println();

    	print("ccdb -c mysql://clas12writer:geom3try@clasdb/clas12 mkvar "+args[0]+";");

    	print(" ccdb -c mysql://clas12writer:geom3try@clasdb/clas12 -v "+args[0]);
    	print(" add /geometry/dc/region region.txt;");

		print(" ccdb -c mysql://clas12writer:geom3try@clasdb/clas12 -v "+args[0]);
    	print(" add /geometry/dc/superlayer superlayer.txt;");   

    	print(" ccdb -c mysql://clas12writer:geom3try@clasdb/clas12 -v "+args[0]);
    	print(" add /geometry/dc/alignment "+args[0]+".txt;"); 

    	print(" cp /work/clas12/thayward/dc_alignment/myClara_instructions/coatjava_5p7p4/");
    	if (alignment_run==0) {
    		print("thayward_test_0087.config");		
    	} else if (alignment_run==1) {
    		print("thayward_test_0085.config");
    	}
    	print(" /work/clas12/thayward/dc_alignment/myClara_instructions/coatjava_5p7p4/");
    	print(args[0]+".config;");

    	print(" cp /work/clas12/thayward/myClara/plugins/clas12/config/");
    	if (alignment_run==0) {
    		print("thayward_test_0087.yaml");		
    	} else if (alignment_run==1) {
    		print("thayward_test_0085.yaml");
    	}
    	print(" /work/clas12/thayward/myClara/plugins/clas12/config/"+args[0]+".yaml;");

    	println(); println();

    	print("gedit /work/clas12/thayward/dc_alignment/myClara_instructions/coatjava_5p7p4/");
    	println(args[0]+".config;");
    	println();
    	println("gedit /work/clas12/thayward/myClara/plugins/clas12/config/"+args[0]+".yaml;");
    	println(); println();

	}

}