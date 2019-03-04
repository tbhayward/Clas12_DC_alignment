#include "Riostream.h"
#include <TROOT.h>
#include "TMinuit.h"
#include "TMath.h"
#include "TRandom1.h"
#include <ctime>

std::clock_t start = std::clock();
int s = (std::clock() - start) / (double) (CLOCKS_PER_SEC );
int sector = 3;


double calculate_residual(int sector, int i, double *par) {
	// nominal shifts by Raffaella de Vita were shifts of 0.2 cm and rotations of 0.1 degrees
	// may need to update these values in future
	double shift = 0.2;
	double rot = 0.1;

	// values in these lists are in micrometers
	// rotations around x and z axis are negligible when compared to other rotations so we do
	// not use them
	double r1_x_list[37] = {-797.4, -800.25, -804.39, -808.61, -813.13, -825.03, -722.85, -717.03, 
		-726.28, -726.51, -739.95, -745.82, 452.51, 455.97, 444.12, 435.27, 424.05, 427.05, 505.8, 
		491.82, 482.43, 482.86, 465.04, 463.05, -118.38, -149.58, -159.96, -164.55, -168.23, 
		-191.08, -115.89, -135.23, -134.44, -162.71, -160.44, -181.26, 6108.69};
	double r1_y_list[37] = {46.02, 69.28, 42.02, 37.28, 42.65, 45.12, -74.08, -69.22, -48.48, 
		-74.72, -46.65, -76.51, -80.23, -85.92, -89.29, -84.83, -98.63, -66.02, 61.85, 70.99, 
		67.09, 64.21, 55.46, 77.85, 7.47, 25.51, 45.91, 45.96, 66.05, 36.95, 10.49, -12.72, 8.87, 
		-35.2, -36.42, -16.63, -60.61};
	double r1_z_list[37] = {321.28, 323.74, 323.04, 326.64, 325.59, 336.33, 296.77, 297.37, 292.73, 
		301.43, 302.07, 301.71, -200.52, -187.34, -184.09, -186.44, -182.38, -170.9, -212.41, 
		-202.75, -200.24, -204.26, -200.15, -196.71, 62.55, 47.82, 65.91, 51.87, 82.21, 80.01, 
		50.1, 54.6, 69.74, 47.65, 77.35, 81.32, -2281.98};
	double r2_x_list[37] = {1572.31, 1573.63, 1560.14, 1550.23, 1544.01, 1548.73, 1418.25, 1417.22, 
		1406.21, 1394.59, 1395.18, 1399.21, -830.23, -834.32, -837.36, -847.54, -861.25, -865.09, 
		-933.13, -934.88, -942.31, -961.91, -966.95, -968.86, 350.88, 332.73, 330.29, 298.93, 
		303.33, 274.75, 319.14, 308.46, 309.54, 286.53, 297.46, 285.36, 7203.42};
	double r2_y_list[37] = {-88.34, -82.37, -82.0, -82.22, -84.11, -79.6, 96.55, 104.44, 95.06, 
		96.62, 95.85, 99.73, 133.89, 139.24, 145.4, 145.97, 139.02, 147.73, -140.5, -137.97, 
		-131.7, -134.62, -125.71, -121.08, -46.23, -52.31, -45.69, -45.86, -36.21, -52.43, 72.68, 
		52.72, 70.22, 51.58, 76.02, 63.74, 53.68};
	// double r2_z_list[37] = {-641.03, -641.54, -636.26, -631.25, -626.12, -624.27, -592.07, -586.99, 
	// 	-576.65, -580.65, -569.9, -569.46, 336.3, 355.66, 345.18, 354.42, 353.41, 363.81, 374.29, 
	// 	390.38, 389.47, 388.52, 397.77, 398.07, -140.4, -142.85, -122.69, -126.1, -101.52, -113.65,
	// 	-119.3, -125.71, -108.13, -129.05, -109.12, -114.54, -2626.56};
	double r2_z_list[37] = {521.971, 548.544, 527.642, 533.472, 582.077, 563.04, 502.122, 
		535.821, 519.498, 508.676, 501.179, 535.36, -269.489, -291.16, 
		-307.244, -273.598, -299.968, -312.822, -323.221, -319.807, -315.168, 
		-362.616, -327.05, -341.402, 105.699, 102.958, 85.6919, 89.0593, 
		103.78, 84.4495, 122.026, 145.611, 107.894, 89.4951, 92.6001, 68.4147};
	double r3_x_list[37] = {-795.18, -774.66, -768.05, -752.33, -740.37, -723.61, -707.76, -694.97, 
		-684.03, -678.54, -677.38, -660.0, 356.21, 381.33, 388.15, 401.35, 435.1, 451.7, 427.57, 
		439.04, 464.17, 476.2, 491.57, 511.14, -219.13, -201.17, -156.67, -157.1, -111.58, -92.44, 
		-208.66, -187.12, -161.97, -141.23, -109.25, -100.41, -7990.72};
	double r3_y_list[37] = {36.39, 40.91, 40.44, 38.8, 35.33, 39.83, -49.6, -49.08, -46.53, -46.41, 
		-48.76, -42.04, -75.21, -73.47, -63.99, -70.67, -74.37, -70.04, 56.96, 53.73, 53.22, 54.96, 
		55.41, 57.63, 39.78, 25.49, 33.28, 28.59, 34.4, 17.61, -21.17, -29.33, -22.26, -26.95, 
		-15.87, -12.43, -113.76};
	double r3_z_list[37] = {303.65, 299.99, 298.87, 292.43, 279.56, 283.19, 276.18, 283.66, 
		275.21, 271.31, 268.49, 268.7, -148.31, -149.01, -166.14, -167.38, -188.56, -185.04, 
		-177.38, -181.15, -184.95, -200.79, -203.71, -199.43, 88.64, 69.25, 77.99, 57.57, 50.29, 
		30.74, 93.8, 76.44, 78.94, 54.99, 45.9, 41.86, 2954.09};
	double r1_cy_list[37] = {38.9, 21.72, -9.59, -27.89, -43.74, -64.49, -1.33, -21.58, -40.09, 
		-65.38, -85.45, -105.26, -22.51, -9.31, -6.03, -5.47, -17.01, -2.62, 44.77, 54.96, 49.08, 
		44.79, 44.47, 49.55, 1.36, -7.28, 9.59, 6.89, 21.71, 5.38, -9.93, -22.23, -2.34, -20.87, 
		-11.66, -25.62, 237.82};
	double r2_cy_list[37] = {39.62, 46.81, 40.17, 27.57, 46.58, 41.33, 181.06, 186.56, 174.16, 
		173.99, 176.5, 187.57, 147.16, 123.0, 73.37, 44.83, 1.08, -7.08, -87.43, -121.46, -153.98, 
		-186.49, -227.14, -262.69, -25.45, -24.94, -21.42, -15.19, -3.57, -21.2, 78.03, 69.78, 
		94.22, 71.62, 90.84, 76.26, 687.95};
	double r3_cy_list[37] = {-73.44, -68.95, -69.77, -66.82, -67.4, -62.19, -134.94, -142.8, 
		-129.23, -125.95, -126.79, -129.27, -43.61, -32.3, -31.2, -35.08, -34.06, -27.42, 83.01, 
		99.52, 101.57, 85.7, 97.66, 103.64, 138.29, 75.28, 44.66, -33.25, -54.55, -98.94, 87.38, 
		11.59, -14.88, -63.26, -104.74, -151.1, -1215.61};

	// positions with may_2018_engineers alignment, calculated by Latif Kabir with coatjava 5c.6.9
	// in October 2018
	double sector_1_list[37] = {893.66, 930.71, 940.31, 923.36, 899.27, 849.3, 892.72, 878.03, 
		840.91, 817.05, 773.47, 782.62, -364.92, -374.89, -398.32, -445.76, -477.66, -508.96, 
		-499.85, -507.88, -535.52, -468.22, -539.05, -525.98, 176.02, 299.1, 119.0, 217.2, 80.61, 
		278.85, 114.55, 274.84, 156.54, 235.11, 85.16, 228.9, 36600.0};
	double sector_2_list[37] = {2031.36, 2223.43, 1923.65, 1748.79, 1570.44, 1550.43, 1433.84, 
		1410.67, 1418.26, 1394.3, 1304.7, 1246.63, -407.3, -405.97, -459.85, -549.73, -578.49, 
		-641.65, -823.21, -810.0, -932.61, -887.72, -940.09, -915.88, 277.22, 369.13, 140.73, 
		251.16, 46.59, 158.02, 346.76, 426.04, 270.24, 354.14, 192.71, 275.81, 39400.0};
	double sector_3_list[37] = {1578.47, 1572.34, 1626.18, 1604.2, 1560.11, 1533.06, 1453.68, 
		1412.44, 1435.55, 1404.79, 1333.29, 1331.6, -701.19, -745.18, -803.38, -843.44, -868.77, 
		-842.84, -874.92, -843.65, -938.23, -890.95, -967.45, -948.65, 304.11, 449.81, 237.04, 
		384.81, 151.14, 303.37, 292.05, 383.88, 247.22, 342.43, 190.97, 314.14, 31000.0};
	double sector_4_list[37] = {1094.01, 1133.17, 1143.98, 1133.1, 1125.19, 1111.25, 1014.19, 
		1041.32, 1003.56, 970.18, 943.05, 966.01, -574.85, -547.84, -552.46, -555.51, -568.25, 
		-549.68, -586.5, -591.15, -627.22, -594.7, -618.49, -613.13, 146.98, 293.15, 119.81, 
		286.71, 146.14, 278.04, 87.35, 237.83, 109.62, 254.31, 166.3, 291.0, 20733.33};
	double sector_5_list[37] = {-509.26, -492.97, -494.57, -512.46, -537.98, -544.52, -402.56, 
		-414.34, -440.96, -443.28, -483.73, -453.98, 270.6, 285.92, 210.89, 255.25, 235.05, 244.23, 
		279.33, 282.66, 176.46, 293.95, 219.29, 241.61, -109.59, 34.3, -110.23, 0.29, -156.26, 
		36.38, -142.88, 19.51, -94.9, 25.55, -69.62, 41.69, 36600.0};
	double sector_6_list[37] = {-116.17, -69.08, -100.33, -106.59, -100.07, -129.02, -38.6, -82.43, 
		-107.55, -143.37, -125.24, -122.67, 17.96, 47.13, 13.52, -73.23, 54.27, 63.98, 40.06, 
		41.09, 4.42, 48.13, 38.96, 63.77, -26.87, 90.12, -24.95, 100.49, -36.38, 93.53, -112.16, 
		88.37, -55.43, 53.51, -40.5, 141.24, 40333.33};

	double current_position;
	if (sector == 1) {
		current_position = sector_1_list[i];
	} else if (sector == 2)	{
		current_position = sector_2_list[i];
	} else if (sector == 3) {
		current_position = sector_3_list[i];
	} else if (sector == 4) {
		current_position = sector_4_list[i];
	} else if (sector == 5) {
		current_position = sector_5_list[i];
	} else if (sector == 6) {
		current_position = sector_6_list[i];
	}

	// desired vertex position is actually the position of the downstream target wall
	double desired_positions[37] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 30600.0};

	double value;
	if (i < 36) {
		// value = ((1/shift)*(par[0]*r1_x_list[i]+par[1]*r1_y_list[i]+par[2]*r1_z_list[i]+
		// par[3]*r2_x_list[i]+par[4]*r2_y_list[i]+par[5]*r2_z_list[i]+par[6]*r3_x_list[i]+
		// par[7]*r3_y_list[i]+par[8]*r3_z_list[i])+(1/rot)*(par[9]+r1_cy_list[i]+
		// par[10]+r2_cy_list[i]+par[11]+r3_cy_list[i]))-current_position;
		value = ((1/shift)*(par[0]*r1_x_list[i]+par[1]*r1_y_list[i]+par[2]*r1_z_list[i]+
		par[3]*r2_x_list[i]+par[4]*r2_y_list[i]+par[5]*r2_z_list[i]+par[6]*r3_x_list[i]+
		par[7]*r3_y_list[i]+par[8]*r3_z_list[i])+(1/rot)*(par[9]*r1_cy_list[i]+
		par[10]*r2_cy_list[i]+par[11]*r3_cy_list[i]))-current_position;
		// value = ((1/shift)*(par[0]*r1_x_list[i]+par[1]*r1_y_list[i]+par[2]*r1_z_list[i]+
		// par[3]*r2_x_list[i]+par[4]*r2_y_list[i]+par[5]*r2_z_list[i]+par[6]*r3_x_list[i]+
		// par[7]*r3_y_list[i]+par[8]*r3_z_list[i])+(cur))-current_position;
	} else { // rotations change vertex position by a negligible amount
		value = ((1/shift)*(par[0]*r1_x_list[i]+par[1]*r1_y_list[i]+par[2]*r1_z_list[i]+
		par[3]*r2_x_list[i]+par[4]*r2_y_list[i]+par[5]*r2_z_list[i]+par[6]*r3_x_list[i]+
		par[7]*r3_y_list[i]+par[8]*r3_z_list[i]))+current_position-desired_positions[i];
	}

	return value;
}

void chi2(int &npar, double *gin, double &f, double *par, int iflag) {
	//calculate chisquare
	double chisq = 0;
	// loop over all layers and vertex position (36 layers + 1 vertex) 
	sector = 6;
	for (int i = 0; i < 37; ++i) {
		// each value weighted equally, so no division by uncertainty, maybe adjust later?
		chisq += pow(calculate_residual(sector, i, par), 2);
	}
	f = chisq; 
}

void minimizer_0003() {

	for (int cur_sector = 1; cur_sector < 7; ++cur_sector) {
		TMinuit *gMinuit = new TMinuit(12);  //initialize TMinuit with a maximum of 12 params
		gMinuit->SetPrintLevel(0);
		gMinuit->SetFCN(chi2); // chi^2 function to be minimized
		double arglist[10]; arglist[0] = 1;
		int ierflg = 0;
		gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
		// Set starting values and step sizes for parameters
		// (param int, name, start value, step size, min limit, upper limit, error flag)
		gMinuit->mnparm(0, "r1x", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(1, "r1y", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(2, "r1z", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(3, "r1cy", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(4, "r2x", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(5, "r2y", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(6, "r2z", 0, 0.1, 0, 0, ierflg);
		gMinuit->mnparm(7, "r2cy", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(8, "r3x", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(9, "r3y", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(10, "r3z", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(11, "r3cy", 0, 0, 0, 0, ierflg);

		arglist[1] = 1.0; arglist[0] = 50000;
		gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

		sector++;

	}
	// TMinuit *gMinuit = new TMinuit(12);  //initialize TMinuit with a maximum of 12 params
	// gMinuit->SetPrintLevel(0);
	// gMinuit->SetFCN(chi2); // chi^2 function to be minimized
	// double arglist[10]; arglist[0] = 1;
	// int ierflg = 0;
	// gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
	// // Set starting values and step sizes for parameters
	// // (param int, name, start value, step size, min limit, upper limit, error flag)
	// gMinuit->mnparm(0, "r1x", 0, 0.001, 0, 0, ierflg);
	// gMinuit->mnparm(1, "r1y", 0, 0.001, 0, 0, ierflg);
	// gMinuit->mnparm(2, "r1z", 0, 0.001, 0, 0, ierflg);
	// gMinuit->mnparm(3, "r1cy", 0, 0.001, 0, 0, ierflg);
	// gMinuit->mnparm(4, "r2x", 0, 0.001, 0, 0, ierflg);
	// gMinuit->mnparm(5, "r2y", 0, 0.001, 0, 0, ierflg);
	// gMinuit->mnparm(6, "r2z", 0, 0.001, 0, 0, ierflg);
	// gMinuit->mnparm(7, "r2cy", 0, 0.001, 0, 0, ierflg);
	// gMinuit->mnparm(8, "r3x", 0, 0.001, 0, 0, ierflg);
	// gMinuit->mnparm(9, "r3y", 0, 0.001, 0, 0, ierflg);
	// gMinuit->mnparm(10, "r3z", 0, 0.001, 0, 0, ierflg);
	// gMinuit->mnparm(11, "r3cy", 0, 0.001, 0, 0, ierflg);

	// arglist[1] = 1.0; arglist[0] = 5000;
	// gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
}