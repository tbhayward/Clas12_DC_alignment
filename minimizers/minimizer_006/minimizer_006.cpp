#include "Riostream.h"
#include <TROOT.h>
#include "TMinuit.h"
#include "TMath.h"
#include "TRandom1.h"
#include <ctime>

std::clock_t start = std::clock();
int s = (std::clock() - start) / (double) (CLOCKS_PER_SEC );
// int sector = 6;

double r1x[6];
double r2x[6];
double r3x[6];
double r1y[6];
double r2y[6];
double r3y[6];
double r1xGlobal[6];
double r2xGlobal[6];
double r3xGlobal[6];
double r1yGlobal[6];
double r2yGlobal[6];
double r3yGlobal[6];
double r1z[6];
double r2z[6];
double r3z[6];
double r1cy[6];
double r2cy[6];
double r3cy[6];
double error;
int cur_sector = 0;


double calculate_residual(int sector, int i, double *par) {
	// nominal shifts by Raffaella de Vita were shifts of 0.2 cm and rotations of 0.1 degrees
	// may need to update these values in future
	double shift = 0.2;
	double rot = 0.2;

	// values in these lists are in micrometers
	// rotations around x and z axis are negligible when compared to other rotations so we do
	// not use them
	
	double r1_x_list[37] = {-666.054, -668.284, -691.604, -677.457, -741.557, -719.981, -629.643, 
		-580.462, -624.003, -621.22, -611.808, -639.731, 363.355, 389.918, 332.532, 369.24, 
		359.263, 334.416, 376.445, 393.833, 390.876, 411.866, 426.337, 376.316, -87.177, 
		-97.426, -121.381, -128.954, -122.852, -160.674, -100.887, -138.387, -154.237, -147.846, 
		-110.104, -179.783, 5.6365};

	double r2_x_list[37] = {1298.442, 1272.766, 1207.232, 1261.492, 1289.027, 1213.194, 1098.59, 
		1129.902, 1151.636, 1170.459, 1179.08, 1098.308, -596.412, -594.838, -647.877, -682.722, 
		-589.469, -685.217, -762.512, -713.423, -792.403, -748.739, -707.787, -721.601, 330.93, 
		255.093, 276.468, 225.857, 232.001, 163.374, 265.465, 199.463, 232.22, 218.428, 
		209.436, 231.782, 12.8148};

	double r3_x_list[37] = {-603.553, -588.693, -595.11, -577.278, -563.347, -561.435, -546.509, 
		-547.068, -535.777, -524.431, -519.895, -518.032, 235.977, 235.96, 263.136, 276.761, 
		293.04, 304.513, 306.244, 317.992, 322.034, 346.032, 376.467, 375.014, -153.686, -124.853, 
		-111.142, -90.178, -72.623, -55.493, -160.009, -152.267, -129.482, -104.302, -83.115, 
		-52.701, -16.013};

	double r1_y_list[37] = {-79.364, -73.756, -74.636, -83.299, -72.356, -66.11, 95.548, 84.328, 
		87.778, 91.426, 85.311, 93.833, 135.17, 133.025, 131.193, 135.31, 127.734, 135.546, 
		-100.726, -107.206, -117.767, -105.739, -90.557, -99.443, -51.872, -38.514, -44.299, 
		-46.808, -41.729, -54.78, 49.148, 30.941, 42.682, 38.276, 37.383, 56.258, -4.03267};

	double r2_y_list[37] = {-83.883, -78.282, -78.71, -88.947, -77.632, -71.39, 89.517, 80.608, 
		82.791, 87.597, 82.503, 90.558, 139.76, 136.713, 131.231, 138.203, 129.013, 136.794, 
		-98.681, -105.06, -114.721, -105.428, -87.054, -98.145, -52.11, -38.765, -44.531, 
		-47.415, -43.203, -57.073, 47.87, 29.193, 42.826, 39.867, 37.652, 54.047, -0.637833};

	double r3_y_list[37] = {19.129, 23.713, 24.664, 10.159, 22.318, 25.914, -45.805, -58.827, 
		-54.549, -48.702, -53.888, -44.001, -40.424, -48.446, -49.387, -45.631, -55.291, 
		-49.54, 56.446, 51.536, 44.053, 62.19, 78.208, 65.492, 21.374, 32.753, 25.116, 
		19.451, 20.921, 2.221, -34.768, -49.566, -35.799, -34.587, -34.223, -10.311, -1.48667};

	double r1_z_list[37] = {183.827, 188.234, 172.546, 172.978, 154.727, 184.707, 155.292, 158.647, 
		171.765, 159.011, 158.367, 161.817, -94.677, -98.421, -80.528, -97.636, -85.628, -78.856, 
		-99.074, -115.746, -101.965, -104.323, -85.398, -91.878, 28.02, 35.573, 43.77, 40.453, 
		37.977, 29.914, 17.144, 18.358, 1.05, 22.932, 22.329, 51.827, -2.691};

	double r2_z_list[37] = {-354.718, -353.724, -352.946, -362.732, -355.786, -345.235, -305.438, 
		-317.293, -319.224, -312.611, -320.29, -303.305, 179.677, 180.733, 181.232, 187.768, 
		184.292, 190.501, 201.387, 199.808, 195.549, 211.934, 228.357, 217.279, -78.077, -61.318, 
		-66.647, -67.756, -56.858, -67.737, -75.76, -90.476, -76.712, -75.592, -72.123, 
		-48.574, -4.03398};

	double r3_z_list[37] = {163.062, 161.205, 160.659, 149.337, 158.006, 160.107, 154.414, 141.671, 
		144.801, 147.094, 139.203, 147.887, -59.75, -68.09, -73.733, -74.148, -82.476, -82.328, 
		-82.604, -93.882, -110.102, -95.312, -84.108, -100.666, 40.399, 46.627, 33.274, 24.734, 
		23.096, 1.374, 44.401, 24.351, 26.155, 17.997, 10.587, 26.148, -0.711833};

	double r1_cy_list[37] = {-87.031, -124.415, -157.712, -205.689, -235.521, -269.584, -131.727, 
		-189.647, -226.276, -267.499, -319.63, -346.232, 60.209, 53.082, 50.901, 54.964, 46.786, 
		53.662, 175.429, 164.819, 155.608, 168.938, 182.275, 164.461, -19.759, -6.447, -13.211, 
		-17.926, -17.976, -33.616, -48.413, -71.951, -65.264, -70.171, -75.982, -60.099, 1.89};

	double r2_cy_list[37] = {531.918, 516.841, 527.076, 500.331, 506.132, 527.554, 741.279, 
		726.974, 731.657, 729.757, 728.375, 733.229, 107.633, 42.134, -32.211, -86.344, 
		-162.457, -226.866, -352.227, -424.155, -509.825, -559.331, -613.26, -697.731, 
		37.667, 49.448, 35.735, 26.107, 23.304, 1.04, 211.918, 192.75, 194.082, 181.161, 
		175.435, 185.414, 9.414};

	double r3_cy_list[37] = {-500.35, -485.149, -501.976, -479.539, -497.403, -477.212, -609.019, 
		-591.206, -579.985, -589.29, -586.815, -562.37, 74.325, 81.648, 102.998, 89.516, 108.82, 
		120.869, 343.392, 343.445, 364.609, 375.536, 408.347, 413.507, 178.875, 95.877, 16.603, 
		-72.542, -158.425, -245.926, 63.804, -13.941, -115.304, -173.488, -253.459, 
		-310.563, -18.9125};

	double std_list[37] = {1388.342, 1377.905, 1364.688, 1357.442, 1367.341, 1371.056, 1354.424, 
		1340.077, 1340.371, 1326.919, 1331.286, 1324.379, 1256.927, 1258.792, 1278.055, 1291.668, 
		1265.617, 1248.795, 1278.81, 1280.334, 1305.279, 1305.465, 1283.226, 1291.915, 1179.257, 
		1191.516, 1216.565, 1213.86, 1226.134, 1225.673, 1225.745, 1245.754, 1254.115, 1258.633, 
		1281.931, 1258.904};
	
	double sector_1_list[37] = {892.407, 907.079, 930.198, 933.175, 917.42, 925.331, 859.265, 
		888.404, 863.707, 866.87, 852.505, 887.079, -359.274, -380.084, -421.774, -440.439, 
		-445.413, -474.137, -548.325, -543.204, -611.047, -513.617, -546.005, -542.592, 165.454, 
		299.502, 141.983, 246.497, 110.85, 201.739, 150.19, 256.275, 165.105, 264.215, 131.152, 
		231.474, 5.462-17.940};
	double sector_2_list[37] = {1904.076, 1886.408, 1872.771, 1792.57, 1716.684, 1766.655, 1573.768, 
		1556.347, 1565.869, 1566.679, 1541.778, 1541.415, -782.514, -766.777, -834.336, -819.533, 
		-872.652, -897.837, -940.321, -953.548, -1014.725, -958.553, -1017.066, -1012.84, 387.862, 
		464.518, 283.797, 385.931, 188.846, 267.987, 340.005, 389.75, 276.475, 396.409, 248.511, 
		347.371, 19.245-21.314};
	double sector_3_list[37] = {285.966, 306.153, 320.823, 348.476, 328.539, 374.66, 249.213, 
		278.841, 262.219, 290.181, 289.423, 327.344, -189.807, -179.265, -208.917, -185.457, 
		-240.264, -182.275, -230.254, -179.13, -217.505, -115.075, -155.316, -135.538, 80.249, 
		210.82, 18.46, 128.456, -7.201, 62.618, 15.998, 78.503, -1.851, 126.251, 25.819, 
		123.921, -14.711-15.480};
	double sector_4_list[37] = {912.768, 926.752, 950.053, 939.979, 964.467, 965.637, 831.808, 
		835.671, 827.861, 817.957, 851.772, 884.039, -483.228, -438.102, -439.901, -424.814, 
		-441.528, -424.817, -518.099, -497.973, -579.429, -443.538, -505.379, -506.364, 
		144.713, 255.697, 135.285, 257.738, 111.04, 218.652, 86.971, 175.935, 88.575, 
		241.813, 133.006, 245.478, -18.594-22.112};
	double sector_5_list[37] = {-951.024, -946.35, -931.842, -930.164, -961.684, -927.236, 
		-842.94, -844.687, -801.039, -802.099, -792.358, -796.042, 439.024, 477.269, 381.585, 
		486.711, 462.863, 505.646, 502.816, 510.173, 454.758, 572.991, 530.852, 504.935, 
		-206.894, -111.648, -221.812, -82.081, -205.836, -118.778, -245.905, -156.089, 
		-195.285, -101.743, -94.173, -50.982, -32.732-30.953};
	double sector_6_list[37] = {-290.199, -255.84, -264.804, -280.676, -270.515, -257.735, 
		-239.106, -217.429, -253.559, -264.269, -252.768, -238.677, 136.857, 129.166, 126.111, 
		103.42, 129.879, 143.063, 154.275, 138.32, 88.369, 188.018, 134.311, 165.137, -54.052, 
		62.589, -58.306, 20.077, -101.737, 1.86, -126.966, 7.046, -137.207, 48.392, -48.497, 
		98.293, -14.989-24.797};

	// double mag_A = 0;
	// double mag_B = 0;
	// double product = 0;
	// for (int i = 0; i < 36; ++i)
	// {
	// 	double value_A = r3_x_list_1[i];
	// 	double value_B = r3_cy_list_1[i];
	// 	product+=value_A*value_B;
	// 	mag_A+=value_A*value_A;
	// 	mag_B+=value_B*value_B;
	// }
	// double dot_product = product / (TMath::Power(mag_A,0.5)*TMath::Power(mag_B,0.5));
	// cout << dot_product << endl;

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
	// double desired_positions[38] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.72891, 5.69459};
	double desired_positions[37] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -5};

	double value = ((1/shift)*(par[0]*r1_x_list[i]+par[1]*r1_y_list[i]+par[2]*r1_z_list[i]+
		par[3]*r2_x_list[i]+par[4]*r2_y_list[i]+par[5]*r2_z_list[i]+par[6]*r3_x_list[i]+
		par[7]*r3_y_list[i]+par[8]*r3_z_list[i])+(1/rot)*(par[9]*r1_cy_list[i]+
		par[10]*r2_cy_list[i]+par[11]*r3_cy_list[i]))-current_position+desired_positions[i];
	value = value*value;
	// if (i<35) {
	// 	value = value / (std_list[i]*std_list[i]);
	// } else	{
	// 	value = value / (20*20);
	// }

	return value;
}

void chi2(int &npar, double *gin, double &f, double *par, int iflag) {
	//calculate chisquare
	double chisq = 0;
	// loop over all layers and vertex position (36 layers + 1 vertex) 
	// sector = 1;
	for (int i = 0; i < 36; ++i) {
		// each value weighted equally, so no division by uncertainty, maybe adjust later?
		chisq += pow(calculate_residual(cur_sector, i, par), 2);
		// cout << i << endl;
	}
	f = chisq; 
}

void minimizer_006() {

	for (int sector = 1; sector < 7; ++sector) {
		cur_sector = sector;
		TMinuit *gMinuit = new TMinuit(12);  //initialize TMinuit with a maximum of 12 params
		gMinuit->SetPrintLevel(-1);
		gMinuit->SetFCN(chi2); // chi^2 function to be minimized
		double arglist[10]; arglist[0] = 1;
		int ierflg = 0;
		gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
		// Set starting values and step sizes for parameters
		// (param int, name, start value, step size, min limit, upper limit, error flag)
		gMinuit->mnparm(0, "r1x", 0, 0.1, 0, 0, ierflg);
		gMinuit->mnparm(1, "r1y", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(2, "r1z", 0, 0.1, 0, 0, ierflg);
		gMinuit->mnparm(3, "r2x", 0, 0.1, 0, 0, ierflg);
		gMinuit->mnparm(4, "r2y", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(5, "r2z", 0, 0.1, 0, 0, ierflg);
		gMinuit->mnparm(6, "r3x", 0, 0.1, 0, 0, ierflg);
		gMinuit->mnparm(7, "r3y", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(8, "r3z", 0, 0.1, 0, 0, ierflg);
		gMinuit->mnparm(9, "r1cy", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(10, "r2cy", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(11, "r3cy", 0, 0, 0, 0, ierflg);

		arglist[1] = 1.0; arglist[0] = 50000;
		gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

		gMinuit->GetParameter(0,r1x[sector-1],error);
		gMinuit->GetParameter(1,r1y[sector-1],error);
		gMinuit->GetParameter(2,r1z[sector-1],error);
		gMinuit->GetParameter(3,r2x[sector-1],error);
		gMinuit->GetParameter(4,r2y[sector-1],error);
		gMinuit->GetParameter(5,r2z[sector-1],error);
		gMinuit->GetParameter(6,r3x[sector-1],error);
		gMinuit->GetParameter(7,r3y[sector-1],error);
		gMinuit->GetParameter(8,r3z[sector-1],error);
		gMinuit->GetParameter(9,r1cy[sector-1],error);
		gMinuit->GetParameter(10,r2cy[sector-1],error);
		gMinuit->GetParameter(11,r3cy[sector-1],error);

		r1xGlobal[sector-1] = r1x[sector-1]*TMath::Cos(0.3333*3.14159*(sector-1)) - 
			r1y[sector-1]*TMath::Sin(0.3333*3.14159*(sector-1));
		r1yGlobal[sector-1] = r1x[sector-1]*TMath::Sin(0.3333*3.14159*(sector-1)) + 
			r1y[sector-1]*TMath::Cos(0.3333*3.14159*(sector-1));

		r2xGlobal[sector-1] = r2x[sector-1]*TMath::Cos(0.3333*3.14159*(sector-1)) - 
			r2y[sector-1]*TMath::Sin(0.3333*3.14159*(sector-1));
		r2yGlobal[sector-1] = r2x[sector-1]*TMath::Sin(0.3333*3.14159*(sector-1)) + 
			r2y[sector-1]*TMath::Cos(0.3333*3.14159*(sector-1));

		r3xGlobal[sector-1] = r3x[sector-1]*TMath::Cos(0.3333*3.14159*(sector-1)) - 
			r3y[sector-1]*TMath::Sin(0.3333*3.14159*(sector-1));
		r3yGlobal[sector-1] = r3x[sector-1]*TMath::Sin(0.3333*3.14159*(sector-1)) + 
			r3y[sector-1]*TMath::Cos(0.3333*3.14159*(sector-1));
	}

	cout<<endl;
	cout<<"#& region sector component dx dy dz dtheta_x dtheta_y dtheta_z"<<endl;
	cout<<"1 1 0 "<<r1xGlobal[0]<<" "<<r1yGlobal[0]<<" "<<r1z[0]<<" 0 "<<r1cy[0]<<" 0"<<endl;
	cout<<"1 2 0 "<<r1xGlobal[1]<<" "<<r1yGlobal[1]<<" "<<r1z[1]<<" 0 "<<r1cy[1]<<" 0"<<endl;
	cout<<"1 3 0 "<<r1xGlobal[2]<<" "<<r1yGlobal[2]<<" "<<r1z[2]<<" 0 "<<r1cy[2]<<" 0"<<endl;
	cout<<"1 4 0 "<<r1xGlobal[3]<<" "<<r1yGlobal[3]<<" "<<r1z[3]<<" 0 "<<r1cy[3]<<" 0"<<endl;
	cout<<"1 5 0 "<<r1xGlobal[4]<<" "<<r1yGlobal[4]<<" "<<r1z[4]<<" 0 "<<r1cy[4]<<" 0"<<endl;
	cout<<"1 6 0 "<<r1xGlobal[5]<<" "<<r1yGlobal[5]<<" "<<r1z[5]<<" 0 "<<r1cy[5]<<" 0"<<endl;
	cout<<"2 1 0 "<<r2xGlobal[0]<<" "<<r2yGlobal[0]<<" "<<r2z[0]<<" 0 "<<r2cy[0]<<" 0"<<endl;
	cout<<"2 2 0 "<<r2xGlobal[1]<<" "<<r2yGlobal[1]<<" "<<r2z[1]<<" 0 "<<r2cy[1]<<" 0"<<endl;
	cout<<"2 3 0 "<<r2xGlobal[2]<<" "<<r2yGlobal[2]<<" "<<r2z[2]<<" 0 "<<r2cy[2]<<" 0"<<endl;
	cout<<"2 4 0 "<<r2xGlobal[3]<<" "<<r2yGlobal[3]<<" "<<r2z[3]<<" 0 "<<r2cy[3]<<" 0"<<endl;
	cout<<"2 5 0 "<<r2xGlobal[4]<<" "<<r2yGlobal[4]<<" "<<r2z[4]<<" 0 "<<r2cy[4]<<" 0"<<endl;
	cout<<"2 6 0 "<<r2xGlobal[5]<<" "<<r2yGlobal[5]<<" "<<r2z[5]<<" 0 "<<r2cy[5]<<" 0"<<endl;
	cout<<"3 1 0 "<<r3xGlobal[0]<<" "<<r3yGlobal[0]<<" "<<r3z[0]<<" 0 "<<r3cy[0]<<" 0"<<endl;
	cout<<"3 2 0 "<<r3xGlobal[1]<<" "<<r3yGlobal[1]<<" "<<r3z[1]<<" 0 "<<r3cy[1]<<" 0"<<endl;
	cout<<"3 3 0 "<<r3xGlobal[2]<<" "<<r3yGlobal[2]<<" "<<r3z[2]<<" 0 "<<r3cy[2]<<" 0"<<endl;
	cout<<"3 4 0 "<<r3xGlobal[3]<<" "<<r3yGlobal[3]<<" "<<r3z[3]<<" 0 "<<r3cy[3]<<" 0"<<endl;
	cout<<"3 5 0 "<<r3xGlobal[4]<<" "<<r3yGlobal[4]<<" "<<r3z[4]<<" 0 "<<r3cy[4]<<" 0"<<endl;
	cout<<"3 6 0 "<<r3xGlobal[5]<<" "<<r3yGlobal[5]<<" "<<r3z[5]<<" 0 "<<r3cy[5]<<" 0"<<endl;
}