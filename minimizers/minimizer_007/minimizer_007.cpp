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
	
	double r1_x_list[37] = {-807.27, -800.889, -809.248, -816.22, -817.14, -826.289, -694.938, 
		-710.092, -716.828, -721.028, -726.41, -719.686, 421.734, 409.551, 394.098, 394.077, 
		382.78, 378.653, 464.451, 448.826, 429.313, 443.476, 446.91, 427.406, -121.641, -118.224, 
		-137.549, -155.2, -158.633, -185.558, -114.223, -139.559, -141.715, -160.093, -160.43, 
		-160.622, 5.6365};

	double r2_x_list[37] = {1298.442, 1272.766, 1207.232, 1261.492, 1289.027, 1213.194, 1098.59, 
		1129.902, 1151.636, 1170.459, 1179.08, 1098.308, -596.412, -594.838, -647.877, -682.722, 
		-589.469, -685.217, -762.512, -713.423, -792.403, -748.739, -707.787, -721.601, 330.93, 
		255.093, 276.468, 225.857, 232.001, 163.374, 265.465, 199.463, 232.22, 218.428, 
		209.436, 231.782, 12.8148};

	double r3_x_list[37] = {-486.79, -491.267, -507.812, -489.659, -469.886, -466.532, -521.681, 
		-454.296, -479.921, -455.671, -473.201, -456.781, 217.859, 199.763, 283.092, 290.095, 
		266.577, 269.469, 274.607, 325.29, 338.966, 281.897, 352.062, 301.629, -145.534, -118.677, 
		-132.98, -58.38, -79.559, -70.447, -167.816, -153.239, -119.16, -119.713, -80.39, 
		-29.873, -16.013};

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

	double r1_z_list[37] = {282.497, 225.032, 250.514, 248.706, 223.828, 264.209, 210.743, 225.042, 
		233.736, 263.851, 211.229, 204.497, -141.965, -144.308, -122.017, -114.628, -141.512, 
		-109.734, -198.558, -141.143, -162.019, -172.812, -115.591, -157.997, 24.373, 6.065, 
		55.833, 63.735, 41.508, 54.137, 25.061, 46.097, 84.234, 26.904, 59.644, 72.763, -2.691};

	double r2_z_list[37] = {-504.067, -541.329, -506.987, -534.182, -560.742, -489.027, -493.206, 
		-491.586, -471.309, -494.97, -465.769, -472.433, 274.282, 285.89, 297.034, 332.356, 264.447, 
		271.331, 312.663, 309.488, 317.684, 295.589, 307.841, 343.811, -100.967, -118.926, -82.107, 
		-78.48, -81.748, -79.395, -108.767, -150.239, -110.116, -135.364, -92.603, -104.234, -4.03398};

	double r3_z_list[37] = {251.771, 226.02, 181.496, 237.72, 180.446, 216.815, 222.397, 239.153, 209.632, 
		248.794, 233.641, 199.876, -109.567, -102.427, -136.911, -129.162, -143.405, -169.966, -124.003, 
		-156.435, -155.814, -167.367, -141.901, -171.413, 48.74, 45.534, 33.899, 57.88, 21.206, 18.338, 
		72.435, 49.951, 48.159, 27.598, 50.422, 57.18, -0.711833};

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
	
	double sector_1_list[37] = {781.544, 784.852, 768.201, 815.635, 775.961, 780.519, 733.896, 
		798.151, 769.419, 754.405, 751.823, 782.142, -336.783, -394.656, -404.438, -415.183, 
		-426.02, -421.522, -505.034, -484.882, -517.292, -498.007, -492.104, -464.73, 148.677, 
		355.361, 161.189, 260.91, 57.437, 195.563, 178.475, 280.678, 138.145, 271.049, 134.447, 
		231.963, 5.462};
	double sector_2_list[37] = {1706.148, 1708.183, 1714.892, 1673.924, 1633.955, 1656.556, 1485.9, 
		1504.824, 1496.487, 1538.15, 1505.756, 1489.233, -811.098, -836.736, -899.383, -883.252, 
		-959.107, -962.557, -959.467, -947.986, -1017.387, -984.072, -1053.639, -1030.426, 380.704, 
		519.994, 316.195, 396.895, 219.935, 327.228, 312.225, 458.438, 317.127, 414.738, 261.939, 
		328.532, 19.245};
	double sector_3_list[37] = {281.72, 245.431, 262.335, 261.809, 267.734, 283.299, 164.392, 
		218.174, 220.533, 236.573, 202.372, 277.184, -136.815, -171.829, -200.377, -153.021, 
		-215.267, -147.416, -167.72, -137.838, -188.078, -97.812, -112.597, -103.468, 38.905, 
		259.512, 33.671, 141.875, -7.392, 102.395, -43.869, 98.423, -18.823, 123.799, 13.642, 
		154.033, -14.711};
	double sector_4_list[37] = {908.78, 907.148, 920.875, 919.955, 946.844, 945.064, 817.655, 
		787.758, 807.899, 810.55, 886.659, 837.446, -504.548, -423.231, -495.319, -512.922, 
		-458.045, -480.55, -562.8, -527.402, -551.445, -468.413, -508.539, -518.825, 154.887, 
		272.235, 128.694, 287.909, 82.545, 223.515, 112.344, 205.871, 130.992, 264.563, 139.48, 
		278.498, -18.594};
	double sector_5_list[37] = {-840.898, -908.342, -812.198, -834.031, -874.087, -837.839, 
		-730.864, -709.266, -717.59, -717.144, -716.098, -683.862, 424.488, 468.599, 453.047, 
		493.865, 392.293, 458.233, 446.528, 480.369, 438.952, 469.268, 535.586, 430.247, 
		-191.583, -78.683, -219.895, -96.63, -235.688, -86.146, -259.596, -149.721, -185.019, 
		-34.312, -99.473, -43.97, -32.732};
	double sector_6_list[37] = {-300.886, -308.052, -300.543, -335.538, -316.738, -318.066, 
		-243.158, -256.605, -310.65, -253.308, -274.449, -253.541, 157.136, 182.841, 130.301, 
		118.23, 153.894, 194.696, 147.946, 172.199, 162.303, 189.335, 166.331, 154.405, 
		-27.971, 92.512, -96.069, 10.485, -106.585, -48.341, -151.276, 14.345, -123.774, 47.62, 
		-13.617, 93.151, -14.989};

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
	for (int i = 0; i < 37; ++i) {
		// each value weighted equally, so no division by uncertainty, maybe adjust later?
		chisq += pow(calculate_residual(cur_sector, i, par), 2);
		// cout << i << endl;
	}
	f = chisq; 
}

void minimizer_007() {

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