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
	
	double r1_x_list_1[37] = {-505.001, -512.293, -500.272, -535.171, -550.866, -539.737, -458.869, -466.902, -437.394, -489.85, -456.276, -479.117, 332.769, 325.055, 309.96, 334.425, 345.505, 355.667, 378.054, 385.396, 310.398, 360.113, 347.002, 297.97, -69.659, -98.622, -113.262, -67.272, -149.987, -176.605, -56.542, -71.203, -140.584, -122.105, -142.023, -154.37, 5.6365};
	double r1_x_list_2[37] = {-191.199, -133.108, -112.802, -197.386, -190.174, -182.158, -323.303, -293.819, -339.853, -334.513, -365.434, -318.816, 354.923, 317.543, 310.6, 291.227, 293.708, 292.124, 328.176, 359.275, 312.506, 262.026, 355.28, 281.367, -128.917, -127.925, -156.551, -114.281, -132.412, -148.386, -114.225, -112.778, -133.661, -85.069, -111.371, -114.05, 5.6365};
	double r1_x_list_3[37] = {-580.81, -627.703, -648.858, -637.952, -645.586, -646.711, -576.717, -554.504, -548.852, -592.196, -594.696, -593.731, 361.492, 279.252, 397.584, 338.094, 316.296, 346.237, 411.372, 387.204, 377.427, 451.718, 374.321, 359.014, -3.385, -94.73, -122.256, -101.9, -139.304, -158.231, -89.602, -80.454, -81.313, -148.348, -107.204, -163.166, 5.6365};
	double r1_x_list_4[37] = {-521.531, -514.146, -551.507, -562.448, -566.854, -658.821, -511.555, -476.729, -551.108, -516.214, -513.194, -530.739, 353.245, 255.688, 358.71, 276.597, 310.108, 261.724, 338.294, 360.217, 335.849, 310.345, 408.098, 325.636, -119.689, -163.87, -60.934, -171.563, -155.159, -206.012, -84.724, -114.634, -148.252, -132.477, -180.708, -90.681, 5.6365};
	double r1_x_list_5[37] = {-680.412, -611.466, -719.723, -642.271, -697.947, -621.999, -585.904, -549.388, -542.159, -535.309, -618.313, -543.761, 352.262, 321.849, 325.547, 358.528, 252.887, 352.792, 318.176, 330.667, 447.205, 404.42, 357.63, 395.534, -143.007, -169.91, -123.235, -182.462, -196.937, -222.53, -14.98, -73.572, -83.634, -15.693, -197.135, -196.205, 5.6365};
	double r1_x_list_6[37] = {-658.809, -650.434, -716.041, -678.166, -698.385, -694.46, -713.444, -723.451, -726.747, -698.728, -655.6, -678.718, 449.951, 419.799, 352.017, 410.711, 438.722, 355.636, 386.403, 465.285, 420.664, 398.834, 492.521, 367.818, -164.483, -117.555, -135.782, -96.843, -59.855, -146.458, -88.392, -121.212, -191.272, -148.33, -104.319, -70.116, 5.6365};

	double r2_x_list_1[37] = {1175.871, 1148.276, 1170.048, 1221.937, 1120.399, 1218.445, 1059.627, 1105.503, 1059.819, 1062.268, 1100.386, 1035.38, -508.528, -538.916, -600.714, -570.475, -553.391, -615.211, -711.738, -690.913, -778.1, -749.605, -731.377, -754.74, 260.044, 287.85, 240.673, 210.527, 177.545, 170.38, 287.739, 259.946, 225.108, 198.572, 221.368, 165.57, 12.8148};
	double r2_x_list_2[37] = {1002.258, 1035.818, 1040.802, 968.433, 895.307, 1028.052, 901.978, 867.451, 936.635, 965.958, 910.075, 902.595, -571.809, -533.893, -589.293, -609.485, -632.967, -559.075, -621.619, -672.977, -630.834, -701.663, -672.238, -705.666, 240.484, 206.901, 217.872, 220.226, 200.192, 200.241, 286.258, 220.327, 265.463, 263.854, 301.487, 262.803, 12.8148};
	double r2_x_list_3[37] = {1117.853, 1010.575, 1034.803, 1086.667, 1065.863, 1092.627, 905.391, 1007.603, 970.211, 1001.867, 995.124, 1022.385, -637.678, -700.165, -590.028, -575.584, -690.881, -592.113, -682.906, -630.4, -638.492, -609.234, -594.97, -709.007, 273.068, 260.601, 238.34, 223.744, 155.541, 166.393, 242.715, 231.73, 245.66, 178.853, 218.341, 118.649, 12.8148};
	double r2_x_list_4[37] = {1311.547, 1259.074, 1274.789, 1232.711, 1302.906, 1289.969, 1117.717, 1070.752, 1107.105, 989.453, 1163.026, 1042.097, -653.961, -681.501, -624.4, -739.636, -609.334, -690.267, -631.38, -665.181, -641.781, -665.08, -724.565, -791.405, 193.842, 275.18, 342.703, 89.099, 256.299, 224.326, 221.539, 244.965, 249.661, 261.316, 215.509, 157.805, 12.8148};
	double r2_x_list_5[37] = {614.502, 632.721, 756.679, 638.375, 621.013, 656.863, 825.837, 823.165, 964.462, 711.769, 855.048, 892.711, -557.049, -557.973, -702.224, -693.781, -563.735, -760.88, -652.111, -703.269, -709.966, -664.678, -593.456, -662.906, 341.102, 311.611, 251.22, 268.302, 247.82, 190.636, 303.771, 191.54, 285.117, 208.876, 275.321, 114.734, 12.8148};
	double r2_x_list_6[37] = {1156.672, 1198.881, 1097.647, 1088.288, 1126.83, 1090.857, 1022.05, 1113.612, 948.374, 897.159, 1007.267, 987.277, -667.42, -699.293, -588.806, -580.12, -695.016, -769.249, -634.162, -622.288, -657.834, -727.643, -602.681, -720.095, 236.834, 243.044, 231.471, 198.905, 246.631, 227.323, 245.486, 243.349, 297.874, 354.882, 268.288, 145.174, 12.8148};

	double r3_x_list_1[37] = {-535.802, -508.419, -502.146, -550.819, -515.607, -484.727, -483.031, -487.569, -489.815, -448.098, -435.975, -487.119, 253.966, 256.587, 255.637, 277.616, 302.784, 320.523, 311.756, 322.253, 304.77, 367.675, 395.346, 395.097, -162.262, -158.823, -98.196, -46.543, -116.856, -48.662, -180.561, -129.992, -165.04, -102.532, -69.656, -27.281, -16.013};
	double r3_x_list_2[37] = {-175.8, -150.74, -145.706, -172.264, -172.396, -171.776, -299.782, -307.12, -301.343, -331.009, -293.903, -285.431, 201.797, 249.383, 279.81, 290.758, 317.969, 313.927, 285.794, 317.361, 300.917, 242.022, 345.079, 321.771, -253.446, -122.724, -141.109, -94.336, -73.757, -30.429, -139.371, -182.223, -134.063, -107.099, -69.596, -34.146, -16.013};
	double r3_x_list_3[37] = {-614.678, -623.864, -662.307, -572.835, -667.449, -645.852, -582.821, -563.088, -548.798, -520.138, -504.261, -522.345, 318.9, 307.157, 314.428, 299.661, 273.22, 338.645, 397.787, 394.502, 381.36, 424.627, 424.045, 413.377, -137.434, -89.12, -71.031, -102.951, -94.668, -101.84, -214.648, -131.916, -146.614, -110.004, -95.989, -98.77, -16.013};
	double r3_x_list_4[37] = {-189.792, -231.733, -274.046, -240.716, -259.714, -279.231, -337.003, -193.474, -227.245, -218.791, -180.169, -214.08, 50.867, 73.741, 138.141, 130.134, 70.587, 95.399, 116.83, 208.681, 92.075, 174.977, 232.738, 202.814, -75.054, -107.158, -19.389, -95.017, -33.708, -20.156, -102.901, -132.11, -59.945, -21.872, -127.022, -21.657, -16.013};
	double r3_x_list_5[37] = {-538.637, -439.81, -433.822, -464.448, -496.338, -453.013, -501.069, -504.238, -502.765, -508.727, -513.346, -527.686, 195.565, 220.305, 223.774, 292.405, 270.036, 317.42, 243.292, 383.084, 299.767, 311.472, 355.891, 370.248, -138.002, -129.024, -131.674, -112.495, -106.741, -103.908, -153.962, -147.129, -133.632, -131.154, -54.229, -119.79, -16.013};
	double r3_x_list_6[37] = {-166.688, -308.654, -244.933, -265.579, -292.446, -152.013, -271.968, -235.537, -252.347, -287.741, -275.489, -275.648, 79.601, 124.527, 78.362, 143.888, 90.66, 131.728, 210.82, 156.631, 149.128, 122.249, 222.612, 109.618, -79.544, -38.404, -37.099, 25.218, 36.96, 42.583, -104.453, -95.205, -130.316, -74.703, -79.64, -0.594, -16.013};

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

	double r1_z_list_1[37] = {315.777, 283.909, 356.601, 326.848, 309.154, 344.948, 199.929, 230.118, 263.306, 242.893, 258.409, 254.86, -103.447, -141.061, -89.879, -136.446, -81.443, -129.315, -143.609, -171.661, -230.904, -145.295, -175.737, -153.869, 69.437, 60.418, 61.272, 100.261, 88.7, 7.338, 9.073, 43.439, 38.052, 42.252, 47.925, 103.815, -2.691};
	double r1_z_list_2[37] = {103.793, 142.149, 166.977, 186.371, 177.374, 144.711, 145.033, 120.649, 115.826, 123.074, 167.042, 195.584, -172.31, -123.849, -193.705, -114.853, -138.78, -146.124, -114.408, -128.303, -92.805, -142.751, -102.036, -102.429, 18.384, 44.307, 72.234, 71.424, 75.753, 73.889, 40.098, 46.907, 25.169, 84.528, 33.316, 71.88, -2.691};
	double r1_z_list_3[37] = {161.6, 151.923, 151.252, 189.174, 176.663, 171.921, 107.518, 209.802, 200.155, 158.905, 168.696, 175.681, -100.858, -140.325, -130.504, -144.171, -101.967, -86.048, -105.013, -87.828, -124.082, -99.259, -105.931, -78.168, 69.759, 79.29, 79.279, 99.725, -30.544, 30.551, -0.68, 92.542, 34.924, 52.595, 25.844, 32.706, -2.691};
	double r1_z_list_4[37] = {228.889, 263.236, 172.687, 275.953, 180.589, 240.916, 183.41, 231.789, 248.235, 214.57, 310.377, 250.93, -186.32, -153.734, -198.938, -156.88, -110.045, -161.935, -115.579, -109.757, -76.596, -110.195, -121.414, -154.565, -14.128, -0.965, 69.829, 43.637, 73.231, 1.239, 56.397, 76.628, -31.798, 55.132, -22.736, 31.412, -2.691};
	double r1_z_list_5[37] = {134.47, 191.109, 286.736, 135.453, 195.719, 224.308, 98.043, 164.499, 188.176, 93.756, 135.855, 212.222, -188.412, -171.621, -200.588, -158.511, -94.708, -49.733, -175.306, -49.103, -57.456, -112.539, -148.531, 35.542, 16.184, -3.804, 108.63, 74.389, 45.205, 21.744, 39.833, -41.912, 77.346, -2.924, 132.562, 23.863, -2.691};
	double r1_z_list_6[37] = {299.532, 207.45, 140.04, 189.224, 243.304, 153.77, 206.147, 256.823, 186.112, 228.935, 187.653, 236.774, -95.07, -188.536, -133.937, -218.584, -152.524, -126.073, -202.094, -206.662, -217.774, -192.587, -85.14, -155.474, 15.752, 133.609, 140.481, 58.969, 33.47, 65.227, -66.05, 92.119, 65.536, 28.13, 106.517, 111.105, -2.691};

	double r2_z_list_1[37] = {-419.949, -420.444, -403.831, -431.648, -446.693, -418.1, -372.577, -370.752, -416.742, -377.275, -359.729, -405.85, 297.296, 276.167, 268.759, 282.141, 303.016, 300.611, 318.717, 303.994, 292.917, 340.51, 345.607, 304.718, -94.135, -123.553, -101.893, -88.305, -103.652, -81.745, -105.515, -106.765, -145.904, -133.724, -106.556, -126.105, -4.03398};
	double r2_z_list_2[37] = {-208.727, -147.295, -200.481, -175.248, -174.922, -161.954, -186.707, -200.93, -247.814, -193.545, -224.945, -243.726, 242.206, 255.376, 276.654, 271.48, 249.535, 270.24, 201.254, 231.131, 256.435, 209.936, 241.823, 241.447, -160.901, -79.382, -89.678, -86.4, -68.534, -71.467, -84.832, -91.819, -68.844, -99.8, -94.796, -31.441, -4.03398};
	double r2_z_list_3[37] = {-494.492, -481.157, -512.39, -488.915, -507.593, -510.222, -450.241, -447.722, -418.233, -449.951, -406.017, -409.419, 303.488, 271.204, 285.444, 285.463, 290.523, 296.538, 331.611, 334.586, 283.175, 367.788, 319.173, 345.368, -84.872, -103.369, -119.657, -137.244, -97.698, -66.244, -143.842, -97.918, -157.063, -137.889, -101.475, -124.278, -4.03398};
	double r2_z_list_4[37] = {-322.776, -358.379, -368.919, -352.864, -397.149, -395.889, -357.829, -288.203, -344.665, -321.211, -293.139, -329.13, 241.801, 228.927, 305.827, 264.983, 260.69, 212.813, 222.457, 268.143, 219.574, 233.679, 307.687, 266.88, -123.319, -191.445, -121.508, -113.978, -27.216, -119.567, -69.493, -93.743, -0.916, 25.669, -105.834, -77.252, -4.03398};
	double r2_z_list_5[37] = {-487.292, -442.771, -449.618, -391.463, -434.198, -455.712, -386.21, -304.709, -409.338, -326.207, -397.23, -392.299, 167.776, 192.658, 271.839, 228.789, 228.016, 215.871, 228.684, 242.581, 316.722, 248.855, 334.78, 281.405, -92.537, -111.961, -118.553, -108.719, -91.357, -134.474, -126.273, -89.9, -122.377, -119.207, -80.245, -153.307, -4.03398};
	double r2_z_list_6[37] = {-423.583, -509.471, -456.71, -376.332, -459.04, -404.075, -401.327, -461.97, -424.842, -383.26, -436.176, -403.01, 168.883, 249.263, 241.404, 315.313, 222.715, 267.241, 276.523, 254.405, 282.378, 265.866, 332.493, 218.818, -110.835, -21.618, -35.312, -71.121, -55.947, -50.819, -129.787, -133.867, -148.839, -148.32, -61.751, -55.093, -4.03398};

	double r3_z_list_1[37] = {218.718, 219.612, 252.029, 241.124, 211.582, 224.508, 202.031, 211.271, 233.219, 226.632, 253.039, 234.401, -68.574, -86.996, -92.43, -92.144, -112.917, -114.465, -99.503, -150.286, -193.637, -139.362, -168.494, -149.907, 68.887, 80.769, 87.837, 48.437, 56.504, -6.346, 65.365, 66.839, -14.338, 8.931, 25.852, 23.802, -0.711833};
	double r3_z_list_2[37] = {127.978, 113.27, 145.673, 144.98, 115.861, 119.335, 188.417, 127.904, 108.278, 169.616, 133.226, 120.023, -109.328, -110.42, -67.235, -125.016, -129.257, -108.607, -102.241, -126.559, -106.387, -143.195, -166.459, -172.571, 67.305, 61.411, 28.801, 76.527, 39.016, 15.494, 93.473, 53.88, 59.18, 70.099, 44.404, 36.916, -0.711833};
	double r3_z_list_3[37] = {134.169, 140.781, 135.091, 183.531, 140.889, 178.838, 188.353, 204.442, 231.321, 220.184, 174.133, 166.211, -93.815, -82.808, -102.687, -101.171, -143.037, -59.357, -116.742, -110.162, -167.422, -66.727, -143.903, -148.564, 126.31, 102.312, 45.655, 49.602, 40.022, 38.04, 24.974, 30.002, -11.161, 39.476, 37.134, -11.718, -0.711833};
	double r3_z_list_4[37] = {266.091, 238.579, 226.784, 236.543, 220.919, 169.548, 205.752, 176.75, 193.199, 225.868, 169.983, 248.271, -177.162, -153.379, -154.101, -189.567, -207.99, -159.187, -52.829, -109.823, -132.906, -127.039, -88.144, -158.686, 64.674, 95.736, 68.739, 47.283, 11.995, 37.41, 5.648, 53.187, 62.458, 95.552, -34.019, 15.511, -0.711833};
	double r3_z_list_5[37] = {210.465, 129.108, 193.553, 77.15, 149.648, 101.4, 178.095, 206.718, 202.796, 153.414, 141.189, 78.883, -139.493, -66.47, -180.927, -47.093, -95.144, -53.041, -151.651, -86.482, -82.494, -75.53, -146.112, -123.562, -0.42, 55.104, 62.055, 22.091, 43.938, 21.749, 115.685, 0.42, 19.9, 67.66, 1.958, -26.047, -0.711833};
	double r3_z_list_6[37] = {240.82, 342.358, 272.32, 299.242, 270.456, 264.417, 167.857, 188.185, 246.783, 186.074, 173.352, 261.657, -34.506, -181.349, -136.889, -182.761, -113.618, -242.975, -180.709, -37.295, -138.266, -202.091, -180.806, -159.209, -1.168, 51.198, 43.11, 47.312, 85.79, 50.369, 7.073, -2.866, 52.082, 83.752, 31.553, 21.902, -0.711833};

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
	
	double sector_1_list[37] = {651.117, 658.16, 660.858, 660.066, 661.961, 675.024, 608.946, 647.163, 
		626.443, 650.644, 628.014, 660.365, -303.922, -324.011, -347.121, -364.902, -371.123, -373.86, 
		-453.43, -430.314, -445.232, -405.797, -435.144, -414.478, 127.981, 271.178, 98.831, 192.185, 
		7.493, 158.859, 101.81, 225.055, 80.446, 188.761, 62.622, 174.989, 5.462};
	double sector_2_list[37] = {1355.748, 1335.636, 1335.087, 1318.046, 1295.666, 1341.279, 1193.768, 
		1245.754, 1225.892, 1224.909, 1214.791, 1232.697, -705.881, -738.611, -789.76, -791.026, 
		-827.147, -818.327, -861.051, -868.23, -895.061, -872.071, -923.808, -893.911, 338.461, 
		461.891, 245.429, 339.012, 157.093, 261.662, 288.534, 396.842, 241.925, 359.396, 209.359, 
		298.492, 19.245};
	double sector_3_list[37] = {202.894, 195.993, 217.491, 229.786, 221.629, 241.384, 138.719, 166.201, 
		174.595, 185.919, 180.182, 231.12, -132.048, -127.85, -139.813, -133.733, -124.724, -132.994, 
		-132.709, -102.478, -167.636, -69.207, -130.671, -99.371, 27.947, 216.224, 6.921, 106.26, 
		-15.972, 68.874, -46.409, 98.525, -29.513, 117.279, 10.974, 96.071, -14.711};
	double sector_4_list[37] = {781.599, 794.351, 787.765, 804.744, 811.201, 841.336, 659.617, 704.157, 
		705.059, 694.973, 744.036, 733.918, -441.353, -446.913, -454.855, -416.721, -429.821, -455.533, 
		-465.331, -461.879, -479.14, -421.054, -452.128, -457.484, 138.004, 252.373, 96.936, 212.289, 
		58.543, 192.383, 53.911, 171.629, 65.734, 240.671, 129.12, 224.256, -18.594};
	double sector_5_list[37] = {-873.28, -869.637, -871.812, -856.111, -851.563, -858.39, -780.904, 
		-746.176, -720.724, -766.859, -744.688, -730.194, 432.709, 490.616, 439.327, 459.279, 485.782, 
		462.429, 492.423, 493.21, 424.286, 467.829, 476.396, 463.897, -191.547, -39.37, -210.639, 
		-114.338, -253.529, -103.263, -233.168, -128.705, -189.408, -79.476, -123.689, -62.929, -32.732};
	double sector_6_list[37] = {-310.708, -336.156, -337.541, -324.305, -330.548, -317.124, -278.92, 
		-322.028, -309.027, -319.611, -310.79, -284.517, 186.586, 172.947, 184.187, 138.775, 145.426, 
		190.361, 191.171, 157.544, 177.735, 220.518, 187.763, 177.405, -34.317, 45.327, -96.254, 
		-39.334, -184.875, -41.436, -111.179, -35.776, -149.125, -19.032, -85.607, 51.309, -14.989};

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

	double value;
	if (sector==1) {
		value = ((1/shift)*(par[0]*r1_x_list_1[j][i]+par[1]*r1_y_list[i]+par[2]*r1_z_list_1[i]+
			par[3]*r2_x_list_1[i]+par[4]*r2_y_list[i]+par[5]*r2_z_list_1[i]+par[6]*r3_x_list_1[i]+
			par[7]*r3_y_list[i]+par[8]*r3_z_list_1[i])+(1/rot)*(par[9]*r1_cy_list[i]+
			par[10]*r2_cy_list[i]+par[11]*r3_cy_list[i]))-current_position+desired_positions[i];
	} else if (sector == 2) {
		value = ((1/shift)*(par[0]*r1_x_list_2[j][i]+par[1]*r1_y_list[i]+par[2]*r1_z_list_2[i]+
			par[3]*r2_x_list_2[i]+par[4]*r2_y_list[i]+par[5]*r2_z_list_2[i]+par[6]*r3_x_list_2[i]+
			par[7]*r3_y_list[i]+par[8]*r3_z_list_2[i])+(1/rot)*(par[9]*r1_cy_list[i]+
			par[10]*r2_cy_list[i]+par[11]*r3_cy_list[i]))-current_position+desired_positions[i];
	} else if (sector == 3) {
		value = ((1/shift)*(par[0]*r1_x_list_3[j][i]+par[1]*r1_y_list[i]+par[2]*r1_z_list_3[i]+
			par[3]*r2_x_list_3[i]+par[4]*r2_y_list[i]+par[5]*r2_z_list_3[i]+par[6]*r3_x_list_3[i]+
			par[7]*r3_y_list[i]+par[8]*r3_z_list_3[i])+(1/rot)*(par[9]*r1_cy_list[i]+
			par[10]*r2_cy_list[i]+par[11]*r3_cy_list[i]))-current_position+desired_positions[i];
	} else if (sector == 4) {
		value = ((1/shift)*(par[0]*r1_x_list_4[j][i]+par[1]*r1_y_list[i]+par[2]*r1_z_list_4[i]+
			par[3]*r2_x_list_4[i]+par[4]*r2_y_list[i]+par[5]*r2_z_list_4[i]+par[6]*r3_x_list_4[i]+
			par[7]*r3_y_list[i]+par[8]*r3_z_list_4[i])+(1/rot)*(par[9]*r1_cy_list[i]+
			par[10]*r2_cy_list[i]+par[11]*r3_cy_list[i]))-current_position+desired_positions[i];
	} else if (sector == 5) {
		value = ((1/shift)*(par[0]*r1_x_list_5[j][i]+par[1]*r1_y_list[i]+par[2]*r1_z_list_5[i]+
			par[3]*r2_x_list_5[i]+par[4]*r2_y_list[i]+par[5]*r2_z_list_5[i]+par[6]*r3_x_list_5[i]+
			par[7]*r3_y_list[i]+par[8]*r3_z_list_5[i])+(1/rot)*(par[9]*r1_cy_list[i]+
			par[10]*r2_cy_list[i]+par[11]*r3_cy_list[i]))-current_position+desired_positions[i];
	} else {
		value = ((1/shift)*(par[0]*r1_x_list_6[j][i]+par[1]*r1_y_list[i]+par[2]*r1_z_list_6[i]+
			par[3]*r2_x_list_6[i]+par[4]*r2_y_list[i]+par[5]*r2_z_list_6[i]+par[6]*r3_x_list_6[i]+
			par[7]*r3_y_list[i]+par[8]*r3_z_list_6[i])+(1/rot)*(par[9]*r1_cy_list[i]+
			par[10]*r2_cy_list[i]+par[11]*r3_cy_list[i]))-current_position+desired_positions[i];
	}
	
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

void minimizer_008() {

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