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


double calculate_residual(int sector, int j, int i, double *par) {
	// nominal shifts by Raffaella de Vita were shifts of 0.2 cm and rotations of 0.1 degrees
	// may need to update these values in future
	double shift = 0.2;
	double rot = 0.2;

	// values in these lists are in micrometers
	// rotations around x and z axis are negligible when compared to other rotations so we do
	// not use them
	
	double r1_x_list[4][37] = {{-954.619, -946.066, -958.763, -966.726, -951.533, -979.925, -841.613, -855.35, -854.834, -859.952, -850.332, -879.696, 458.156, 448.408, 410.323, 424.585, 413.646, 415.68, 515.651, 478.013, 483.028, 459.474, 485.225, 460.979, -129.066, -127.258, -156.379, -165.307, -174.931, -206.653, -99.211, -138.069, -158.912, -151.556, -152.214, -200.257, 5.6365},
		{-835.119, -851.249, -869.089, -850.566, -859.652, -893.704, -733.454, -756.395, -756.504, -741.182, -761.921, -759.207, 431.097, 411.593, 424.437, 389.955, 402.05, 383.078, 472.672, 487.659, 467.957, 464.693, 451.853, 434.817, -119.847, -108.029, -148.784, -151.313, -162.622, -190.526, -107.451, -123.792, -153.13, -154.169, -159.471, -172.637, 5.6365},
		{-679.974, -671.921, -689.795, -706.853, -697.021, -702.561, -590.818, -587.028, -605.312, -630.691, -606.522, -621.403, 350.841, 367.543, 350.748, 361.976, 330.637, 343.74, 405.639, 439.375, 387.635, 432.939, 448.223, 390.309, -124.323, -117.305, -114.513, -153.395, -145.28, -170.069, -125.149, -139.816, -164.092, -154.657, -165.726, -146.252, 5.6365},
		{-609.734, -684.869, -661.838, -683.353, -712.052, -717.693, -570.575, -601.073, -650.44, -579.502, -585.315, -568.259, 375.559, 353.174, 326.237, 349.94, 357.829, 381.432, 408.32, 415.716, 404.281, 431.751, 419.68, 387.663, -57.015, -111.091, -117.321, -101.739, -145.602, -125.235, -45.887, -88.981, -109.711, -171.037, -79.915, -113.318, 5.6365}};
	// double r1_x_list[4][37] = {{-954.619, -946.066, -958.763, -966.726, -951.533, -979.925, -841.613, -855.35, -854.834, -859.952, -850.332, -879.696, 458.156, 448.408, 410.323, 424.585, 413.646, 415.68, 515.651, 478.013, 483.028, 459.474, 485.225, 460.979, -129.066, -127.258, -156.379, -165.307, -174.931, -206.653, -99.211, -138.069, -158.912, -151.556, -152.214, -200.257, 5.6365},
	// 	{-835.119, -851.249, -869.089, -850.566, -859.652, -893.704, -733.454, -756.395, -756.504, -741.182, -761.921, -759.207, 431.097, 411.593, 424.437, 389.955, 402.05, 383.078, 472.672, 487.659, 467.957, 464.693, 451.853, 434.817, -119.847, -108.029, -148.784, -151.313, -162.622, -190.526, -107.451, -123.792, -153.13, -154.169, -159.471, -172.637, 5.6365},
	// 	{-679.974, -671.921, -689.795, -706.853, -697.021, -702.561, -590.818, -587.028, -605.312, -630.691, -606.522, -621.403, 350.841, 367.543, 350.748, 361.976, 330.637, 343.74, 405.639, 439.375, 387.635, 432.939, 448.223, 390.309, -124.323, -117.305, -114.513, -153.395, -145.28, -170.069, -125.149, -139.816, -164.092, -154.657, -165.726, -146.252, 5.6365},
	// 	{-609.734, -684.869, -661.838, -683.353, -712.052, -717.693, -570.575, -601.073, -650.44, -579.502, -585.315, -568.259, 375.559, 353.174, 326.237, 349.94, 357.829, 381.432, 408.32, 415.716, 404.281, 431.751, 419.68, 387.663, -57.015, -111.091, -117.321, -101.739, -145.602, -125.235, -45.887, -88.981, -109.711, -171.037, -79.915, -113.318, 5.6365}};
	// double r1_x_list[4][37] = {{-954.619, -946.066, -958.763, -966.726, -951.533, -979.925, -841.613, -855.35, -854.834, -859.952, -850.332, -879.696, 458.156, 448.408, 410.323, 424.585, 413.646, 415.68, 515.651, 478.013, 483.028, 459.474, 485.225, 460.979, -129.066, -127.258, -156.379, -165.307, -174.931, -206.653, -99.211, -138.069, -158.912, -151.556, -152.214, -200.257, 5.6365},
	// 	{-835.119, -851.249, -869.089, -850.566, -859.652, -893.704, -733.454, -756.395, -756.504, -741.182, -761.921, -759.207, 431.097, 411.593, 424.437, 389.955, 402.05, 383.078, 472.672, 487.659, 467.957, 464.693, 451.853, 434.817, -119.847, -108.029, -148.784, -151.313, -162.622, -190.526, -107.451, -123.792, -153.13, -154.169, -159.471, -172.637, 5.6365},
	// 	{-679.974, -671.921, -689.795, -706.853, -697.021, -702.561, -590.818, -587.028, -605.312, -630.691, -606.522, -621.403, 350.841, 367.543, 350.748, 361.976, 330.637, 343.74, 405.639, 439.375, 387.635, 432.939, 448.223, 390.309, -124.323, -117.305, -114.513, -153.395, -145.28, -170.069, -125.149, -139.816, -164.092, -154.657, -165.726, -146.252, 5.6365},
	// 	{-609.734, -684.869, -661.838, -683.353, -712.052, -717.693, -570.575, -601.073, -650.44, -579.502, -585.315, -568.259, 375.559, 353.174, 326.237, 349.94, 357.829, 381.432, 408.32, 415.716, 404.281, 431.751, 419.68, 387.663, -57.015, -111.091, -117.321, -101.739, -145.602, -125.235, -45.887, -88.981, -109.711, -171.037, -79.915, -113.318, 5.6365}};
	// double r1_x_list[4][37] = {{-954.619, -946.066, -958.763, -966.726, -951.533, -979.925, -841.613, -855.35, -854.834, -859.952, -850.332, -879.696, 458.156, 448.408, 410.323, 424.585, 413.646, 415.68, 515.651, 478.013, 483.028, 459.474, 485.225, 460.979, -129.066, -127.258, -156.379, -165.307, -174.931, -206.653, -99.211, -138.069, -158.912, -151.556, -152.214, -200.257, 5.6365},
	// 	{-835.119, -851.249, -869.089, -850.566, -859.652, -893.704, -733.454, -756.395, -756.504, -741.182, -761.921, -759.207, 431.097, 411.593, 424.437, 389.955, 402.05, 383.078, 472.672, 487.659, 467.957, 464.693, 451.853, 434.817, -119.847, -108.029, -148.784, -151.313, -162.622, -190.526, -107.451, -123.792, -153.13, -154.169, -159.471, -172.637, 5.6365},
	// 	{-679.974, -671.921, -689.795, -706.853, -697.021, -702.561, -590.818, -587.028, -605.312, -630.691, -606.522, -621.403, 350.841, 367.543, 350.748, 361.976, 330.637, 343.74, 405.639, 439.375, 387.635, 432.939, 448.223, 390.309, -124.323, -117.305, -114.513, -153.395, -145.28, -170.069, -125.149, -139.816, -164.092, -154.657, -165.726, -146.252, 5.6365},
	// 	{-609.734, -684.869, -661.838, -683.353, -712.052, -717.693, -570.575, -601.073, -650.44, -579.502, -585.315, -568.259, 375.559, 353.174, 326.237, 349.94, 357.829, 381.432, 408.32, 415.716, 404.281, 431.751, 419.68, 387.663, -57.015, -111.091, -117.321, -101.739, -145.602, -125.235, -45.887, -88.981, -109.711, -171.037, -79.915, -113.318, 5.6365}};
	// double r1_x_list[4][37] = {{-954.619, -946.066, -958.763, -966.726, -951.533, -979.925, -841.613, -855.35, -854.834, -859.952, -850.332, -879.696, 458.156, 448.408, 410.323, 424.585, 413.646, 415.68, 515.651, 478.013, 483.028, 459.474, 485.225, 460.979, -129.066, -127.258, -156.379, -165.307, -174.931, -206.653, -99.211, -138.069, -158.912, -151.556, -152.214, -200.257, 5.6365},
	// 	{-835.119, -851.249, -869.089, -850.566, -859.652, -893.704, -733.454, -756.395, -756.504, -741.182, -761.921, -759.207, 431.097, 411.593, 424.437, 389.955, 402.05, 383.078, 472.672, 487.659, 467.957, 464.693, 451.853, 434.817, -119.847, -108.029, -148.784, -151.313, -162.622, -190.526, -107.451, -123.792, -153.13, -154.169, -159.471, -172.637, 5.6365},
	// 	{-679.974, -671.921, -689.795, -706.853, -697.021, -702.561, -590.818, -587.028, -605.312, -630.691, -606.522, -621.403, 350.841, 367.543, 350.748, 361.976, 330.637, 343.74, 405.639, 439.375, 387.635, 432.939, 448.223, 390.309, -124.323, -117.305, -114.513, -153.395, -145.28, -170.069, -125.149, -139.816, -164.092, -154.657, -165.726, -146.252, 5.6365},
	// 	{-609.734, -684.869, -661.838, -683.353, -712.052, -717.693, -570.575, -601.073, -650.44, -579.502, -585.315, -568.259, 375.559, 353.174, 326.237, 349.94, 357.829, 381.432, 408.32, 415.716, 404.281, 431.751, 419.68, 387.663, -57.015, -111.091, -117.321, -101.739, -145.602, -125.235, -45.887, -88.981, -109.711, -171.037, -79.915, -113.318, 5.6365}};
	// double r1_x_list[4][37] = {{-954.619, -946.066, -958.763, -966.726, -951.533, -979.925, -841.613, -855.35, -854.834, -859.952, -850.332, -879.696, 458.156, 448.408, 410.323, 424.585, 413.646, 415.68, 515.651, 478.013, 483.028, 459.474, 485.225, 460.979, -129.066, -127.258, -156.379, -165.307, -174.931, -206.653, -99.211, -138.069, -158.912, -151.556, -152.214, -200.257, 5.6365},
	// 	{-835.119, -851.249, -869.089, -850.566, -859.652, -893.704, -733.454, -756.395, -756.504, -741.182, -761.921, -759.207, 431.097, 411.593, 424.437, 389.955, 402.05, 383.078, 472.672, 487.659, 467.957, 464.693, 451.853, 434.817, -119.847, -108.029, -148.784, -151.313, -162.622, -190.526, -107.451, -123.792, -153.13, -154.169, -159.471, -172.637, 5.6365},
	// 	{-679.974, -671.921, -689.795, -706.853, -697.021, -702.561, -590.818, -587.028, -605.312, -630.691, -606.522, -621.403, 350.841, 367.543, 350.748, 361.976, 330.637, 343.74, 405.639, 439.375, 387.635, 432.939, 448.223, 390.309, -124.323, -117.305, -114.513, -153.395, -145.28, -170.069, -125.149, -139.816, -164.092, -154.657, -165.726, -146.252, 5.6365},
	// 	{-609.734, -684.869, -661.838, -683.353, -712.052, -717.693, -570.575, -601.073, -650.44, -579.502, -585.315, -568.259, 375.559, 353.174, 326.237, 349.94, 357.829, 381.432, 408.32, 415.716, 404.281, 431.751, 419.68, 387.663, -57.015, -111.091, -117.321, -101.739, -145.602, -125.235, -45.887, -88.981, -109.711, -171.037, -79.915, -113.318, 5.6365}};
	
	double r2_x_list[4][37] = {{1822.005, 1812.622, 1830.316, 1802.596, 1790.081, 1802.511, 1629.155, 1597.752, 1612.797, 1613.637, 1601.37, 1609.771, -770.93, -795.685, -815.673, -794.888, -812.73, -822.369, -882.928, -936.239, -914.534, -939.066, -923.29, -945.192, 377.432, 336.611, 316.469, 314.177, 282.116, 239.545, 364.097, 301.658, 291.391, 311.87, 273.379, 219.102, 12.8148},
		{1680.437, 1669.602, 1689.141, 1638.262, 1625.403, 1648.99, 1454.672, 1452.766, 1440.58, 1449.577, 1455.612, 1452.268, -736.566, -773.226, -744.278, -765.707, -769.373, -795.325, -917.532, -874.861, -896.982, -881.625, -905.966, -929.949, 346.408, 333.258, 281.19, 265.259, 261.472, 223.008, 333.97, 305.304, 291.263, 289.343, 250.482, 241.105, 12.8148},
		{1411.155, 1361.925, 1408.787, 1371.051, 1437.127, 1371.537, 1232.541, 1254.849, 1278.782, 1264.753, 1262.215, 1240.227, -745.475, -711.633, -754.973, -735.859, -766.885, -728.333, -833.725, -808.871, -867.248, -819.324, -820.443, -851.233, 309.446, 299.803, 314.509, 251.913, 276.635, 227.27, 285.872, 245.449, 234.248, 277.991, 245.626, 208.567, 12.8148},
		{1344.898, 1335.937, 1263.116, 1306.399, 1278.172, 1228.855, 1051.966, 1014.691, 1088.449, 1098.686, 1168.394, 1095.702, -704.193, -712.628, -765.082, -751.244, -740.63, -701.51, -722.089, -765.779, -743.789, -767.104, -741.392, -790.385, 294.261, 205.236, 224.518, 222.913, 220.757, 202.173, 361.227, 360.782, 293.728, 299.488, 276.171, 236.517, 12.8148}};
	
	double r3_x_list[4][37] = {{-696.578, -666.28, -667.467, -672.494, -644.678, -654.338, -634.573, -641.736, -613.775, -614.326, -602.036, -609.393, 253.022, 271.229, 253.612, 296.388, 301.79, 331.202, 343.342, 325.961, 349.957, 353.846, 401.524, 413.582, -164.063, -122.871, -117.01, -89.03, -80.054, -71.815, -149.583, -157.321, -135.397, -105.722, -74.695, -78.252, -16.013},
		{-620.178, -626.019, -615.078, -596.538, -584.093, -601.987, -566.302, -570.657, -556.958, -536.938, -536.864, -528.777, 238.466, 237.446, 279.3, 263.589, 296.416, 306.069, 305.987, 350.317, 346.263, 377.848, 381.707, 385.792, -162.313, -114.174, -119.767, -78.424, -54.85, -65.465, -164.209, -150.663, -137.636, -115.014, -81.577, -62.892, -16.013},
		{-532.392, -512.554, -509.376, -512.301, -480.486, -491.624, -493.79, -476.444, -461.419, -486.233, -444.969, -446.163, 190.488, 227.87, 219.127, 267.123, 257.808, 296.648, 276.721, 322.445, 313.606, 355.088, 404.847, 362.59, -143.297, -115.58, -90.648, -96.992, -66.313, -47.582, -171.648, -164.81, -153.772, -105.633, -82.285, -45.811, -16.013},
		{-451.785, -537.528, -483.214, -509.628, -502.994, -497.914, -453.539, -500.636, -527.828, -422.416, -427.005, -406.065, 212.338, 215.411, 203.996, 254.714, 291.768, 314.783, 278.497, 304.25, 300.024, 367.623, 399.788, 364.48, -96.567, -120.817, -95.137, -70.759, -85.258, -34.854, -119.792, -77.416, -86.314, -146.161, -21.752, -41.656, -16.013}};
	
	double r1_z_list[4][37] = {{92.782, 103.972, 95.057, 88.389, 99.859, 102.086, 80.175, 75.939, 89.045, 82.064, 88.194, 89.0, -40.7, -45.663, -61.298, -32.49, -43.682, -34.562, -43.947, -73.046, -42.327, -46.891, -32.907, -41.434, 28.134, 31.078, 16.812, 25.118, 14.463, 2.397, 31.179, 5.756, 6.776, 36.237, 30.463, -2.672, -2.691},
		{174.737, 160.575, 164.624, 160.653, 156.654, 155.167, 122.92, 125.578, 126.881, 133.503, 132.22, 137.32, -63.768, -85.213, -51.202, -74.626, -55.545, -71.562, -102.638, -62.663, -77.563, -60.05, -74.309, -71.601, 18.444, 36.249, 14.388, 31.118, 33.785, 24.695, 28.778, 24.656, 21.082, 28.667, 20.292, 26.771, -2.691},
		{285.61, 275.937, 296.446, 286.365, 313.831, 283.779, 255.74, 270.097, 281.152, 261.693, 271.516, 260.576, -194.794, -173.812, -172.401, -145.17, -162.513, -126.8, -203.479, -173.461, -182.556, -148.376, -119.719, -158.97, 46.329, 50.755, 71.321, 53.023, 71.721, 67.323, 16.673, 18.745, 12.876, 53.43, 40.746, 39.711, -2.691},
		{438.04, 398.915, 352.112, 403.141, 358.912, 389.517, 281.897, 249.391, 298.966, 306.174, 369.411, 335.71, -265.815, -264.349, -263.298, -207.83, -233.177, -180.921, -242.609, -264.709, -236.808, -202.897, -173.87, -213.984, 62.402, 40.284, 47.11, 74.843, 59.535, 76.197, 107.321, 126.359, 80.917, 126.402, 149.139, 131.007, -2.691}};
	
	double r2_z_list[4][37] = {{-224.934, -216.606, -223.442, -219.306, -214.427, -224.308, -196.288, -213.17, -200.41, -201.207, -191.506, -207.772, 102.649, 100.846, 80.92, 103.86, 96.353, 110.071, 121.464, 91.551, 111.815, 99.328, 134.976, 117.546, -39.277, -27.296, -37.426, -34.856, -34.618, -50.681, -26.316, -44.812, -55.722, -31.712, -20.25, -49.172, -4.03398},
		{-283.463, -298.562, -302.04, -299.476, -296.374, -315.801, -264.914, -275.426, -269.054, -250.743, -261.842, -258.699, 143.102, 135.77, 156.306, 134.49, 151.814, 144.838, 148.694, 177.264, 161.564, 174.057, 169.779, 160.466, -57.936, -28.731, -56.099, -47.085, -41.192, -55.114, -54.644, -58.615, -69.439, -53.235, -50.588, -47.953, -4.03398},
		{-500.664, -489.085, -497.53, -497.024, -491.502, -493.519, -445.702, -427.398, -446.01, -461.994, -427.16, -437.643, 221.633, 244.317, 234.522, 265.506, 242.061, 271.715, 274.411, 323.074, 283.98, 333.264, 365.503, 312.523, -124.731, -102.543, -81.339, -110.925, -88.522, -95.54, -126.044, -137.44, -135.143, -113.05, -109.421, -76.096, -4.03398},
		{-697.799, -768.422, -732.703, -753.391, -771.643, -762.723, -656.186, -679.921, -724.348, -654.174, -642.033, -618.106, 380.102, 374.25, 356.301, 390.269, 417.136, 442.437, 435.847, 448.618, 461.855, 483.044, 494.644, 480.004, -112.863, -149.321, -136.851, -116.849, -140.075, -109.712, -105.847, -136.603, -140.22, -181.263, -79.641, -82.917, -4.03398}};
	
	double r3_z_list[4][37] = {{92.917, 95.974, 85.962, 91.201, 97.207, 85.247, 89.436, 70.949, 80.561, 80.897, 84.835, 72.473, -28.383, -29.585, -55.637, -33.723, -48.933, -36.372, -32.108, -65.053, -48.847, -61.797, -36.171, -50.81, 26.819, 32.693, 17.953, 16.838, 10.96, -9.588, 38.998, 15.754, 1.937, 16.275, 22.138, -12.588, -0.711833},
		{165.164, 147.727, 153.768, 149.402, 141.032, 113.95, 135.439, 122.703, 126.394, 141.726, 122.33, 123.806, -51.757, -67.352, -48.494, -72.89, -60.223, -74.972, -90.497, -60.253, -87.12, -73.173, -89.004, -97.519, 39.562, 62.427, 24.079, 29.737, 25.469, 1.862, 42.175, 31.152, 15.957, 20.418, 13.053, 13.083, -0.711833},
		{284.277, 285.488, 287.627, 274.387, 282.775, 269.099, 246.417, 274.183, 264.943, 240.068, 257.962, 234.2, -156.219, -139.428, -161.336, -143.899, -167.129, -153.305, -174.617, -136.665, -184.451, -140.84, -127.788, -175.231, 55.426, 66.79, 71.687, 34.346, 32.214, 17.445, 46.654, 35.91, 27.575, 27.618, 17.186, 36.061, -0.711833},
		{418.438, 355.772, 349.288, 374.911, 341.979, 328.242, 321.96, 285.434, 274.612, 318.172, 344.265, 340.11, -177.046, -187.313, -229.768, -184.382, -203.622, -157.131, -204.76, -190.393, -205.025, -191.65, -172.309, -220.803, 132.975, 80.933, 87.059, 78.808, 41.695, 34.962, 160.549, 111.942, 89.791, 21.05, 103.938, 70.415, -0.711833}};
	 
	double std_list[4][36] = {{867.167, 880.263, 888.86, 874.585, 876.224, 855.317, 853.715, 871.423, 879.88, 875.802, 862.114, 839.255, 909.378, 916.313, 979.415, 951.352, 907.235, 893.843, 913.828, 930.735, 963.857, 971.757, 952.0, 949.319, 809.99, 819.438, 805.775, 809.79, 829.962, 819.193, 849.644, 855.646, 843.238, 849.339, 855.564, 857.103},
		{976.555, 988.253, 982.436, 983.015, 993.242, 972.815, 962.621, 967.637, 970.766, 961.01, 959.687, 934.643, 958.652, 959.844, 969.651, 980.697, 950.587, 940.724, 973.212, 973.853, 988.009, 1001.173, 978.429, 996.365, 889.388, 868.027, 885.791, 888.341, 895.088, 903.408, 921.084, 911.139, 929.352, 928.399, 928.766, 928.496},
		{1185.64, 1155.565, 1148.124, 1132.392, 1134.627, 1124.644, 1169.668, 1135.294, 1126.233, 1111.944, 1116.786, 1112.242, 1019.245, 1024.129, 1009.509, 1001.275, 1022.527, 1031.493, 1042.552, 1059.005, 1030.908, 1035.614, 1059.558, 1076.572, 962.036, 972.379, 962.586, 976.886, 1003.19, 988.396, 976.651, 1036.189, 981.144, 991.046, 1032.591, 992.296},
		{1270.448, 1252.078, 1257.065, 1236.989, 1236.054, 1245.067, 1251.264, 1236.132, 1227.845, 1226.259, 1223.823, 1229.722, 1098.736, 1108.079, 1092.661, 1095.931, 1112.502, 1099.792, 1149.622, 1153.55, 1147.899, 1146.946, 1173.448, 1156.003, 1011.36, 1038.371, 1052.357, 1059.083, 1084.751, 1074.506, 1053.898, 1092.856, 1067.656, 1075.843, 1089.171, 1077.991}};
	
	// {0,8}, {8, 16}, {16, 28}, {28, 90}
	double sector_1_array[4][37] = {{862.6, 855.342, 910.22, 897.588, 885.541, 862.528, 886.555, 861.584, 838.279, 782.628, 775.963, 773.584, -293.254, -348.62, -365.715, -384.868, -424.085, -442.6, -460.007, -492.778, -570.365, -426.941, -476.993, -486.92, 132.907, 228.429, 69.525, 252.258, 55.168, 148.501, 85.424, 200.526, 101.004, 236.928, 131.025, 249.138, 5.462},
		{885.877, 883.522, 908.855, 958.622, 931.014, 897.326, 860.224, 870.677, 868.103, 865.182, 857.594, 862.696, -357.13, -346.809, -425.034, -375.885, -434.389, -466.135, -521.797, -487.08, -603.523, -500.625, -548.981, -558.018, 157.101, 264.905, 98.743, 235.157, 113.5, 159.062, 153.352, 254.075, 98.784, 215.803, 100.827, 198.552, 5.462},
		{787.773, 806.208, 830.813, 853.18, 873.871, 902.027, 726.345, 787.947, 800.032, 850.723, 853.826, 967.267, -408.537, -365.196, -374.938, -433.476, -419.797, -369.189, -574.597, -561.766, -528.27, -550.366, -562.422, -500.919, 75.825, 247.466, 91.538, 143.356, 18.099, 215.837, 58.418, 193.531, 120.542, 251.134, 238.874, 292.452, 5.462},
		{823.064, 662.296, 797.065, 696.437, 636.029, 639.968, 803.896, 720.908, 697.749, 711.962, 778.349, 646.286, -263.828, -305.857, -358.083, -328.862, -380.462, -441.147, -474.137, -477.395, -524.717, -497.282, -448.314, -497.063, 359.617, 428.741, 152.111, 154.417, -158.603, -114.583, 309.878, 327.119, 111.467, 164.424, -20.581, -52.918, 5.462}};
	double sector_2_array[4][37] = {{1873.667, 1844.226, 1872.63, 1794.236, 1763.495, 1784.168, 1565.937, 1546.769, 1567.328, 1497.077, 1492.284, 1482.599, -686.211, -693.556, -706.316, -738.561, -798.948, -824.844, -874.788, -882.567, -954.811, -847.264, -918.467, -896.33, 399.547, 420.778, 225.101, 333.3, 137.57, 209.478, 293.222, 325.751, 237.456, 324.43, 201.718, 292.137, 19.245},
		{2096.121, 2061.84, 2041.047, 1845.949, 1744.187, 1849.915, 1591.623, 1561.102, 1585.463, 1580.668, 1550.717, 1546.78, -695.16, -701.159, -771.842, -769.122, -851.408, -860.798, -931.803, -899.03, -1038.299, -921.164, -1029.362, -1015.349, 362.523, 380.252, 231.113, 335.138, 200.578, 235.613, 339.581, 368.651, 227.537, 380.262, 235.512, 295.704, 19.245},
		{1697.317, 1745.762, 1710.123, 1776.783, 1750.397, 1771.954, 1500.529, 1547.441, 1552.337, 1594.477, 1507.121, 1591.313, -926.596, -902.854, -944.598, -1001.333, -1008.805, -976.716, -1037.382, -989.585, -1010.091, -1052.927, -1039.603, -1021.228, 340.657, 460.395, 273.196, 377.823, 166.214, 286.8, 237.929, 392.205, 265.217, 393.973, 288.385, 390.459, 19.245},
		{1839.525, 1749.051, 1747.695, 1612.714, 1671.674, 1501.324, 1595.113, 1640.0, 1567.26, 1618.027, 1608.415, 1547.562, -780.34, -836.83, -944.467, -852.586, -1005.445, -1052.944, -890.628, -912.156, -1064.936, -927.489, -1060.153, -1149.025, 586.724, 691.976, 344.371, 412.081, 63.793, 114.422, 420.346, 481.692, 332.862, 355.761, 167.852, 81.941, 19.245}};
	double sector_3_array[4][37] = {{231.36, 212.983, 241.335, 258.747, 236.619, 282.258, 205.339, 191.409, 185.142, 179.512, 194.836, 229.217, -110.87, -114.333, -167.118, -134.997, -159.371, -154.312, -144.95, -158.443, -225.034, -28.553, -64.373, -77.593, 79.568, 138.545, -30.224, 154.02, -29.131, 67.069, -13.136, 53.462, -32.276, 100.052, -11.043, 88.653, -14.711},
		{288.124, 310.583, 296.182, 350.924, 360.798, 367.444, 260.0, 254.498, 243.318, 310.48, 289.973, 333.488, -133.051, -166.404, -192.166, -138.08, -191.424, -145.889, -201.175, -169.123, -241.564, -63.005, -198.46, -164.22, 70.633, 123.556, -1.657, 128.848, -14.683, 59.98, 18.626, 94.351, -38.165, 100.029, 11.486, 96.414, -14.711},
		{260.798, 318.689, 328.13, 370.613, 362.692, 425.505, 177.723, 255.324, 272.934, 342.022, 339.097, 434.18, -277.992, -237.93, -184.025, -220.588, -231.049, -79.158, -266.698, -163.044, -187.872, -152.801, -164.564, -78.319, -38.704, 183.905, 3.779, 92.584, -15.128, 158.354, -81.175, 28.345, 3.899, 133.49, 75.89, 196.386, -14.711},
		{381.49, 245.795, 255.758, 220.284, 200.405, 147.159, 207.983, 160.371, 61.41, 157.146, 109.095, 119.833, -133.387, -172.308, -225.703, -77.115, -199.151, -254.495, -18.471, -105.418, -136.934, -67.325, -120.4, -262.24, 260.69, 389.453, 56.33, 74.995, -223.777, -141.586, 170.75, 174.188, -83.666, 111.246, -181.11, -86.299, -14.711}};
	double sector_4_array[4][37] = {{809.856, 854.794, 883.45, 830.16, 826.329, 841.603, 760.076, 758.083, 752.929, 733.337, 743.293, 773.725, -409.225, -400.995, -392.736, -336.104, -346.202, -349.043, -449.739, -465.05, -503.344, -365.155, -409.145, -403.27, 94.167, 216.006, 67.526, 221.749, 86.204, 195.391, 53.881, 127.417, 45.038, 206.543, 103.717, 230.805, -18.594},
		{831.92, 844.587, 856.007, 883.224, 876.667, 894.295, 727.869, 743.082, 750.988, 774.355, 768.836, 782.043, -434.102, -413.42, -444.826, -387.744, -416.718, -381.642, -456.678, -406.235, -530.955, -372.751, -450.027, -439.408, 100.273, 179.405, 61.864, 268.82, 102.814, 200.703, 84.81, 163.311, 65.46, 213.562, 130.95, 206.492, -18.594},
		{887.791, 939.687, 942.85, 967.339, 963.0, 1061.577, 720.121, 857.885, 857.538, 936.649, 946.212, 1016.816, -568.236, -507.309, -500.5, -541.094, -510.051, -435.415, -657.304, -587.173, -598.55, -607.559, -603.821, -596.329, 173.735, 315.837, 157.406, 262.972, 76.718, 264.176, 70.15, 180.376, 111.203, 288.766, 269.511, 354.049, -18.594},
		{1165.861, 1134.373, 1128.773, 1059.206, 1003.297, 962.332, 942.079, 867.21, 892.463, 902.015, 907.013, 802.331, -429.531, -521.67, -598.896, -550.603, -559.165, -636.294, -542.27, -570.323, -612.357, -554.468, -644.71, -650.467, 348.752, 418.415, 126.905, 310.951, -80.314, 32.042, 266.189, 273.597, 181.222, 418.701, 102.063, 119.344, -18.594}};
	double sector_5_array[4][37] = {{-1215.051, -1147.962, -1156.989, -1160.146, -1157.742, -1172.134, -1064.775, -1021.591, -1037.024, -1067.182, -1040.073, -1069.433, 563.092, 578.339, 246.218, 530.125, 556.987, 553.293, 596.874, 553.292, 490.766, 671.207, 619.026, 625.633, -213.011, -128.27, -279.007, -94.856, -236.873, -125.127, -291.838, -180.596, -259.86, -107.702, -181.508, -88.533, -32.732},
		{-1047.348, -1025.573, -1062.716, -1010.348, -1012.486, -1005.55, -894.864, -887.829, -895.572, -854.535, -870.258, -850.336, 486.114, 524.481, 466.888, 519.311, 530.292, 529.262, 533.603, 558.704, 459.664, 608.392, 514.588, 522.644, -217.95, -133.361, -253.444, -92.854, -237.742, -144.884, -236.089, -149.231, -215.409, -87.384, -107.696, -82.075, -32.732},
		{-973.594, -961.777, -962.01, -903.112, -892.808, -819.765, -811.527, -784.93, -762.528, -693.205, -698.904, -640.54, 399.988, 501.751, 515.707, 534.79, 510.086, 597.582, 437.637, 479.479, 502.035, 427.681, 433.837, 467.496, -299.019, -77.313, -246.175, -161.585, -274.395, -77.569, -305.538, -178.837, -229.242, -79.163, -97.898, 18.613, -32.732},
		{-982.628, -1030.053, -965.898, -988.525, -970.148, -1091.929, -649.514, -758.388, -844.496, -873.838, -791.843, -823.106, 588.729, 526.081, 443.59, 543.949, 443.683, 382.621, 591.79, 564.595, 516.058, 532.004, 573.244, 414.27, -43.553, 13.922, -255.325, -154.16, -491.578, -337.128, -30.541, 56.232, -115.477, -70.609, -164.144, -136.38, -32.732}};
	double sector_6_array[4][37] = {{-427.561, -409.68, -396.912, -400.34, -383.518, -392.765, -323.127, -340.963, -363.986, -429.454, -415.239, -396.198, 179.241, 150.377, 170.922, 160.737, 167.546, 168.845, 192.682, 177.383, 148.277, 255.897, 238.807, 279.357, -80.53, -0.776, -101.984, 23.987, -112.224, -18.37, -184.203, -70.483, -159.915, -12.557, -87.854, 40.054, -14.989},
		{-247.442, -245.478, -253.695, -201.191, -208.663, -232.025, -189.108, -199.347, -221.814, -201.032, -218.448, -195.221, 106.102, 126.382, 91.106, 75.266, 138.811, 156.993, 84.213, 115.966, 50.824, 186.988, 95.966, 143.266, -41.403, 19.713, -109.152, 43.851, -81.828, 4.242, -116.31, -37.011, -118.102, 46.702, -36.036, 66.951, -14.989},
		{-406.414, -373.312, -369.928, -340.374, -338.826, -249.943, -344.715, -292.223, -242.484, -224.972, -230.391, -145.458, 68.066, 123.913, 199.391, 133.727, 182.08, 229.643, 78.367, 168.964, 154.53, 132.968, 175.075, 191.424, -117.604, 3.357, -78.511, -7.727, -135.173, 65.409, -226.054, -80.354, -108.034, -3.693, 1.288, 141.666, -14.989},
		{-495.762, -555.07, -536.828, -569.286, -589.709, -614.637, -426.423, -560.098, -563.852, -583.757, -594.78, -640.328, 313.749, 218.749, 152.097, 242.523, 212.532, 135.282, 450.234, 330.282, 310.113, 398.064, 299.94, 236.778, 140.015, 137.276, -135.521, -101.159, -378.931, -282.567, 15.048, 12.77, -124.22, -87.245, -225.211, -180.094, -14.989}};

	double current_position;
	if (sector == 1) {
		current_position = sector_1_array[j][i];
	} else if (sector == 2)	{
		current_position = sector_2_array[j][i];
	} else if (sector == 3) {
		current_position = sector_3_array[j][i];
	} else if (sector == 4) {
		current_position = sector_4_array[j][i];
	} else if (sector == 5) {
		current_position = sector_5_array[j][i];
	} else if (sector == 6) {
		current_position = sector_6_array[j][i];
	}

	// desired vertex position is actually the position of the downstream target wall
	// double desired_positions[38] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.72891, 5.69459};
	double desired_positions[37] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -5};

	double value = ((1/shift)*(par[0]*r1_x_list[j][i]+par[1]*r1_z_list[j][i]+
		par[2]*r2_x_list[j][i]+par[3]*r2_z_list[j][i]+par[4]*r3_x_list[j][i]+
		par[5]*r3_z_list[j][i]))-current_position+desired_positions[i];
	
	value = pow(value,2) / pow(std_list[j][i],2);

	return value;
}

void chi2(int &npar, double *gin, double &f, double *par, int iflag) {
	//calculate chisquare
	double chisq = 0;
	// loop over all layers and vertex position (36 layers + 1 vertex)
	for (int j = 0; j < 4; ++j) { 
		for (int i = 0; i < 37; ++i) {
			// each value weighted equally, so no division by uncertainty, maybe adjust later?
			chisq += calculate_residual(cur_sector, j, i, par) ;
			// cout << i << endl;
		}
	}
	f = chisq; 
}

void minimizer_009() {

	for (int sector = 1; sector < 7; ++sector) {
		cur_sector = sector;
		TMinuit *gMinuit = new TMinuit(6);  //initialize TMinuit with a maximum of 12 params
		gMinuit->SetPrintLevel(-1);
		gMinuit->SetFCN(chi2); // chi^2 function to be minimized
		double arglist[10]; arglist[0] = 1;
		int ierflg = 0;
		gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
		// Set starting values and step sizes for parameters
		// (param int, name, start value, step size, min limit, upper limit, error flag)
		gMinuit->mnparm(0, "r1x", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(1, "r1z", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(2, "r2x", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(3, "r2z", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(4, "r3x", 0, 0, 0, 0, ierflg);
		gMinuit->mnparm(5, "r3z", 0, 0.1, 0, 0, ierflg);

		arglist[1] = 1.0; arglist[0] = 50000;
		gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

		gMinuit->GetParameter(0,r1x[sector-1],error);
		gMinuit->GetParameter(1,r1z[sector-1],error);
		gMinuit->GetParameter(2,r2x[sector-1],error);
		gMinuit->GetParameter(3,r2z[sector-1],error);
		gMinuit->GetParameter(4,r3x[sector-1],error);
		gMinuit->GetParameter(5,r3z[sector-1],error);

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
	cout<<"1 1 0 "<<r1xGlobal[0]<<" "<<r1yGlobal[0]<<" "<<r1z[0]<<" 0 "<<0<<" 0"<<endl;
	cout<<"1 2 0 "<<r1xGlobal[1]<<" "<<r1yGlobal[1]<<" "<<r1z[1]<<" 0 "<<0<<" 0"<<endl;
	cout<<"1 3 0 "<<r1xGlobal[2]<<" "<<r1yGlobal[2]<<" "<<r1z[2]<<" 0 "<<0<<" 0"<<endl;
	cout<<"1 4 0 "<<r1xGlobal[3]<<" "<<r1yGlobal[3]<<" "<<r1z[3]<<" 0 "<<0<<" 0"<<endl;
	cout<<"1 5 0 "<<r1xGlobal[4]<<" "<<r1yGlobal[4]<<" "<<r1z[4]<<" 0 "<<0<<" 0"<<endl;
	cout<<"1 6 0 "<<r1xGlobal[5]<<" "<<r1yGlobal[5]<<" "<<r1z[5]<<" 0 "<<0<<" 0"<<endl;
	cout<<"2 1 0 "<<r2xGlobal[0]<<" "<<r2yGlobal[0]<<" "<<r2z[0]<<" 0 "<<0<<" 0"<<endl;
	cout<<"2 2 0 "<<r2xGlobal[1]<<" "<<r2yGlobal[1]<<" "<<r2z[1]<<" 0 "<<0<<" 0"<<endl;
	cout<<"2 3 0 "<<r2xGlobal[2]<<" "<<r2yGlobal[2]<<" "<<r2z[2]<<" 0 "<<0<<" 0"<<endl;
	cout<<"2 4 0 "<<r2xGlobal[3]<<" "<<r2yGlobal[3]<<" "<<r2z[3]<<" 0 "<<0<<" 0"<<endl;
	cout<<"2 5 0 "<<r2xGlobal[4]<<" "<<r2yGlobal[4]<<" "<<r2z[4]<<" 0 "<<0<<" 0"<<endl;
	cout<<"2 6 0 "<<r2xGlobal[5]<<" "<<r2yGlobal[5]<<" "<<r2z[5]<<" 0 "<<0<<" 0"<<endl;
	cout<<"3 1 0 "<<r3xGlobal[0]<<" "<<r3yGlobal[0]<<" "<<r3z[0]<<" 0 "<<0<<" 0"<<endl;
	cout<<"3 2 0 "<<r3xGlobal[1]<<" "<<r3yGlobal[1]<<" "<<r3z[1]<<" 0 "<<0<<" 0"<<endl;
	cout<<"3 3 0 "<<r3xGlobal[2]<<" "<<r3yGlobal[2]<<" "<<r3z[2]<<" 0 "<<0<<" 0"<<endl;
	cout<<"3 4 0 "<<r3xGlobal[3]<<" "<<r3yGlobal[3]<<" "<<r3z[3]<<" 0 "<<0<<" 0"<<endl;
	cout<<"3 5 0 "<<r3xGlobal[4]<<" "<<r3yGlobal[4]<<" "<<r3z[4]<<" 0 "<<0<<" 0"<<endl;
	cout<<"3 6 0 "<<r3xGlobal[5]<<" "<<r3yGlobal[5]<<" "<<r3z[5]<<" 0 "<<0<<" 0"<<endl;
}