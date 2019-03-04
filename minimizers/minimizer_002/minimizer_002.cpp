#include "Riostream.h"
#include <TROOT.h>
#include "TMinuit.h"
#include "TMath.h"
#include "TRandom1.h"
#include <ctime>
#include <algorithm>    // std::min_element, std::max_element

std::clock_t start = std::clock();
int s = (std::clock() - start) / (double) (CLOCKS_PER_SEC );
static const int num_residuals = 36*6;
float run_2467[num_residuals]; 
float r1_x[num_residuals]; float r1_y[num_residuals]; float r1_z[num_residuals]; 
float r1_cx[num_residuals];
float r2_x[num_residuals]; float r2_y[num_residuals]; float r2_z[num_residuals]; 
float r2_cx[num_residuals];
float r3_x[num_residuals]; float r3_y[num_residuals]; float r3_z[num_residuals]; 
float r3_cx[num_residuals];

float r1_x_diff[num_residuals]; float r1_y_diff[num_residuals]; float r1_z_diff[num_residuals]; 
float r1_cx_diff[num_residuals];
float r2_x_diff[num_residuals]; float r2_y_diff[num_residuals]; float r2_z_diff[num_residuals]; 
float r2_cx_diff[num_residuals];
float r3_x_diff[num_residuals]; float r3_y_diff[num_residuals]; float r3_z_diff[num_residuals]; 
float r3_cx_diff[num_residuals];

// main program area
void minimizer_002() {

	// Read in number of data points from external file
	ifstream dataTable;
	string firstPart = "/Users/tbhayward/Dropbox/Personal_Files/Academics/2018/Research/";
	string secondPart = "Drift_Chambers/minimizers/minimizer_002/Residuals/";
	string thirdPart = "run_2467_copy";
	string fourthPart = ".txt";
	string file = firstPart+secondPart+thirdPart+fourthPart;
	dataTable.open( file );
	for (int i = 0; i < num_residuals; ++i) {
		dataTable >> run_2467[i];
	}
	dataTable.close();

	thirdPart = "r1_x";
	file = firstPart+secondPart+thirdPart+fourthPart;
	dataTable.open( file );
	for (int i = 0; i < num_residuals; ++i) {
		dataTable >> r1_x[i];
	}
	dataTable.close();

	thirdPart = "r1_y";
	file = firstPart+secondPart+thirdPart+fourthPart;
	dataTable.open( file );
	for (int i = 0; i < num_residuals; ++i) {
		dataTable >> r1_y[i];
	}
	dataTable.close();

	thirdPart = "r1_z";
	file = firstPart+secondPart+thirdPart+fourthPart;
	dataTable.open( file );
	for (int i = 0; i < num_residuals; ++i) {
		dataTable >> r1_z[i];
	}
	dataTable.close();

	thirdPart = "r1_cx";
	file = firstPart+secondPart+thirdPart+fourthPart;
	dataTable.open( file );
	for (int i = 0; i < num_residuals; ++i) {
		dataTable >> r1_cx[i];
	}
	dataTable.close();

	thirdPart = "r2_x";
	file = firstPart+secondPart+thirdPart+fourthPart;
	dataTable.open( file );
	for (int i = 0; i < num_residuals; ++i) {
		dataTable >> r2_x[i];
	}
	dataTable.close();

	thirdPart = "r2_y";
	file = firstPart+secondPart+thirdPart+fourthPart;
	dataTable.open( file );
	for (int i = 0; i < num_residuals; ++i) {
		dataTable >> r2_y[i];
	}
	dataTable.close();

	thirdPart = "r2_z";
	file = firstPart+secondPart+thirdPart+fourthPart;
	dataTable.open( file );
	for (int i = 0; i < num_residuals; ++i) {
		dataTable >> r2_z[i];
	}
	dataTable.close();

	thirdPart = "r2_cx";
	file = firstPart+secondPart+thirdPart+fourthPart;
	dataTable.open( file );
	for (int i = 0; i < num_residuals; ++i) {
		dataTable >> r2_cx[i];
	}
	dataTable.close();

	thirdPart = "r3_x";
	file = firstPart+secondPart+thirdPart+fourthPart;
	dataTable.open( file );
	for (int i = 0; i < num_residuals; ++i) {
		dataTable >> r3_x[i];
	}
	dataTable.close();

	thirdPart = "r3_y";
	file = firstPart+secondPart+thirdPart+fourthPart;
	dataTable.open( file );
	for (int i = 0; i < num_residuals; ++i) {
		dataTable >> r3_y[i];
	}
	dataTable.close();

	thirdPart = "r3_z";
	file = firstPart+secondPart+thirdPart+fourthPart;
	dataTable.open( file );
	for (int i = 0; i < num_residuals; ++i) {
		dataTable >> r3_z[i];
	}
	dataTable.close();

	thirdPart = "r3_cx";
	file = firstPart+secondPart+thirdPart+fourthPart;
	dataTable.open( file );
	for (int i = 0; i < num_residuals; ++i) {
		dataTable >> r3_cx[i];
	}
	dataTable.close();

	static const int num_scales = 60000; double scales[num_scales]; 
	for (int i = 0; i < num_scales; ++i) {
		double j = i;
		scales[i] = 8*j/(num_scales-1) - 4;
		// cout << scales[i] << endl;
	}
	double run_2467_chi2[num_scales] = {0}; 
	double r1_x_chi2[num_scales] = {0}; 
	double r1_y_chi2[num_scales] ={0}; 
	double r1_z_chi2[num_scales] = {0}; 
	double r1_cx_chi2[num_scales] = {0};
	double r2_x_chi2[num_scales] = {0}; 
	double r2_y_chi2[num_scales] = {0}; 
	double r2_z_chi2[num_scales] = {0}; 
	double r2_cx_chi2[num_scales] = {0};
	double r3_x_chi2[num_scales] = {0}; 
	double r3_y_chi2[num_scales] = {0}; 
	double r3_z_chi2[num_scales] = {0}; 
	double r3_cx_chi2[num_scales] = {0}; 


	for (int j = 0; j < num_scales; ++j) {
		for (int i = 0; i < num_residuals; ++i) {
			run_2467_chi2[j] += TMath::Abs(run_2467[i]);

			r1_x_diff[i] = run_2467[i] - scales[j]*r1_x[i]; 
			r1_x_chi2[j] += TMath::Abs(r1_x_diff[i]); 

			r1_y_diff[i] = run_2467[i] - scales[j]*r1_y[i]; 
			r1_y_chi2[j] += TMath::Abs(r1_y_diff[i]);

			r1_z_diff[i] = run_2467[i] - scales[j]*r1_z[i]; 
			r1_z_chi2[j] += TMath::Abs(r1_z_diff[i]);

			r1_cx_diff[i] = run_2467[i] - scales[j]*r1_cx[i]; 
			r1_cx_chi2[j] += TMath::Abs(r1_cx_diff[i]);

			r2_x_diff[i] = run_2467[i] - scales[j]*r2_x[i]; 
			r2_x_chi2[j] += TMath::Abs(r2_x_diff[i]);

			r2_y_diff[i] = run_2467[i] - scales[j]*r2_y[i]; 
			r2_y_chi2[j] += TMath::Abs(r2_y_diff[i]);

			r2_z_diff[i] = run_2467[i] - scales[j]*r2_z[i]; 
			r2_z_chi2[j] += TMath::Abs(r2_z_diff[i]);

			r2_cx_diff[i] = run_2467[i] - scales[j]*r2_cx[i]; 
			r2_cx_chi2[j] += TMath::Abs(r2_cx_diff[i]);

			r3_x_diff[i] = run_2467[i] - scales[j]*r3_x[i]; 
			r3_x_chi2[j] += TMath::Abs(r3_x_diff[i]);

			r3_y_diff[i] = run_2467[i] - scales[j]*r3_y[i]; 
			r3_y_chi2[j] += TMath::Abs(r3_y_diff[i]);

			r3_z_diff[i] = run_2467[i] - scales[j]*r3_z[i]; 
			r3_z_chi2[j] += TMath::Abs(r3_z_diff[i]);

			r3_cx_diff[i] = run_2467[i] - scales[j]*r3_cx[i]; 
			r3_cx_chi2[j] += TMath::Abs(r3_cx_diff[i]);
		}
	}

	int min_r1_z_index = 0;
	for (int i = 1; i < num_scales; ++i) {
		if (r1_z_chi2[i]<r1_z_chi2[min_r1_z_index]) {
			min_r1_z_index = i;
		}
	}
	cout<<min_r1_z_index<<" "<<scales[min_r1_z_index]<<" "<<r1_z_chi2[min_r1_z_index]/
		run_2467_chi2[0]<<endl;

	int min_r3_z_index = 0;
	for (int i = 1; i < num_scales; ++i) {
		if (r3_z_chi2[i]<r3_z_chi2[min_r3_z_index]) {
			min_r3_z_index = i;
		}
	}
	cout<<min_r3_z_index<<" "<<scales[min_r3_z_index]<<" "<<r3_z_chi2[min_r3_z_index]/
		run_2467_chi2[0]<<endl;

	int min_r1_cx_index = 0;
	for (int i = num_scales/2; i < num_scales; ++i) {
		if (r1_cx_chi2[i]<r1_cx_chi2[min_r1_cx_index]) {
			min_r1_cx_index = i;
		}
	}
	cout<<min_r1_cx_index<<" "<<scales[min_r1_cx_index]<<" "<<r1_cx_chi2[min_r1_cx_index]/
		run_2467_chi2[0]<<endl;

	int min_r2_cx_index = 0;
	for (int i = num_scales/2; i < num_scales; ++i) {
		if (r2_cx_chi2[i]<r2_cx_chi2[min_r2_cx_index]) {
			min_r2_cx_index = i;
		}
	}
	cout<<min_r2_cx_index<<" "<<scales[min_r2_cx_index]<<" "<<r2_cx_chi2[min_r2_cx_index]/
		run_2467_chi2[0]<<endl;

	int min_r3_cx_index = 0;
	for (int i = num_scales/2; i < num_scales; ++i) {
		if (r3_cx_chi2[i]<r3_cx_chi2[min_r3_cx_index]) {
			min_r3_cx_index = i;
		}
	}
	cout<<min_r3_cx_index<<" "<<scales[min_r3_cx_index]<<" "<<r3_cx_chi2[min_r3_cx_index]/
		run_2467_chi2[0]<<endl;




}