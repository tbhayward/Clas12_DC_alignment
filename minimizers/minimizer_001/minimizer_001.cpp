/*
root function designed to find the optimal maximum Q2 and function to use for extracting
the proton radius
*/

#include "Riostream.h"
#include <TROOT.h>
#include "TMinuit.h"
#include "TMath.h"
#include "TRandom1.h"
#include <ctime>

std::clock_t start = std::clock();
int s = (std::clock() - start) / (double) (CLOCKS_PER_SEC );

float r1_x[18],r1_y[18],r1_z[18],r1_cx[18],r1_cy[18],r1_cz[18],
	r2_x[18],r2_y[18],r2_z[18],r2_cx[18],r2_cy[18],r2_cz[18],
	r3_x[18],r3_y[18],r3_z[18],r3_cx[18],r3_cy[18],r3_cz[18], data[18];

double fnc(float r1_x, float r1_y, float r1_z, float r1_cx, float r1_cy, float r1_cz,
	float r2_x, float r2_y, float r2_z, float r2_cx, float r2_cy, float r2_cz,
	float r3_x, float r3_y, float r3_z, float r3_cx, float r3_cy, float r3_cz, double *par) {

	double value; // value to be returned by evaluated function
	value = (par[0]+r1_x)+(par[1]+r1_y)+(par[2]+r1_z)+(par[3]+r1_cx)+(par[4]+r1_cy)+(par[5]+r1_cz)+
		(par[6]+r2_x)+(par[7]+r2_y)+(par[8]+r2_z)+(par[9]+r2_cx)+(par[10]+r2_cy)+(par[11]+r2_cz)+
		(par[12]+r3_x)+(par[13]+r3_y)+(par[14]+r3_z)+(par[15]+r3_cx)+(par[16]+r3_cy)+
		(par[17]+r3_cz);
	return value;
}

void fit_minimizer(int &npar, double *gin, double &f, double *par, int iflag) {
	//calculate chisquare
	double chisq = 0;
	double delta;
	for (int i=0;i<18; i++) {
    	delta  = (data[i]-fnc(r1_x[i], r1_y[i], r1_z[i], r1_cx[i], r1_cy[i], r1_cz[i],
    		r2_x[i], r2_y[i], r2_z[i], r2_cx[i], r2_cy[i], r2_cz[i],
    		r3_x[i], r3_y[i], r3_z[i], r3_cx[i], r3_cy[i], r3_cz[i], par))/1;
    	chisq += delta*delta;
    }
    f = chisq;
}

void minimizer_001() {
	r1_x[0]=1;r1_y[0]=0;r1_z[0]=0;r1_cx[0]=0;r1_cy[0]=0;r1_cz[0]=0;
	r2_x[0]=0;r2_y[0]=0;r2_z[0]=0;r2_cx[0]=0;r2_cy[0]=0;r2_cz[0]=0;
	r3_x[0]=0;r3_y[0]=0;r3_z[0]=0;r3_cx[0]=0;r3_cy[0]=0;r3_cz[0]=0;

	r1_x[1]=0;r1_y[1]=1;r1_z[1]=0;r1_cx[1]=0;r1_cy[1]=0;r1_cz[1]=0;
	r2_x[1]=0;r2_y[1]=0;r2_z[1]=0;r2_cx[1]=0;r2_cy[1]=0;r2_cz[1]=0;
	r3_x[1]=0;r3_y[1]=0;r3_z[1]=0;r3_cx[1]=0;r3_cy[1]=0;r3_cz[1]=0;

	r1_x[2]=0;r1_y[2]=0;r1_z[2]=1;r1_cx[2]=0;r1_cy[2]=0;r1_cz[2]=0;
	r2_x[2]=0;r2_y[2]=0;r2_z[2]=0;r2_cx[2]=0;r2_cy[2]=0;r2_cz[2]=0;
	r3_x[2]=0;r3_y[2]=0;r3_z[2]=0;r3_cx[2]=0;r3_cy[2]=0;r3_cz[2]=0;

	r1_x[3]=0;r1_y[3]=0;r1_z[3]=0;r1_cx[3]=0.5;r1_cy[3]=0;r1_cz[3]=0;
	r2_x[3]=0;r2_y[3]=0;r2_z[3]=0;r2_cx[3]=0;r2_cy[3]=0;r2_cz[3]=0;
	r3_x[3]=0;r3_y[3]=0;r3_z[3]=0;r3_cx[3]=0;r3_cy[3]=0;r3_cz[3]=0;

	r1_x[4]=0;r1_y[4]=0;r1_z[4]=0;r1_cx[4]=0;r1_cy[4]=0.5;r1_cz[4]=0;
	r2_x[4]=0;r2_y[4]=0;r2_z[4]=0;r2_cx[4]=0;r2_cy[4]=0;r2_cz[4]=0;
	r3_x[4]=0;r3_y[4]=0;r3_z[4]=0;r3_cx[4]=0;r3_cy[4]=0;r3_cz[4]=0;

	r1_x[5]=0;r1_y[5]=0;r1_z[5]=0;r1_cx[5]=0;r1_cy[5]=0;r1_cz[5]=0.5;
	r2_x[5]=0;r2_y[5]=0;r2_z[5]=0;r2_cx[5]=0;r2_cy[5]=0;r2_cz[5]=0;
	r3_x[5]=0;r3_y[5]=0;r3_z[5]=0;r3_cx[5]=0;r3_cy[5]=0;r3_cz[5]=0;

	r1_x[6]=0;r1_y[6]=0;r1_z[6]=0;r1_cx[6]=0;r1_cy[6]=0;r1_cz[6]=0;
	r2_x[6]=1;r2_y[6]=0;r2_z[6]=0;r2_cx[6]=0;r2_cy[6]=0;r2_cz[6]=0;
	r3_x[6]=0;r3_y[6]=0;r3_z[6]=0;r3_cx[6]=0;r3_cy[6]=0;r3_cz[6]=0;

	r1_x[7]=0;r1_y[7]=0;r1_z[7]=0;r1_cx[7]=0;r1_cy[7]=0;r1_cz[7]=0;
	r2_x[7]=0;r2_y[7]=1;r2_z[7]=0;r2_cx[7]=0;r2_cy[7]=0;r2_cz[7]=0;
	r3_x[7]=0;r3_y[7]=0;r3_z[7]=0;r3_cx[7]=0;r3_cy[7]=0;r3_cz[7]=0;

	r1_x[8]=0;r1_y[8]=0;r1_z[8]=0;r1_cx[8]=0;r1_cy[8]=0;r1_cz[8]=0;
	r2_x[8]=0;r2_y[8]=0;r2_z[8]=1;r2_cx[8]=0;r2_cy[8]=0;r2_cz[8]=0;
	r3_x[8]=0;r3_y[8]=0;r3_z[8]=0;r3_cx[8]=0;r3_cy[8]=0;r3_cz[8]=0;

	r1_x[9]=0;r1_y[9]=0;r1_z[9]=0;r1_cx[9]=0;r1_cy[9]=0;r1_cz[9]=0;
	r2_x[9]=0;r2_y[9]=0;r2_z[9]=0;r2_cx[9]=0.5;r2_cy[9]=0;r2_cz[9]=0;
	r3_x[9]=0;r3_y[9]=0;r3_z[9]=0;r3_cx[9]=0;r3_cy[9]=0;r3_cz[9]=0;

	r1_x[10]=0;r1_y[10]=0;r1_z[10]=0;r1_cx[10]=0;r1_cy[10]=0;r1_cz[10]=0;
	r2_x[10]=0;r2_y[10]=0;r2_z[10]=0;r2_cx[10]=0;r2_cy[10]=0.5;r2_cz[10]=0;
	r3_x[10]=0;r3_y[10]=0;r3_z[10]=0;r3_cx[10]=0;r3_cy[10]=0;r3_cz[10]=0;

	r1_x[11]=0;r1_y[11]=0;r1_z[11]=0;r1_cx[11]=0;r1_cy[11]=0;r1_cz[11]=0;
	r2_x[11]=0;r2_y[11]=0;r2_z[11]=0;r2_cx[11]=0;r2_cy[11]=0;r2_cz[11]=0.5;
	r3_x[11]=0;r3_y[11]=0;r3_z[11]=0;r3_cx[11]=0;r3_cy[11]=0;r3_cz[11]=0;

	r1_x[12]=0;r1_y[12]=0;r1_z[12]=0;r1_cx[12]=0;r1_cy[12]=0;r1_cz[12]=0;
	r2_x[12]=0;r2_y[12]=0;r2_z[12]=0;r2_cx[12]=0;r2_cy[12]=0;r2_cz[12]=0;
	r3_x[12]=1;r3_y[12]=0;r3_z[12]=0;r3_cx[12]=0;r3_cy[12]=0;r3_cz[12]=0;

	r1_x[13]=0;r1_y[13]=0;r1_z[13]=0;r1_cx[13]=0;r1_cy[13]=0;r1_cz[13]=0;
	r2_x[13]=0;r2_y[13]=0;r2_z[13]=0;r2_cx[13]=0;r2_cy[13]=0;r2_cz[13]=0;
	r3_x[13]=0;r3_y[13]=1;r3_z[13]=0;r3_cx[13]=0;r3_cy[13]=0;r3_cz[13]=0;

	r1_x[14]=0;r1_y[14]=0;r1_z[14]=0;r1_cx[14]=0;r1_cy[14]=0;r1_cz[14]=0;
	r2_x[14]=0;r2_y[14]=0;r2_z[14]=0;r2_cx[14]=0;r2_cy[14]=0;r2_cz[14]=0;
	r3_x[14]=0;r3_y[14]=0;r3_z[14]=1;r3_cx[14]=0;r3_cy[14]=0;r3_cz[14]=0;

	r1_x[15]=0;r1_y[15]=0;r1_z[15]=0;r1_cx[15]=0;r1_cy[15]=0;r1_cz[15]=0;
	r2_x[15]=0;r2_y[15]=0;r2_z[15]=0;r2_cx[15]=0;r2_cy[15]=0;r2_cz[15]=0;
	r3_x[15]=0;r3_y[15]=0;r3_z[15]=0;r3_cx[15]=0.5;r3_cy[15]=0;r3_cz[15]=0;

	r1_x[16]=0;r1_y[16]=0;r1_z[16]=0;r1_cx[16]=0;r1_cy[16]=0;r1_cz[16]=0;
	r2_x[16]=0;r2_y[16]=0;r2_z[16]=0;r2_cx[16]=0;r2_cy[16]=0;r2_cz[16]=0;
	r3_x[16]=0;r3_y[16]=0;r3_z[16]=0;r3_cx[16]=0;r3_cy[16]=0.5;r3_cz[16]=0;

	r1_x[17]=0;r1_y[17]=0;r1_z[17]=0;r1_cx[17]=0;r1_cy[17]=0;r1_cz[17]=0;
	r2_x[17]=0;r2_y[17]=0;r2_z[17]=0;r2_cx[17]=0;r2_cy[17]=0;r2_cz[17]=0;
	r3_x[17]=0;r3_y[17]=0;r3_z[17]=0;r3_cx[17]=0;r3_cy[17]=0;r3_cz[17]=0.5;

	data[0] = 3957.4967041660325; data[1] = 2759.270966965731; data[2] = 2773.235004723582;
	data[3] = 2772.8131008564446; data[4] = 2741.588545132006; data[5] = 2816.3199541504896;
	data[6] = 4385.567036789882; data[7] = 2775.3726715217; data[8] = 3699.783137022425;
	data[9] = 2764.200840165788; data[10] = 2736.0524635301886; data[11] = 3266.8758851886123;
	data[12] = 3835.7094582405907; data[13] = 2752.5241146659528; data[14] = 2752.759208149603;
	data[15] = 2813.9879180668117; data[16] = 2747.254089435589; data[17] = 3010.5130546230594;

	TMinuit *gMinuit = new TMinuit(18);  //initialize TMinuit with a maximum of 5 params
	gMinuit->SetFCN(fit_minimizer);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

	// Set starting values and step sizes for parameters
	static Double_t vstart[4] = {0, 0 , 0 , 0};
	static Double_t step[4] = {0.1 , 0.1 , 0.1 , 0.1};
	gMinuit->mnparm(0, "r1_x", vstart[0], step[0], 0,0,ierflg);
	gMinuit->mnparm(1, "r1_y", vstart[1], step[1], 0,0,ierflg);
	gMinuit->mnparm(2, "r1_z", vstart[2], step[2], 0,0,ierflg);
	gMinuit->mnparm(3, "r1_cx", vstart[3], step[3], 0,0,ierflg);
	gMinuit->mnparm(4, "r1_cy", vstart[0], step[0], 0,0,ierflg);
	gMinuit->mnparm(5, "r1_cz", vstart[1], step[1], 0,0,ierflg);
	gMinuit->mnparm(6, "r2_x", vstart[0], step[0], 0,0,ierflg);
	gMinuit->mnparm(7, "r2_y", vstart[1], step[1], 0,0,ierflg);
	gMinuit->mnparm(8, "r2_z", vstart[2], step[2], 0,0,ierflg);
	gMinuit->mnparm(9, "r2_cx", vstart[3], step[3], 0,0,ierflg);
	gMinuit->mnparm(10, "r2_cy", vstart[0], step[0], 0,0,ierflg);
	gMinuit->mnparm(11, "r2_cz", vstart[1], step[1], 0,0,ierflg);
	gMinuit->mnparm(12, "r3_x", vstart[0], step[0], 0,0,ierflg);
	gMinuit->mnparm(13, "r3_y", vstart[1], step[1], 0,0,ierflg);
	gMinuit->mnparm(14, "r3_z", vstart[2], step[2], 0,0,ierflg);
	gMinuit->mnparm(15, "r3_cx", vstart[3], step[3], 0,0,ierflg);
	gMinuit->mnparm(16, "r3_cy", vstart[0], step[0], 0,0,ierflg);
	gMinuit->mnparm(17, "r3_cz", vstart[1], step[1], 0,0,ierflg);

	// Now ready for minimization step
	arglist[0] = 500;
	arglist[1] = 1.;
	gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
	gMinuit->SetPrintLevel(3);

	// Print results
	double amin,edm,errdef;
	int nvpar,nparx,icstat;
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	//gMinuit->mnprin(3,amin);
}