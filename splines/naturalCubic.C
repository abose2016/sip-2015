#include "TCanvas.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include "Math/Polynomial.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"

// How to plot an interpolator? -- the ROOT way, with a global variable.
ROOT::Math::Interpolator *gInter;

double interp_wrap(double *x, double *par) {
	double s;
	s = gInter->Eval(x[0]);
	return s;
}

//Here I define two different types of test data that we can use
enum tfun_t { trig, sort_rand };

//make some vectors and fill them with either random data or points from sinx
void fill_test_data(vector<double> &xv, vector<double> &yv, tfun_t t){
	switch (t)// the switch function is a somewhat different version of if and else, maybe I should change it
	{
	case sort_rand:
	{
		TRandom3 jrand(250);
		const int n = 8;
		for (int i = 0; i < n; i = i + 1) {
			xv.push_back(jrand.Uniform(20.0));
			yv.push_back(jrand.Uniform(6.5, 36.5));
		}
		sort(xv.begin(), xv.end());
		sort(yv.begin(), yv.end());
	}
	return;
	case trig:
	{
		for (int i = 0; i < 8; i++){
			double x = i*0.55 + .3;
			xv.push_back(x);
			yv.push_back(sin(x));
		}
	}
	return;
	}
}

// Prints vector(used for testing)
void print_vector(vector<double> myvec) {
	cout << "[ ";
	for (vector<double>::iterator it = myvec.begin(); it != myvec.end(); it++) {
		cout << *it << ", ";   // NO endl
	}
	cout << "]" << endl;
}

// main of the program
void naturalCubic(){

	// initialize data
	vector<double> xv;
	vector<double> yv;
	fill_test_data(xv, yv, sort_rand);// here I chose random data

	// print vectors to check that they were properly filled
	print_vector(xv);
	print_vector(yv);

	//Load graph from data
	TGraph* gr = new TGraph(xv.size(), &xv[0], &yv[0]);
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(21);

	//make canvas and draw graph
	Double_t w = 600;
	Double_t h = 600;
	TCanvas *c1 = new TCanvas("c hello", "Test data and interpolation", w, h);
	c1->SetWindowSize(w + (w - c1->GetWw()), h + (h - c1->GetWh()));
	c1->SetGrid();
	gr->Draw("AP");

	// Create and interpolator, and superimpose the curve on top of the plotted data points.
	// Here you can chose a fit function to use from the Gnu GSL, documented in section 27.3 of the manual
	//Posible fits to use: kLINEAR, kPOLYNOMIAL, kCSPLINE, kCSPLINE_PERIODIC, kAKIMA, kAKIMA_PERIODIC
	ROOT::Math::Interpolation::Type itype;
	itype = ROOT::Math::Interpolation::kCSPLINE;
	ROOT::Math::Interpolator *inter = new ROOT::Math::Interpolator(xv, yv, itype);
	gInter = inter;

	//draw the fit function 
	TF1 *gr2 = new TF1("interp", interp_wrap, xv[0], *(xv.end() - 1), 0);
	gr2->Draw("SAME");

}

