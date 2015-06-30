#include <vector>
#include <iostream>
#include "TMinuit.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include <string>
#include "TRandom3.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TF1.h"
#include "TStyle.h"
#include <sstream>

using namespace std;

//Data definition
const int nbins = 5;
vector <double> xVector, yVector, xErrorVector, yErrorVector;

//______________________________________________________________________________
double func(float x, double *par)
{
	double value = (((par[0] * x*x) + (par[1] * x) + par[2]));
	return value;
}

//______________________________________________________________________________
double func2(float x, double par0, double par1, double par2, double epar0, double epar1, double epar2)
{
	double value = (((par0 * x*x) + (par1 * x) + par2));
	return value;
}

//______________________________________________________________________________
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
	//calculate chisquare
	double chisq = 0;
	for (int i = 0; i < nbins; i++) {
		double delta = (yVector.at(i) - func(xVector.at(i), par)) / yErrorVector.at(i);
		chisq += delta*delta;
	}
	f = chisq;
}

//______________________________________________________________________________
// Construct a graph from vectors and error vectors
TGraphErrors *LoadGraphFromVectorsWithError(vector<double> xVector, vector<double> yVector, vector<double> xErrorVector, vector<double> yErrorVector, string xTitle, string yTitle)
{
	int n = xVector.size();

	if ((xVector.size() == yVector.size()) &&
		(yVector.size() == yErrorVector.size()) &&
		(xErrorVector.size() == yErrorVector.size()))
	{
		//Create a graph
		TGraphErrors *gr = new TGraphErrors(n, &xVector[0], &yVector[0], &xErrorVector[0], &yErrorVector[0]);
		gr->SetTitle("");
		gr->SetMarkerStyle(20);
		gr->SetMarkerSize(1.2);
		gr->SetLineWidth(2);
		gr->GetXaxis()->SetTitle(xTitle.c_str());
		gr->GetXaxis()->CenterTitle();
		gr->GetYaxis()->SetTitle(yTitle.c_str());
		gr->GetYaxis()->CenterTitle();
		return gr;
		delete gr;
	}
	else
	{
		TGraphErrors *gr0 = new TGraphErrors();
		return gr0;
		delete gr0;
	}
}

//______________________________________________________________________________
void tMinuitFit()
{
	// x values
	xVector.push_back(1);
	xVector.push_back(2.5);
	xVector.push_back(4.67);
	xVector.push_back(8);
	xVector.push_back(9.1);
	// y values
	yVector.push_back(10);
	yVector.push_back(3.5);
	yVector.push_back(8.67);
	yVector.push_back(25);
	yVector.push_back(2);
	// xError values
	xErrorVector.push_back(0);
	xErrorVector.push_back(0);
	xErrorVector.push_back(0);
	xErrorVector.push_back(0);
	xErrorVector.push_back(0);
	// yError values
	yErrorVector.push_back(1.2);
	yErrorVector.push_back(2.3);
	yErrorVector.push_back(.38);
	yErrorVector.push_back(.21);
	yErrorVector.push_back(1);

	TCanvas *c1 = new TCanvas("c1", "interpolation", 0, 0, 1000, 800);
	TGraphErrors *g1 = LoadGraphFromVectorsWithError(xVector, yVector, xErrorVector, yErrorVector, "X Axis (arbitrary units)", "Y Axis (arbitrary units)");
	gStyle->SetOptFit(1111);
	g1->Draw("ap");

	// Set starting values and step sizes for parameters
	vector<double> vstart, step;
	vstart.push_back(3);
	vstart.push_back(1);
	vstart.push_back(.1);
	vstart.push_back(.01);

	step.push_back(.1);
	step.push_back(.1);
	step.push_back(.01);
	step.push_back(.01);

	int npar = 3;
	TMinuit *gMinuit = new TMinuit(npar);  //initialize TMinuit with a maximum of 5 params
	gMinuit->SetFCN(fcn);

	double arglist[10];
	int ierflg = 0;
	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

	for (int i = 0; i < 3; i++) {
		stringstream ss;
		ss<<"a"<<i;
		gMinuit->mnparm(i, ss.str().c_str(), vstart.at(i), step.at(i), 0, 0, ierflg);
	}

	// Now ready for minimization step
	arglist[0] = 500;
	arglist[1] = 1.;
	gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

	// Print results
	double amin, edm, errdef;
	int nvpar, nparx, icstat;
	gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

	double par0, par1, par2, epar0, epar1, epar2;
	gMinuit->GetParameter(0, par0, epar0);
	gMinuit->GetParameter(1, par1, epar1);
	gMinuit->GetParameter(2, par2, epar2);

	stringstream form;
	form << "(" << par0 << "* x*x) + (" << par1 << "* x) + " << par2;
	TF1 *minFunc = new TF1("minFunc", form.str().c_str(), 0.0, 20.0);
	minFunc->Draw("SAME");
}

