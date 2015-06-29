#include <vector>
#include <iostream>
#include "TMinuit.h"

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
void minuitFit()
{
	// x values
	xVector.push_back(1);
	xVector.push_back(2.5);
	xVector.push_back(7.67);
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
	yErrorVector.push_back(.2);
	yErrorVector.push_back(.3);
	yErrorVector.push_back(.8);
	yErrorVector.push_back(.1);
	yErrorVector.push_back(1);

	TCanvas *c1 = new TCanvas("c1", "interpolation", 0, 0, 1000, 800);
	TGraphErrors *g1 = LoadGraphFromVectorsWithError(xVector, yVector, xErrorVector, yErrorVector, "X Axis (arbitrary units)", "Y Axis (arbitrary units)");
	gStyle->SetOptFit(1111);
	g1->Draw("ap");

	// Set starting values and step sizes for parameters
	static vector <double> vstart, step;
	vstart.push_back(3);
	vstart.push_back(1);
	vstart.push_back(.1);
	vstart.push_back(.01);

	step.push_back(.1);
	step.push_back(.1);
	step.push_back(.01);
	step.push_back(.001);

	int npar = 4;
	TMinuit *gMinuit = new TMinuit(npar);  //initialize TMinuit with a maximum of 5 params
	gMinuit->SetFCN(fcn);

	double arglist[10];
	int ierflg = 0;
	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

	//stringstream
	//sstream ss;
	//ss<<"a"<<i;
	//-> char*: 
	//ss.str().c_str()

	gMinuit->mnparm(0, "a1", vstart.at(0), step.at(0), 0, 0, ierflg);
	gMinuit->mnparm(1, "a2", vstart.at(1), step.at(1), 0, 0, ierflg);
	gMinuit->mnparm(2, "a3", vstart.at(2), step.at(2), 0, 0, ierflg);
	gMinuit->mnparm(3, "a4", vstart.at(3), step.at(3), 0, 0, ierflg);

	// Now ready for minimization step
	arglist[0] = 500;
	arglist[1] = 1.;
	gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

	// Print results
	double amin, edm, errdef;
	int nvpar, nparx, icstat;
	gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
	//gMinuit->mnprin(3,amin);
	//graph?

	//double par0, epar0;
	//gMinuit->GetParameter(0, par0, epar0);

}

