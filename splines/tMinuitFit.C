#include <vector>
#include <iostream>
#include <iomanip>
#include "TMinuit.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include <string>
#include "TRandom3.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TF1.h"
#include "TLatex.h"
#include "TStyle.h"
#include <sstream>

using namespace std;

//Data definition
vector <double> xVector, yVector, xErrorVector, yErrorVector;
	

//______________________________________________________________________________
double func(float x, vector< double > vpar)
{
	double value = vpar[0] + vpar[1]*x + vpar[2]*x*x;
	return value;
}

//______________________________________________________________________________
double ComputeChi2(vector< double > vpar)
{
	//calculate chisquare
	double chisq = 0;
	for (int i = 0; i < (int)yVector.size(); i++) {
		double delta = (yVector.at(i) - func(xVector.at(i), vpar)) / yErrorVector.at(i);
		chisq += delta*delta;
	}
	return chisq;
}

//______________________________________________________________________________
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
	vector< double > vpar;
	vpar.push_back(par[0]);
	vpar.push_back(par[1]);
	vpar.push_back(par[2]);

	f = ComputeChi2(vpar);
}



//______________________________________________________________________________
// Fill x,y and error vectors with random points from TRandom#
void FillRandVectors(vector<double> &xVector, vector< double > &yVector, vector< double > &xErrorVector, vector< double > &yErrorVector, int n, int seed = 250, double lowerBound= 5, double upperBound= 20, double lowerErrorBound = .5, double upperErrorBound= 5)
{
	//Call TRandom3
	TRandom3 *jrand = new TRandom3(seed);
	for (int i = 0; i<n; i = i + 1)
	{
		xVector.push_back(jrand->Uniform(upperBound));
		yVector.push_back(jrand->Uniform(lowerBound, upperBound));
		xErrorVector.push_back(0);
		yErrorVector.push_back(jrand->Uniform(lowerErrorBound, upperErrorBound));
	}
	delete jrand;
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
		TGraphErrors *gr = new TGraphErrors(n, &xVector[0], &yVector[0], 0, &yErrorVector[0]);
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

  
	double seed = 59050830;
	int npts = 5;
	FillRandVectors(xVector, yVector, xErrorVector, yErrorVector, npts, seed);
	TGraphErrors *g1 = LoadGraphFromVectorsWithError(xVector, yVector, xErrorVector, yErrorVector, "X Axis (arbitrary units)", "Y Axis (arbitrary units)");

	// Set starting values and step sizes for parameters
	int npar = 3;
	vector<double> vstart, vstep;
	for(int i=0; i<npar; i++)
	{
		vstart.push_back(0.);
		vstep.push_back(.1);
	}

	TMinuit *myMinuit = new TMinuit(npar);  //initialize TMinuit with a maximum of 5 params
	myMinuit->SetFCN(fcn);
	double arglist[10];
	int ierflg = 0;
	arglist[0] = 1;
	myMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

	for (int i = 0; i < npar; i++) {
		stringstream ss;
		ss<<"a"<<i;
		myMinuit->mnparm(i, ss.str().c_str(), vstart.at(i), vstep.at(i), 0, 0, ierflg);
	}

	// Now ready for minimization step
	arglist[0] = 500;
	arglist[1] = 1.;
	myMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

	// Print results
	double amin, edm, errdef;
	int nvpar, nparx, icstat;
	myMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

	vector< double > vpar;
	for(int i=0; i<npar; i++)
	{
		double par, epar;
		myMinuit->GetParameter(i, par, epar);
		vpar.push_back(par);
	}

	stringstream form;
	form << vpar[0] << "+" << vpar[1] << "*x + " << vpar[2] <<"*x*x";
	TF1 *minFunc = new TF1("minFunc", form.str().c_str(), 0.0, 20.0);

	//Retrieve the chi2 and number of degrees of freedom
	double chi2 = ComputeChi2(vpar);
	int ndf = npts-npar;
	double prob = TMath::Prob(chi2,ndf);

	//Drawing
	TLatex tex;
	tex.SetTextAlign(12);
	tex.SetTextSize(0.04);
	tex.SetNDC();

	stringstream slatex1;
	slatex1<<setprecision(3)<<"#chi^{2} / ndf = "<<chi2<<" / "<<ndf; 

	stringstream slatex2;
	slatex2<<setprecision(2)<<"P(#chi^{2},ndf) = "<<prob*100<<" %"; 

	TCanvas *c1 = new TCanvas("c1", "interpolation", 0, 0, 1000, 800);
	c1->cd();
	g1->Draw("ap");
	minFunc->Draw("SAME");
	tex.DrawLatex(0.18,0.21,slatex1.str().c_str());
	tex.DrawLatex(0.18,0.16,slatex2.str().c_str());


}

