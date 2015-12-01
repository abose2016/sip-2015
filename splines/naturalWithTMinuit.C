#include <iomanip>
#include "TMinuit.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include <string>
#include "TSystem.h"
#include "TF1.h"
#include "TLatex.h"
#include "TStyle.h"
#include <sstream>
#include <gsl/gsl_spline.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TRandom3.h"

using namespace std;

//Global variables
int nControl = 15;
vector<double> xControlCubic, yControlCubic;
vector<double> xVector, yVector;

gsl_interp_accel *acc = gsl_interp_accel_alloc ();
gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, nControl);	

//______________________________________________________________________________
double ComputeChi2(vector< double > vpar)
{
	//calculate chisquare
	double chisq = 0;
	for (int i = 0; i < (int)yVector.size(); i++) {
//		double delta = (yVector.at(i) - func(xVector.at(i), vpar)) / yErrorVector.at(i);
		chisq += delta*delta;
	}
	return chisq;
}

//______________________________________________________________________________
void fcn(int &npar, double *gin, double &f, double *par, int iflag, int seed = 250, double lowerBound= 5, double upperBound= 20, double lowerErrorBound = .5, double upperErrorBound= 5)
{
	//Populate data
	TRandom3 *jrand = new TRandom3(seed);
	for (int i = 0; i < nControl; ++i)
	{
		xi = ((double) nControl / (nControl - 1)) * i;
	 	yi = jrand->Uniform(lowerBound,upperBound);
		yErrorVector.push_back(jrand->Uniform(lowerErrorBound, upperErrorBound));
		xControlCubic.push_back(xi);
		yControlCubic.push_back(yi);
	}

   gsl_spline_init (spline, &xControlCubic[0], &yControlCubic[0], nControl);

    int nSpline = 1+int((xControlCubic.back()-xControlCubic.front())/stepSpline);
    vector <double> xSpline, ySpline;
    for (int i= 0; i < nSpline; i++)
    {
        xSpline.push_back(xControlCubic.front()+i*stepSpline);
        ySpline.push_back(gsl_spline_eval (spline, xSpline.at(i), acc));
    }

	//free memory
	gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
	delete jrand;
    
	vector< double > vpar;
	for(int i = 0; i < nControl; i++) {
		vpar.push_back(yControlCubic[i]);	
	}
	f = ComputeChi2(vpar);
}


//______________________________________________________________________________
// Construct a graph from vectors and error vectors
TGraphErrors *LoadGraphFromVectorsWithError(vector<double> xVector, vector<double> yVector, vector<double> yErrorVector, string xTitle, string yTitle)
{
	int n = xVector.size();
	if ((xVector.size() == yVector.size())&&(yVector.size() == yErrorVector.size())
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
	FillRandVectors(xVector, yVector, yErrorVector, npts, seed);
	TGraphErrors *g1 = LoadGraphFromVectorsWithError(xVector, yVector, yErrorVector, "X Axis (arbitrary units)", "Y Axis (arbitrary units)");

	// Set starting values and step sizes for parameters
	int npar = npts;
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

	TGraph *cGraph = naturalCubic(stepSpline, xVector, yVector);
	cGraph->SetLineColor(kRed);


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
	cGraph->Draw("SAME");
	tex.DrawLatex(0.18,0.21,slatex1.str().c_str());
	tex.DrawLatex(0.18,0.16,slatex2.str().c_str());


}

