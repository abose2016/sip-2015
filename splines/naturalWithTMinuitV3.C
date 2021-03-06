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
#include <gsl/gsl_spline.h>

using namespace std;

//Data definition
vector <double> xVector, yVector, yErrorVector; //represents the original set of data to fit with a spline
gsl_interp_accel *acc;
gsl_spline *spline;

//______________________________________________________________________________
double ComputeChi2(vector< double > ySplineVector)
{
	//calculate chisquare
	double chisq = 0;
	for (int i = 0; i < (int)yVector.size(); i++) 
	{
		double delta = (yVector.at(i) - ySplineVector.at(i) )/ yErrorVector.at(i);
		chisq += delta*delta;
	}
	return chisq;
}

//______________________________________________________________________________
void fcn(int &, double *, double &f, double *par, int )
{
	vector< double > vpar;
	for(int i = 0; i < (int)xVector.size(); i++)	vpar.push_back(par[i]); //used to store the values of the par array

	gsl_spline_init (spline, &xVector[0], &vpar[0], xVector.size()); //initialize the spline

	vector <double> ySpline;
	for (int i= 0; i < (int)xVector.size(); i++) ySpline.push_back(gsl_spline_eval (spline, xVector[i], acc));

	f = ComputeChi2(ySpline); //calculate chi square - to be minimized by the TMinuit function
}



//______________________________________________________________________________
// Fill y and error vectors with random points from TRandom#
void FillRandVectors(int nPoints, vector<double> &xVector, vector< double > &yVector, vector< double > &yErrorVector, int seed = 250, double lowerBound= 10, double upperBound= 20, double lowerErrorBound = 1, double upperErrorBound= 2)
{
	//Call TRandom3
	TRandom3 *jrand = new TRandom3(seed);
	for (int i = 0; i < nPoints; i++)
	{
		xVector.push_back(i);
		yVector.push_back(jrand->Uniform(lowerBound, upperBound));
		yErrorVector.push_back(jrand->Uniform(lowerErrorBound, upperErrorBound));

		cout<<yVector[i]<<" "<<yErrorVector[i]<<endl;
	}
	delete jrand;
}
//______________________________________________________________________________
// Construct a graph from vectors and y error vector
TGraphErrors *LoadGraphFromVectorsWithError(vector<double> xVector, vector<double> yVector, vector<double> yErrorVector, string xTitle, string yTitle)
{
	int n = xVector.size();

	if ((xVector.size() == yVector.size()) &&
		(yVector.size() == yErrorVector.size()))
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
	//Load the data
	int nPoints = 8;
	double seed = 250;
	FillRandVectors(nPoints, xVector, yVector, yErrorVector, seed); //create random vectors for y values and y error vector

	//Intialization of the global variables
	acc = gsl_interp_accel_alloc ();
	spline = gsl_spline_alloc (gsl_interp_cspline, nPoints);	

	//Initialize Minuit
	int npar = nPoints;
	vector<double> vstart, vstep;
	for(int i=0; i<npar; i++) 	//set starting values and step sizes for parameters
	{
		vstart.push_back(10.);
		vstep.push_back(0.01);
	}

	TMinuit *myMinuit = new TMinuit(npar);  //initialize TMinuit with a maximum of npar (5)
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

	//Perform the Minuit fit
	arglist[0] = 500;
	arglist[1] = 1.;
	myMinuit->mnexcm("MIGRAD", arglist, 2, ierflg); //minimization

	//Retrieve best-fit parameters
	vector< double > bestFitParams, e_bestFitParams;
	for(int i=0; i<npar; i++)
	{
		double par, epar;
		myMinuit->GetParameter(i, par, epar); //retrieve best fit parameters
		bestFitParams.push_back(par);
		e_bestFitParams.push_back(epar);
	}

	//Store the best-fit spline in a TGraph
	gsl_spline_init (spline, &xVector[0], &bestFitParams[0], xVector.size()); //initialize the spline

	double stepSpline =.01;
	int nPointsSpline = int ((xVector[nPoints-1]-xVector[0])/stepSpline); //calculate the number of points of the spline
	vector< double > xSpline, ySpline;
	for (int i= 0; i < nPointsSpline; i++){
		xSpline.push_back(xVector[0]+i*stepSpline);
		ySpline.push_back(gsl_spline_eval (spline, xSpline.back(), acc));
	}

	TGraph *grSpline = new TGraph(nPointsSpline, &xSpline[0], &ySpline[0]);
	grSpline->SetTitle("");
	grSpline->GetXaxis()->SetTitle("Control Points");
	grSpline->GetYaxis()->SetTitle("Arbitrary Values");

	//Free the memory for the spline
	gsl_spline_free (spline); //frees the memory used by the spline
 	gsl_interp_accel_free (acc);

	//Graph of the data
	TGraph *g1 = new TGraph(xVector.size(), &xVector[0], &yVector[0]);
	g1->SetMarkerStyle(20);
	g1->SetTitle("");
	g1->GetXaxis()->SetTitle("Control Points");
	g1->GetXaxis()->SetLimits(0, 7);
	g1->GetYaxis()->SetTitle("Arbitrary Values");

	//Draw
	TCanvas *c1 = new TCanvas("c1", "interpolation", 0, 0, 1000, 800);
	c1->cd();
	g1->Draw("ap");
	grSpline->Draw("SAME l");
}

