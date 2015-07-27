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
#include <time.h>

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
// Loads data from Biteau and Williams 2015
void FillVectorsWithData(vector <double> &xVectorData, vector <double> &yVectorData, vector <double> &yErrorVectorData) 
{
	xVectorData.push_back(0.38);
	xVectorData.push_back(0.80);
	xVectorData.push_back(1.7);
	xVectorData.push_back(3.6);
	xVectorData.push_back(7.6);
	xVectorData.push_back(16.);
	xVectorData.push_back(34.);
	xVectorData.push_back(72.);

	yVectorData.push_back(7.4);
	yVectorData.push_back(12.4);
	yVectorData.push_back(10.8);
	yVectorData.push_back(6.7);
	yVectorData.push_back(0.95);
	yVectorData.push_back(3.79);
	yVectorData.push_back(3.14);
	yVectorData.push_back(9.8);

	yErrorVectorData.push_back(3.7);
	yErrorVectorData.push_back(1.9);
	yErrorVectorData.push_back(0.8);
	yErrorVectorData.push_back(0.7);
	yErrorVectorData.push_back(0.53);
	yErrorVectorData.push_back(0.26);
	yErrorVectorData.push_back(0.25);
	yErrorVectorData.push_back(1.1);
}

//______________________________________________________________________________
// Construct a graph from vectors and error vectors
TGraphErrors *LoadGraphFromVectorsWithError(std::vector< double > xVector, std::vector< double > yVector, std::vector< double > yErrorVector, string xTitle, string yTitle)
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
	} else
	{
		TGraphErrors *gr0 = new TGraphErrors();
		return gr0;
		delete gr0;
	}
}

//______________________________________________________________________________
void tMinuitFit()
{
	clock_t tcstart = clock();

	//Load the data
	int nPoints = 8;
	FillVectorsWithData(xVector, yVector, yErrorVector);

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
	myMinuit->SetPrintLevel(-1);//No output: -1, output:1

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
	std::cout<<std::endl;
	vector< double > bestFitParams, e_bestFitParams;
	for(int i=0; i<npar; i++)
	{
		double par, epar;
		myMinuit->GetParameter(i, par, epar); //retrieve best fit parameters
		bestFitParams.push_back(par);
		e_bestFitParams.push_back(epar);
		std::cout<<par<<" "<<epar<<std::endl;
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

	// Graph the spline
	TGraph *grSpline = new TGraph(nPointsSpline, &xSpline[0], &ySpline[0]);
	clock_t tcstop = clock();
	grSpline->SetTitle("");
	grSpline->GetXaxis()->SetTitle("X-axis  [A.U.]");
	grSpline->GetYaxis()->SetTitle("Y-axis  [A.U.]");

	//Free the memory for the spline
	gsl_spline_free (spline); //frees the memory used by the spline
 	gsl_interp_accel_free (acc);

	//Graph of the data
	TGraphErrors *g1 = LoadGraphFromVectorsWithError(xVector, yVector, yErrorVector, "X Axis (arbitrary units)", "Y Axis (arbitrary units)");

	//Draw
	TCanvas *c1 = new TCanvas("c1", "interpolation", 0, 0, 1000, 800);
	c1->cd();
	grSpline->Draw("al");
	g1->Draw("SAME p");

	//Prints the computation time
	std::cout<<"Natural cubic splines: "<<((float)tcstop-(float)tcstart)/ (CLOCKS_PER_SEC)<<"s"<<std::endl;
}

