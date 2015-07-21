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
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <sstream>
#include <string>
#include <time.h>
#include "TH1.h"
#include "TLegend.h"

using namespace std;

//Data definition
gsl_interp_accel *acc_GLOB;
gsl_spline *spline_GLOB;
gsl_bspline_workspace *bw_GLOB;

vector< vector<double> > xEvents_GLOB, yEvents_GLOB, yErrorEvents_GLOB; //each of these is a vector of vectors (of randomized data points)
vector< vector<double> > xBSplineValues_GLOB, yBSplineValues_GLOB, yBSplineErrorValues_GLOB; //each of these is a vector of vectors (of the appropriate spline values at the same index as the original data points)
vector< vector<double> > xCSplineValues_GLOB, yCSplineValues_GLOB; 
vector<double> xData_GLOB, yData_GLOB, yErrorData_GLOB; //temporary global vectors that allow all of the functions to see the chosen data set at any given time

//B-spline specific________________________________________________________________________________
gsl_vector *xControl, *yControl, *w, *B, *c;
gsl_matrix *X, *cov;
gsl_multifit_linear_workspace *mw

//C-spline specific________________________________________________________________________________
vector<double> vstart, vstep;
TMinuit *myMinuit;

//______________________________________________________________________________
double ComputeChi2(vector< double > ySplineVector)
{
	//calculate chisquare
	double chisq = 0;
	for (int i = 0; i < (int)yData_GLOB.size(); i++) 
	{
		double delta = (yData_GLOB.at(i) - ySplineVector.at(i) )/ yErrorData_GLOB.at(i);
		chisq += delta*delta;
	}
	return chisq; 
}

//______________________________________________________________________________
void fcn(int &, double *, double &f, double *par, int )
{
	vector< double > vpar;
	for(int i = 0; i < (int)xData_GLOB.size(); i++)	vpar.push_back(par[i]); //used to store the values of the par array

	gsl_spline_init (spline_GLOB, &xData_GLOB[0], &vpar[0], xData_GLOB.size()); //initialize the spline

	vector <double> ySpline;
	for (int i= 0; i < (int)xData_GLOB.size(); i++) ySpline.push_back(gsl_spline_eval (spline_GLOB, xData_GLOB[i], acc_GLOB));

	f = ComputeChi2(ySpline); //calculate chi square - to be minimized by the TMinuit function
}

//______________________________________________________________________________
// Fill y and error vectors with random points from TRandom3
void FillRandVectors(int nPoints, vector<double> &xVector, vector< double > &yVector, vector< double > &yErrorVector, int seed = 250, double lowerBound= 10, double upperBound= 20, double lowerErrorBound = 1, double upperErrorBound= 2)
{
	//Call TRandom3
	TRandom3 *jrand = new TRandom3(seed);
	xVector.clear();
	yVector.clear();
	yErrorVector.clear();

	for (int i = 0; i < nPoints; i++)
	{
		xVector.push_back(i);
		yVector.push_back(jrand->Uniform(lowerBound, upperBound));
		yErrorVector.push_back(jrand->Uniform(lowerErrorBound, upperErrorBound));
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
vector< vector<double> > cSpline(int nPoints, int npar, vector <double> xData, vector <double> yData, vector <double> yErrorData, double stepSpline =.01, double start = 10., double step = 0.01)
{
	//Populate the global variables
	xData_GLOB = xData;
	yData_GLOB = yData;
	yErrorData_GLOB = yErrorData;

	//Initialize Minuit
	for(int i=0; i<npar; i++) 	//set starting values and step sizes for parameters
	{
		vstart.push_back(start);
		vstep.push_back(step);
	}

	myMinuit = new TMinuit(npar);  //initialize TMinuit with a maximum of npar
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
	vector< double > bestFitParams, e_bestFitParams;
	for(int i=0; i<npar; i++)
	{
		double par, epar;
		myMinuit->GetParameter(i, par, epar); //retrieve best fit parameters
		bestFitParams.push_back(par);
		e_bestFitParams.push_back(epar);
	}

	//Store the best-fit spline in a TGraph
	gsl_spline_init (spline_GLOB, &xData[0], &bestFitParams[0], xData.size()); //initialize the spline

	int nPointsSpline = int ((xData[nPoints-1]-xData[0])/stepSpline); //calculate the number of points of the spline
	vector< double > xSpline, ySpline;
	for (int i= 0; i < nPointsSpline; i++){
		xSpline.push_back(xData[0]+i*stepSpline);
		ySpline.push_back(gsl_spline_eval (spline_GLOB, xSpline.back(), acc_GLOB));
	}

	//Construct a vector of vectors that will store the xSpline values and the ySpline values
	vector< vector<double> > cSplineValues;
	cSplineValues.push_back(xSpline);
	cSplineValues.push_back(ySpline);

	return cSplineValues;
}
//_________________________________________________________________________
vector< vector<double> > bSpline(int nControl, int npar, vector <double> xDataB, vector <double> yDataB, vector <double> yErrorDataB, double stepSpline = 0.01, double xmin = 0, double xmax = 9)
{
	xControl = gsl_vector_alloc(nControl);
	yControl = gsl_vector_alloc(nControl);
	w = gsl_vector_alloc(nControl);

//Populate gsl vectors with the appropriate data
	for (int i = 0; i < nControl; ++i)
	{
		double xi = xDataB.at(i);
		double yi = yDataB.at(i);
		double sigma = yErrorDataB.at(i);
		gsl_vector_set(xControl, i, xi);
		gsl_vector_set(yControl, i, yi);
		gsl_vector_set(w, i, 1.0/(sigma*sigma));
	}

	//Create a b spline workspace and allocate its memory
	gsl_bspline_knots_uniform(xmin, xmax, bw_GLOB);

	//Set up the variables for the fit matrix
	B = gsl_vector_alloc(npar);
	X = gsl_matrix_alloc(nPoints, npar);

	for (int i = 0; i < nControl; ++i)
	{
		gsl_bspline_eval(gsl_vector_get(xControl, i), B, bw_GLOB);//Compute B_j(xi) for all j 
		for (int j = 0; j < npar; ++j)	gsl_matrix_set(X, i, j, gsl_vector_get(B, j));//Fill in row i of X
	}

	//Declare variables for the fit and allocate their memory
	double chisq;	
	c = gsl_vector_alloc(npar);
	mw = gsl_multifit_linear_alloc(nPoints, npar);

	//Do the fit
	cov = gsl_matrix_alloc(npar, npar);
	gsl_multifit_linear(X, yControl, c, cov, &chisq, mw);	

	//Output the curve and store the values of the spline in two vectors
	int nValues = 1+int((xDataB.back() - xDataB.front())/stepSpline);
	vector<double> xValues, yValues, splineError;
	for (int i = 0; i < nValues; i++)
	{
		double xi = xDataB.front()+i*stepSpline;
		gsl_bspline_eval(xi, B, bw_GLOB);
		double yerr, yi;
		gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
		xValues.push_back(xi);
		yValues.push_back(yi);
		splineError.push_back(yerr);
	}

	//Free the memory used
	gsl_vector_free(xControl);
	gsl_vector_free(yControl);

	//Construct a vector of vectors that will store the xSpline, ySpline, and ySplineError values
	vector< vector<double> > bSplineValues;
	bSplineValues.push_back(xValues);
	bSplineValues.push_back(yValues);
	bSplineValues.push_back(splineError); //this probably does not work because when the number of parameters is the same as the number of points, yerr is NAN

	return bSplineValues;
} 

//_________________________________________________________________________________
void integratedSplinesV3(double seed = 231) 
{
	int nEvents = 1000; //number of times the data will be randomized
	int nPoints = 9;
	vector <double> xData, yData, yErrorData; //temporary vectors that are only used to get random values from FillRand function

	for(int i = 0; i < nEvents; i++) 
	{
		FillRandVectors(nPoints, xData, yData, yErrorData, seed*(i+1)); //populates random vectors for y values and y error vector
		xEvents_GLOB.push_back(xData);
		yEvents_GLOB.push_back(yData);
		yErrorEvents_GLOB.push_back(yErrorData);
	}

	//Intialization of the variables
	const int npar = nPoints; //must be one less than nPoints in order for linearity histogram to work
	double xmin = 0;
	double xmax = nPoints;
	const int orderSpline = 4;
	const int nbreak = npar+2-orderSpline;
	double stepSpline = 0.01;

	acc_GLOB = gsl_interp_accel_alloc ();
	spline_GLOB = gsl_spline_alloc (gsl_interp_cspline, nPoints);	
	bw_GLOB = gsl_bspline_alloc(orderSpline, nbreak);

	//B-spline
	clock_t tbstart;
	clock_t tbstop;
	vector <double> timeb, interpb, linearb;

	for(int i = 0; i < (int)xEvents_GLOB.size(); i++)
	{
		tbstart = clock();
		xData_GLOB = xEvents_GLOB.at(i); //assigning the global variables to the current vector value in the events vector
		yData_GLOB = yEvents_GLOB.at(i);
		yErrorData_GLOB = yErrorEvents_GLOB.at(i);
		vector< vector<double> > bSplineValues = bSpline(nPoints, npar, xData_GLOB, yData_GLOB, yErrorData_GLOB, stepSpline, xmin, xmax);
		xBSplineValues_GLOB.push_back(bSplineValues.at(0));
		yBSplineValues_GLOB.push_back(bSplineValues.at(1));
		yBSplineErrorValues_GLOB.push_back(bSplineValues.at(2));
		xData_GLOB.clear();
		yData_GLOB.clear();
		yErrorData_GLOB.clear();	
		tbstop = clock();
		timeb.push_back(((float)tbstop-(float)tbstart)/ CLOCKS_PER_SEC);		
	}

	//C-spline
	clock_t tcstart;
	clock_t tcstop;
	vector <double> timec;

	for(int i = 0; i < (int)xEvents_GLOB.size(); i++) //loop through each event
	{
		tcstart = clock();
		xData_GLOB = xEvents_GLOB.at(i); //assigning the global variables to the current vector value in the events vector
		yData_GLOB = yEvents_GLOB.at(i);
		yErrorData_GLOB = yErrorEvents_GLOB.at(i);
		vector< vector<double> > cSplineValues = cSpline(nPoints, npar, xData_GLOB, yData_GLOB, yErrorData_GLOB, stepSpline);
		xCSplineValues_GLOB.push_back(cSplineValues.at(0));
		yCSplineValues_GLOB.push_back(cSplineValues.at(1));
		xData_GLOB.clear();
		yData_GLOB.clear();
		yErrorData_GLOB.clear();	
		tcstop = clock();
		timec.push_back(((float)tcstop-(float)tcstart)/ CLOCKS_PER_SEC);						
	}

//Histograms______________________________________________________________________________________

	// time
	int nbinsT = 100;
	double xlowT = 0;
	double xupT = 0.001;

	TH1D *hello2 = new TH1D("Time","Time; time (secs/10,000); number of runs", nbinsT, xlowT, xupT); 
	TH1D *hello3 = new TH1D("TimeC","TimeC; time; number of runs", nbinsT, xlowT, xupT); 
	hello3->SetLineColor(kRed);

	for(int i=0; i<(int)timec.size(); i++) 
	{
		hello2->Fill(timeb.at(i));
		hello3->Fill(timec.at(i));	
	}

/*	// linearity - only works when npar is less than npoints
	int nbinsL = 8;
	int xlowL = 0;
	int xupL = 8;
	vector <double> linear;
	for(int i = 0; i < (int)yEvents_GLOB.size(); i++) 
	{
		for(int j = 1; j < (int)yEvents_GLOB[i].size(); j++)
		{
			linear.push_back(yBSplineErrorValues_GLOB[i][j]/yEvents_GLOB[i][j]);
		}
	}	
	
	TH1D *linearHist = new TH1D("Linearity check","Linearity; xAxis; yAxis", nbinsL, xlowL, xupL); 

	for (int i = 0; i < (int)linear.size(); i++) linearHist->Fill(linear.at(i));
	TCanvas *canLinear = new TCanvas("c5", "Linearity Histogram");
	canLinear->cd();
	linearHist->Draw(""); */

	// interpolation for B spline
	int nbinsI = 40;
	int xlowI = -4;
	int xupI = 4;
	vector <double> interpB;
	for(int i = 0; i < (int)yEvents_GLOB.size(); i++) 
	{
		for(int j = 1; j < (int)yEvents_GLOB[i].size(); j++)
		{
			interpB.push_back(yEvents_GLOB[i][j]-yBSplineValues_GLOB[i][j]);
		}
	}	

	TH1D *interpHistB = new TH1D("Interp B","Interp; xAxis; yAxis", nbinsI, xlowI, xupI); 
	for ( int i=0; i<(int)interpB.size(); i++) interpHistB->Fill(interpB.at(i));

	// interpolation for C spline
	vector <double> interpC;
	for(int i = 0; i < (int)yEvents_GLOB.size(); i++) 
	{
		for(int j = 1; j < (int)yEvents_GLOB[i].size(); j++)
		{
			interpC.push_back(yEvents_GLOB[i][j]-yCSplineValues_GLOB[i][j]);
		}
	}	

	TH1D *interpHistC = new TH1D("Interp C","Interp; xAxis; yAxis", nbinsI, xlowI, xupI); 
	for ( int i=0; i<(int)interpC.size(); i++) interpHistC->Fill(interpC.at(i));
	interpHistC->SetLineColor(kRed);
	
	TCanvas *canInterp = new TCanvas("c6", "Interpolation Histogram");
	canInterp->cd();
	interpHistB->Draw("");
	interpHistC->Draw("same");

	//Legends
	TLegend *leg = new TLegend(0.75,0.70,0.4,0.85);
	leg->SetLineColor(kWhite); 
	leg->SetFillColor(kWhite);
	leg->SetMargin(0.3); 
	leg->AddEntry(hello2,"b-spline","l");
	leg->AddEntry(hello3,"c-spline","l");

	//Draw to canvas
	TCanvas *can1 = new TCanvas("c2", "Timing");
	can1->cd();
	hello2->Draw("");
	hello3->Draw("same");
	leg-> Draw();

	//Free the memory for the spline
	gsl_spline_free (spline_GLOB); //frees the memory used by the spline
 	gsl_interp_accel_free (acc_GLOB);
	gsl_bspline_free(bw_GLOB);

	gsl_vector_free(B);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_multifit_linear_free(mw);
	gsl_matrix_free(cov);


}

