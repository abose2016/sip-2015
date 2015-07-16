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
#include <time.h>

using namespace std;

//Data definition
gsl_interp_accel *acc;
gsl_spline *spline;
gsl_bspline_workspace *bw;

vector<double> xData, yData, yErrorData;

//______________________________________________________________________________
double ComputeChi2(vector< double > ySplineVector)
{
	//calculate chisquare
	double chisq = 0;
	for (int i = 0; i < (int)yData.size(); i++) 
	{
		double delta = (yData.at(i) - ySplineVector.at(i) )/ yErrorData.at(i);
		chisq += delta*delta;
	}
	return chisq; 
}

//______________________________________________________________________________
void fcn(int &, double *, double &f, double *par, int )
{
	vector< double > vpar;
	for(int i = 0; i < (int)xData.size(); i++)	vpar.push_back(par[i]); //used to store the values of the par array

	gsl_spline_init (spline, &xData[0], &vpar[0], xData.size()); //initialize the spline

	vector <double> ySpline;
	for (int i= 0; i < (int)xData.size(); i++) ySpline.push_back(gsl_spline_eval (spline, xData[i], acc));

	f = ComputeChi2(ySpline); //calculate chi square - to be minimized by the TMinuit function
}

//______________________________________________________________________________
// Fill y and error vectors with random points from TRandom3
void FillRandVectors(int nPoints, vector<double> &xVector, vector< double > &yVector, vector< double > &yErrorVector, int seed = 250, double lowerBound= 10, double upperBound= 20, double lowerErrorBound = 1, double upperErrorBound= 2)
{
	//Call TRandom3
	TRandom3 *jrand = new TRandom3(seed);
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
 TGraph *cSpline(int nPoints, int npar, vector <double> xDataC, vector <double> yDataC, vector <double> yErrorDataC, double stepSpline =.01)
{
	//Initialize Minuit
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
	vector< double > bestFitParams, e_bestFitParams;
	for(int i=0; i<npar; i++)
	{
		double par, epar;
		myMinuit->GetParameter(i, par, epar); //retrieve best fit parameters
		bestFitParams.push_back(par);
		e_bestFitParams.push_back(epar);
	}

	//Store the best-fit spline in a TGraph
	gsl_spline_init (spline, &xDataC[0], &bestFitParams[0], xDataC.size()); //initialize the spline

	int nPointsSpline = int ((xDataC[nPoints-1]-xDataC[0])/stepSpline); //calculate the number of points of the spline
	vector< double > xSpline, ySpline;
	for (int i= 0; i < nPointsSpline; i++){
		xSpline.push_back(xDataC[0]+i*stepSpline);
		ySpline.push_back(gsl_spline_eval (spline, xSpline.back(), acc));
	}

	
	// Graph the spline
	TGraph *grSplineC = new TGraph(nPointsSpline, &xSpline[0], &ySpline[0]);
	grSplineC->SetTitle("");
	grSplineC->GetXaxis()->SetTitle("X-axis  [A.U.]");
	grSplineC->GetYaxis()->SetTitle("Y-axis  [A.U.]");
	
	return grSplineC;

}
//_________________________________________________________________________
TGraph *bSpline( int nControl, int npar, vector <double> xDataB, vector <double> yDataB, vector <double> yErrorDataB, double stepSpline = 0.01)
{
	//Initialize variables
	double xmin = 0;
	double xmax = 15;

	//Declare and allocate memory to compose data set of control points
	gsl_vector *xControl = gsl_vector_alloc(nControl);
	gsl_vector *yControl = gsl_vector_alloc(nControl);
	gsl_vector *w = gsl_vector_alloc(nControl);

//Populate data set with monotonically increasing, uniform x values and randomized y values
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
	gsl_bspline_knots_uniform(xmin, xmax, bw);

	//Set up the variables for the fit matrix
	gsl_vector *B = gsl_vector_alloc(npar);
	gsl_matrix *X = gsl_matrix_alloc(nControl, npar);
	for (int i = 0; i < nControl; ++i)
	{
		gsl_bspline_eval(gsl_vector_get(xControl, i), B, bw);//Compute B_j(xi) for all j 
		for (int j = 0; j < npar; ++j)	gsl_matrix_set(X, i, j, gsl_vector_get(B, j));//Fill in row i of X
	}

	//Declare variables for the fit and allocate their memory
	double chisq;	
	gsl_vector *c = gsl_vector_alloc(npar);
	gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(nControl, npar);

	//Do the fit
	gsl_matrix *cov = gsl_matrix_alloc(npar, npar);
	gsl_multifit_linear(X, yControl, c, cov, &chisq, mw);	

	//Output the curve and store the values of the spline in two vectors
	int nValues = 1+int((xDataB.back() - xDataB.front())/stepSpline);
	vector<double> xValues, yValues, splineError;
	for (int i = 0; i < nValues; i++)
	{
		double xi = xDataB.front()+i*stepSpline;
		gsl_bspline_eval(xi, B, bw);
		double yerr, yi;
		gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
		xValues.push_back(xi);
		yValues.push_back(yi);
		splineError.push_back(yerr);
	}

	//Free the memory used
	gsl_vector_free(xControl);
	gsl_vector_free(yControl);
	gsl_vector_free(B);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_multifit_linear_free(mw);
	gsl_matrix_free(cov);

	TGraph *grSplineB = new TGraph(xValues.size(), &xValues[0], &yValues[0]);
	return grSplineB;
} 

//_________________________________________________________________________________

void integratedSplines() 
{
	//Load the data
	int nPoints = 9;
	vector <double> xData, yData, yErrorData; 
	double seed = 99846895709;
	FillRandVectors(nPoints, xData, yData, yErrorData, seed); //create random vectors for y values and y error vector

	//Intialization of the variables
	const int npar =  nPoints;
	const int orderSpline = 4;
	const int nbreak = npar+2-orderSpline;
	double stepSpline = 0.01;
	acc = gsl_interp_accel_alloc ();
	spline = gsl_spline_alloc (gsl_interp_cspline, nPoints);	
	bw = gsl_bspline_alloc(4, nbreak);
	
	//B-spline
	clock_t tbstart = clock();
	TGraph *bGraph = bSpline(nPoints, npar, xData, yData, yErrorData, stepSpline);
	clock_t tbstop = clock();
	bGraph->SetTitle("");
  bGraph->GetXaxis()->SetTitle("X-axis  [A.U.]");
  bGraph->GetYaxis()->SetTitle("Y-axis  [A.U.]");

	//C-spline
	clock_t tcstart = clock();
	TGraph *cGraph = cSpline(nPoints, npar, xData, yData, yErrorData, stepSpline);
	clock_t tcstop = clock();
	cGraph->SetLineColor(kRed);

	//Control points
	TGraph *pGraph = new TGraph (nPoints, &xData[0], &yData[0]);
	pGraph-> SetMarkerStyle(20);
	pGraph->SetMarkerColor(kBlue);

	//Free the memory for the spline
	gsl_spline_free (spline); //frees the memory used by the spline
 	gsl_interp_accel_free (acc);
	gsl_bspline_free(bw);

	//Draw to canvas
	TCanvas *c1 = new TCanvas("c1", "Graph", 200, 10, 700, 500);
 	c1->cd();
	bGraph->Draw("al");
	cGraph->Draw("same l");
	pGraph->Draw("same p");

	//Prints the computation time
	std::cout<<"B-splines: "<<((float)tbstop-(float)tbstart)/ (CLOCKS_PER_SEC)<<"s"<<std::endl;
	std::cout<<"Natural cubic splines: "<<((float)tcstop-(float)tcstart)/ (CLOCKS_PER_SEC)<<"s"<<std::endl;
}














