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

vector<double> xData_GLOB, yData_GLOB, yErrorData_GLOB; //temporary global vectors that allow all of the functions to see the chosen data set at any given time

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
int binarySearch(vector<double> vector, double key)
{
	int start=1, end=(int)vector.size();
	int mid=(start+end)/2;

	while(start<=end && vector[mid]!=key)
	{
		if(vector[mid] < key)
		{
			start = mid + 1;
		}
		else 
		{
			end = mid - 1;
		}
		mid = (start+end)/2;
	}

	if(vector[mid] == key)
		return mid; 
	else
		return -1;
}

//______________________________________________________________________________
vector< vector<double> > cSpline(int nPoints, int npar, vector <double> xData, vector <double> yData, vector <double> yErrorData, double stepSpline, double start, double step)
{
	//Populate the global variables//-> DONE IN THE MAIN
	xData_GLOB = xData;//-> DONE IN THE MAIN
	yData_GLOB = yData;//-> DONE IN THE MAIN
	yErrorData_GLOB = yErrorData;//-> DONE IN THE MAIN

	//Initialize Minuit
	vector<double> vstart, vstep;
	for(int i=0; i<npar; i++) 	//set starting values and step sizes for parameters
	{
		vstart.push_back(start);
		vstep.push_back(step);
	}

	TMinuit *myMinuit = new TMinuit(npar);  //initialize TMinuit with a maximum of npar (5) -> PASSED AS ARGUMENT
	myMinuit->SetFCN(fcn);//-> DONE IN THE MAIN
	myMinuit->SetPrintLevel(-1);//No output: -1, output:1//-> DONE IN THE MAIN

	double arglist[10];//-> DONE IN THE MAIN
	int ierflg = 0;//-> DONE IN THE MAIN
	arglist[0] = 1;//-> DONE IN THE MAIN
	myMinuit->mnexcm("SET ERR", arglist, 1, ierflg);//-> DONE IN THE MAIN

	for (int i = 0; i < npar; i++) {//-> DONE IN THE MAIN
		stringstream ss;//-> DONE IN THE MAIN
		ss<<"a"<<i;//-> DONE IN THE MAIN
		myMinuit->mnparm(i, ss.str().c_str(), vstart.at(i), vstep.at(i), 0, 0, ierflg);//-> DONE IN THE MAIN
	}//-> DONE IN THE MAIN

	//Perform the Minuit fit
	arglist[0] = 500;//-> DONE IN THE MAIN
	arglist[1] = 1.;//-> DONE IN THE MAIN
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

	vector< vector<double> > cSplineValues;
	cSplineValues.push_back( vector< double > ());//x
	cSplineValues.push_back( vector< double > ());//y
	for (int i= 0; i < nPointsSpline; i++){
		cSplineValues[0].push_back(xData[0]+i*stepSpline);
		cSplineValues[1].push_back(gsl_spline_eval (spline_GLOB, cSplineValues[0].back(), acc_GLOB));
	}

	return cSplineValues;
}
//_________________________________________________________________________
vector< vector<double> > bSpline(int nControl, int npar, vector <double> xDataB, vector <double> yDataB, vector <double> yErrorDataB, double stepSpline, double xmin, double xmax, gsl_bspline_workspace *bw, gsl_vector *xControl, gsl_vector *yControl, gsl_vector *w /*, gsl_vector *B, gsl_matrix *X, gsl_vector *c, gsl_multifit_linear_workspace *mw, gsl_matrix *cov */)
{
	//Populate gsl vectors with the appropriate data
	for (int i = 0; i < nControl; ++i)
	{//-> DONE IN THE MAIN
		double xi = xDataB.at(i);//-> DONE IN THE MAIN
		double yi = yDataB.at(i);//-> DONE IN THE MAIN
		double sigma = yErrorDataB.at(i);//-> DONE IN THE MAIN
		gsl_vector_set(xControl, i, xi);//-> DONE IN THE MAIN
		gsl_vector_set(yControl, i, yi);//-> DONE IN THE MAIN
		gsl_vector_set(w, i, 1.0/(sigma*sigma));//-> DONE IN THE MAIN
	}//-> DONE IN THE MAIN

	//Create a b spline workspace and allocate its memory
	gsl_bspline_knots_uniform(xmin, xmax, bw);//-> DONE IN THE MAIN

	//Set up the variables for the fit matrix
	gsl_vector *B = gsl_vector_alloc(npar);// -> PASSED AS ARGUMENT 
	gsl_matrix *X = gsl_matrix_alloc(nControl, npar);// -> PASSED AS ARGUMENT 
	for (int i = 0; i < nControl; ++i)//-> DONE IN THE MAIN
	{//-> DONE IN THE MAIN
		gsl_bspline_eval(gsl_vector_get(xControl, i), B, bw);//Compute B_j(xi) for all j //-> DONE IN THE MAIN
		for (int j = 0; j < npar; ++j)	gsl_matrix_set(X, i, j, gsl_vector_get(B, j));//Fill in row i of X//-> DONE IN THE MAIN
	}//-> DONE IN THE MAIN

	//Declare variables for the fit and allocate their memory

	gsl_vector *c = gsl_vector_alloc(npar);// -> PASSED AS ARGUMENT 
	gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(nControl, npar);// -> PASSED AS ARGUMENT 

	//Do the fit
	gsl_matrix *cov = gsl_matrix_alloc(npar, npar);// -> PASSED AS ARGUMENT 

	double chisq;
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
	gsl_vector_free(B);//-> DONE IN THE MAIN
	gsl_matrix_free(X);//-> DONE IN THE MAIN
	gsl_vector_free(c);//-> DONE IN THE MAIN
	gsl_multifit_linear_free(mw);//-> DONE IN THE MAIN
	gsl_matrix_free(cov);//-> DONE IN THE MAIN

	//Construct a vector of vectors that will store the xSpline, ySpline, and ySplineError values
	vector< vector<double> > bSplineValues;
	bSplineValues.push_back(xValues);
	bSplineValues.push_back(yValues);
	bSplineValues.push_back(splineError); //this probably does not work because when the number of parameters is the same as the number of points, yerr is NAN

	return bSplineValues;
} 

//_________________________________________________________________________________
void integratedSplinesV4b(double seed = 231) 
{
	//Load the data
	int nEvents = 1000; //number of times the data will be randomized
	int nPoints = 9;

	vector< vector<double> > xEvents, yEvents, yErrorEvents; //each of these is a vector of vectors (of randomized data points)
	for(int i = 0; i < nEvents; i++) 
	{
		vector <double> xData, yData, yErrorData; //temporary vectors that are only used to get random values from FillRand function
		FillRandVectors(nPoints, xData, yData, yErrorData, seed*(i+1)); //populates random vectors for y values and y error vector
		xEvents.push_back(xData);
		yEvents.push_back(yData);
		yErrorEvents.push_back(yErrorData);
	}

	//Intialization of the variables
	const int npar = nPoints;
	const int orderSpline = 4;
	const int nbreak = npar+2-orderSpline;
	double stepSpline = 0.01;
	double xminBSplineWorkspace = 0;
	double xmaxBSplineWorkspace = 9;
	double startCSplineWorkspace = 15.;
	double stepCSplineWorkspace = 1.5;

	acc_GLOB = gsl_interp_accel_alloc ();
	spline_GLOB = gsl_spline_alloc (gsl_interp_cspline, nPoints);	
	gsl_bspline_workspace *bw = gsl_bspline_alloc(orderSpline, nbreak); // -> PASSED AS ARGUMENT
	
	//B- and C-splines
	clock_t tbstart, tbstop, tcstart, tcstop;
	vector <double> timeb, timec;

	vector< vector<double> > xBSplineValues, yBSplineValues;
	vector< vector<double> > xCSplineValues, yCSplineValues; 

	//Setup for the B-spline
	//Declare and allocate memory to compose data set of control points
	gsl_vector *xControl = gsl_vector_alloc(nPoints);
	gsl_vector *yControl = gsl_vector_alloc(nPoints);
	gsl_vector *w = gsl_vector_alloc(nPoints);

	for(int i = 0; i < (int)xEvents.size(); i++)
	{
		tbstart = clock();
		vector< vector<double> > bSplineValues = bSpline(nPoints, npar, xEvents.at(i), yEvents.at(i), yErrorEvents.at(i), stepSpline, xminBSplineWorkspace, xmaxBSplineWorkspace, bw, xControl, yControl, w);
		tbstop = clock();
		timeb.push_back(((float)tbstop-(float)tbstart)/ (CLOCKS_PER_SEC/1000.) );		

		xBSplineValues.push_back(bSplineValues.at(0));
		yBSplineValues.push_back(bSplineValues.at(1));

		tcstart = clock();
		vector< vector<double> > cSplineValues = cSpline(nPoints, npar, xEvents.at(i), yEvents.at(i), yErrorEvents.at(i), stepSpline, startCSplineWorkspace, stepCSplineWorkspace);
		tcstop = clock();
		timec.push_back(((float)tcstop-(float)tcstart)/ (CLOCKS_PER_SEC/1000.) );		

		xCSplineValues.push_back(cSplineValues.at(0));
		yCSplineValues.push_back(cSplineValues.at(1));
	}

	//Histograms______________________________________________________________________________________

	//Time
	int nbins = 100;
	double xlow = 0;
	double xup = 1.;

	TH1D *hTimeB = new TH1D("Time","; time [ms]; number of runs", nbins, xlow, xup); 
	hTimeB->SetStats(0);
	TH1D *hTimeC = new TH1D("TimeC","; time [ms]; number of runs", nbins, xlow, xup); 
	hTimeC->SetLineColor(kRed);
	hTimeC->SetStats(0);

	for(int i=0; i<(int)timec.size(); i++) 
	{
		hTimeB->Fill(timeb.at(i));
		hTimeC->Fill(timec.at(i));	
	}

	//Interpolation distance -> Should FIND THE RIGHT J FOR SPLINE
	vector <double> interpB, interpC;
	for(int i = 0; i < (int)yEvents.size(); i++)
	{
		for(int j = 0; j < (int)yEvents[i].size(); j++)
		{
			double key = xEvents[i][j];
			int indexForB = binarySearch(xBSplineValues[i], key);
			int indexForC = binarySearch(xCSplineValues[i], key);
			std::cout << "B: " << indexForB << " C: " << indexForC << endl;
			if(indexForB != -1) 
				interpB.push_back( (yEvents[i][j]-yBSplineValues[i][indexForB])/yErrorEvents[i][indexForB] );
			if(indexForC != -1)
				interpC.push_back( (yEvents[i][j]-yCSplineValues[i][indexForC])/yErrorEvents[i][indexForC] );
		}
	}	

	int nbinsI = 40;
	int xlowI = -4;
	int xupI = 4;
	TH1D *hInterpB = new TH1D("Interp B","; xAxis; yAxis", nbinsI, xlowI, xupI); 
	for(int i=0; i<(int)interpB.size(); i++) hInterpB->Fill(interpB.at(i));
	hInterpB->SetStats(0);

	TH1D *hInterpC = new TH1D("Interp C","; xAxis; yAxis", nbinsI, xlowI, xupI); 
	for (int i=0; i<(int)interpC.size(); i++) hInterpC->Fill(interpC.at(i));
	hInterpC->SetLineColor(kRed);
	hInterpC->SetStats(0);	

	//Draws______________________________________________________________________________________

	//Interpolation distance
	TCanvas *c1 = new TCanvas("c1", "Interpolation distance");
	c1->cd();
	hInterpB->Draw("");
	hInterpC->Draw("same");

	//Time
	TLegend *leg = new TLegend(0.75,0.70,0.4,0.85);
	leg->SetLineColor(kWhite); 
	leg->SetFillColor(kWhite);
	leg->SetMargin(0.3); 
	leg->AddEntry(hTimeB,"b-spline","l");
	leg->AddEntry(hTimeC,"c-spline","l");

	TCanvas *c2 = new TCanvas("c2", "Computation time");
	c2->cd();
	hTimeB->Draw("");
	hTimeC->Draw("same");
	leg-> Draw();

	//Free the memory for the spline
	gsl_spline_free (spline_GLOB); //frees the memory used by the spline
 	gsl_interp_accel_free (acc_GLOB);
	gsl_bspline_free(bw);

	gsl_vector_free(xControl);//-> DONE IN THE MAIN
	gsl_vector_free(yControl);//-> DONE IN THE MAIN

}

