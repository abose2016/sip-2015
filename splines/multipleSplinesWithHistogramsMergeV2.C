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
// Returns the index of value 'key' in param vector
int binarySearch(vector<double> vector, double key)
{
	int start=0, end=(int)vector.size()-1;
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
//Makes a linear interpolator and returns a vector of vectors that plot the interpolation
vector <vector<double> > linearInterpolation( vector<double> xData, vector<double> yData, double stepSpline)
{
	vector <vector<double> > linearInterp;
	linearInterp.push_back( vector<double> () );//x
	linearInterp.push_back( vector<double> () );//y
	for(int i = 0; i < (int)xData.size(); i++)
	{
		if( i == (int)xData.size()-1)
		{
			linearInterp[0].push_back(i);
			linearInterp[1].push_back(yData.at(i));		
		}
		else
		{
			double x0 = xData.at(i);
			double x1 = xData.at(i+1);
			double y0 = yData.at(i);
			double y1 = yData.at(i+1);
			vector< double > xSpline, ySpline;

			for (double j = xData.at(i); j < xData.at(i+1); j += stepSpline)
			{
					//temp variables for interp
					double m = (y1 - y0)/(x1 - x0);
					double b = -m*x0 + y0;
					double y = m * (j) + b; 
					linearInterp[0].push_back(j);
					linearInterp[1].push_back(y); 
			}
		} 
	} 

	return linearInterp;
}

//___________________________________________________________________________
//Returns a vector of the average distances between two vectors of vectors
vector<double> diff(vector <vector<double> > yInterp1, vector <vector<double> > yInterp2) 
{
	double sum = 0;
	vector<double> differences;

	for(int i = 0; i < (int)yInterp1.size(); i++) 
	{
		for(int j = 0; j < (int) yInterp1[0].size(); j++)
		{
			sum += (yInterp2[i][j] - yInterp1[i][j])*(yInterp2[i][j] - yInterp1[i][j]);
		}
		differences.push_back(sqrt(sum/(int)yInterp1[0].size()));
		sum = 0;
	}

	return differences;
}

//______________________________________________________________________________
// Uses TMinuit to conduct minimization and calculation of natural cubic spline
vector< vector<double> > cSpline(int nPoints, int npar, vector <double> xData, double stepSpline, TMinuit *myMinuit)
{
	double arglist[10];
	int ierflg = 0;

	//Perform the Minuit fit
	arglist[0] = 500;
	arglist[1] = 1.;
	myMinuit->mnexcm("MIGRAD", arglist, 2, ierflg); //minimization

	//Retrieve best-fit parameters
	vector< double > bestFitParams, e_bestFitParams;
	for(int i=0; i<npar; i++)
	{
		double par, epar;
		myMinuit->GetParameter(i, par, epar);
		bestFitParams.push_back(par);
		e_bestFitParams.push_back(epar);
	}

	gsl_spline_init (spline_GLOB, &xData[0], &bestFitParams[0], xData.size()); //initialize the spline

	int nPointsSpline = int (1.+(xData[nPoints-1]-xData[0])/stepSpline); //calculate the number of points of the spline
	//1+ -> IMPORTANT OTHERWISE YOU MISS THE LAST POINT

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
// Uses GSL libraries to conduct minimization and fit for B-spline
vector< vector<double> > bSpline(vector <double> xDataB, double stepSpline, gsl_bspline_workspace *bw, gsl_vector *yControl, gsl_vector *B, gsl_matrix *X, gsl_vector *c, gsl_multifit_linear_workspace *mw, gsl_matrix *cov)
{	
	//Do the fit
	double chisq;
	gsl_multifit_linear(X, yControl, c, cov, &chisq, mw);

	//Output the curve and store the values of the spline in two vectors
	int nValues = 1+int((xDataB.back() - xDataB.front())/stepSpline);
	vector< vector<double> > bSplineValues;
	bSplineValues.push_back( vector<double> () );//x
	bSplineValues.push_back( vector<double> () );//y
	for (int i = 0; i < nValues; i++)
	{
		double xi = xDataB.front()+i*stepSpline;
		gsl_bspline_eval(xi, B, bw);
		double yerr, yi;
		gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
		bSplineValues[0].push_back(xi);
		bSplineValues[1].push_back(yi);
	}

	return bSplineValues;
} 

//_________________________________________________________________________________
void multipleSplinesWithHistogramsMerge(int iEventLook = 163, int nEvents = 10000, int nPoints = 9, double seed = 231) 
{
	double lowerBound = 10; //bounds for random vector function
	double upperBound = 20;
	double lowerErrorBound = 1;
	double upperErrorBound = 2;

	//Load the data
	vector< vector<double> > xEvents, yEvents, yErrorEvents; //each of these is a vector of vectors (of randomized data points)
	for(int i = 0; i < nEvents; i++) 
	{
		vector <double> xData, yData, yErrorData; //temporary vectors that are only used to get random values from FillRand function
		FillRandVectors(nPoints, xData, yData, yErrorData, seed*(i+1), lowerBound, upperBound, lowerErrorBound, upperErrorBound); //populates random vectors for y values and y error vector
		xEvents.push_back(xData);
		yEvents.push_back(yData);
		yErrorEvents.push_back(yErrorData);
	}

	//Intialization of variables
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
	
	//B- and C-splines
	clock_t tbstart, tbstop, tcstart, tcstop;
	vector <double> timeb, timec;

	vector< vector<double> > xBSplineValues, yBSplineValues;
	vector< vector<double> > xCSplineValues, yCSplineValues; 
	vector< vector<double> > xLinearInterpValues, yLinearInterpValues; 
//Setup for the C-spline_________________________________________________________________________
	TMinuit *myMinuit = new TMinuit(npar);  //initialize TMinuit with a maximum of npar 
	myMinuit->SetFCN(fcn);
	myMinuit->SetPrintLevel(-1);//No output: -1, output:1

	double arglist[10];
	int ierflg = 0;
	arglist[0] = 1;
	myMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

	//Initialize Minuit
	vector<double> vstart, vstep;
	for(int i=0; i<npar; i++) 	//set starting values and step sizes for parameters
	{
		vstart.push_back(startCSplineWorkspace);
		vstep.push_back(stepCSplineWorkspace);
	}

	for (int i = 0; i < npar; i++) {
		stringstream ss;
		ss<<"a"<<i;
		myMinuit->mnparm(i, ss.str().c_str(), vstart.at(i), vstep.at(i), 0, 0, ierflg);
	}

//Setup for the B-spline_________________________________________________________________________
	gsl_bspline_workspace *bw = gsl_bspline_alloc(orderSpline, nbreak);
	gsl_vector *xControl = gsl_vector_alloc(nPoints);
	gsl_vector *yControl = gsl_vector_alloc(nPoints);
	gsl_vector *w = gsl_vector_alloc(nPoints);

	//Create a b spline workspace and allocate its memory
	gsl_bspline_knots_uniform(xminBSplineWorkspace, xmaxBSplineWorkspace, bw);

	//Set up the variables for the fit matrix
	gsl_vector *B = gsl_vector_alloc(npar); 
	gsl_vector *c = gsl_vector_alloc(npar);
	gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(nPoints, npar);  
	gsl_matrix *cov = gsl_matrix_alloc(npar, npar);
	gsl_matrix *X = gsl_matrix_alloc(nPoints, npar);

//Looping begins for the calculations of the B and C-splines for each event_____________________
	for(int i = 0; i < (int)xEvents.size(); i++)
	{
		//Populate the global variables
		xData_GLOB = xEvents.at(i);
		yData_GLOB = yEvents.at(i);
		yErrorData_GLOB = yErrorEvents.at(i);

		//Populate gsl vectors with the appropriate data
		for (int j = 0; j < nPoints; ++j)
		{
			gsl_vector_set(xControl, j, xEvents[i][j]);
			gsl_vector_set(yControl, j, yEvents[i][j]);
			gsl_vector_set(w, j, 1.0/(yErrorEvents[i][j]*yErrorEvents[i][j]));
		}
	
		for (int n = 0; n < nPoints; ++n)
		{
			gsl_bspline_eval(gsl_vector_get(xControl, n), B, bw);//Compute B_j(xi) for all j 
			for (int l = 0; l < npar; ++l)	gsl_matrix_set(X, n, l, gsl_vector_get(B, l));//Fill in row i of X
		}

		tbstart = clock();
		vector< vector<double> > bSplineValues = bSpline(xEvents.at(i), stepSpline, bw, yControl, B, X, c, mw, cov);
		tbstop = clock();
		timeb.push_back(((float)tbstop-(float)tbstart)/ (CLOCKS_PER_SEC/1000.) );		

		xBSplineValues.push_back(bSplineValues.at(0));
		yBSplineValues.push_back(bSplineValues.at(1));

		tcstart = clock();
		vector< vector<double> > cSplineValues = cSpline(nPoints, npar, xEvents.at(i), stepSpline, myMinuit);
		tcstop = clock();
		timec.push_back(((float)tcstop-(float)tcstart)/ (CLOCKS_PER_SEC/1000.) );		

		xCSplineValues.push_back(cSplineValues.at(0));
		yCSplineValues.push_back(cSplineValues.at(1));
		
	//Linear interpolation
		vector <vector<double> > interpPoints = linearInterpolation(xData_GLOB, yData_GLOB, stepSpline);
		xLinearInterpValues.push_back(interpPoints.at(0));
		yLinearInterpValues.push_back(interpPoints.at(1));
	}

//Histograms______________________________________________________________________________________
	//Time
	int nbins = 300;
	double xlow = 0;
	double xup = 0.3;

	TH1D *hTimeB = new TH1D("Time","Timing; time [ms]; Number of Events", nbins, xlow, xup); 
	hTimeB->SetStats(0);
	TH1D *hTimeC = new TH1D("TimeC","Timing; time [ms]; Number of Events", nbins, xlow, xup); 
	hTimeC->SetLineColor(kRed);
	hTimeC->SetStats(0);

	for(int i=0; i<(int)timec.size(); i++) 
	{
		hTimeB->Fill(timeb.at(i));
		hTimeC->Fill(timec.at(i));	
	}

	//Interpolation
	vector <double> interpB, interpC;
	for(int i = 0; i < (int)yEvents.size(); i++)
	{
		for(int j = 0; j < (int)yEvents[i].size(); j++)
		{
			int indexForB = binarySearch(xBSplineValues[i], xEvents[i][j]);
			int indexForC = binarySearch(xCSplineValues[i], xEvents[i][j]);

			interpB.push_back( (yEvents[i][j]-yBSplineValues[i][indexForB])/yErrorEvents[i][j] );
			interpC.push_back( (yEvents[i][j]-yCSplineValues[i][indexForC])/yErrorEvents[i][j] );
		}
	}

	int nbinsI = 101;
	double xlowI = -0.1;
	double xupI = 0.1;
	TH1D *hInterpB = new TH1D("Interp B","Interpolation; Distance between spline and data normalized by error; Number of Events", nbinsI, xlowI, xupI); 
	for(int i=0; i<(int)interpB.size(); i++) hInterpB->Fill(interpB.at(i));
	hInterpB->SetStats(0);

	TH1D *hInterpC = new TH1D("Interp C","Interpolation; Distance between spline and data normalized by error; Number of Events", nbinsI, xlowI, xupI); 
	for (int i=0; i<(int)interpC.size(); i++) hInterpC->Fill(interpC.at(i));
	hInterpC->SetLineColor(kGreen);
	hInterpC->SetStats(0);	

	//Differences
	int nbinsd = 1000;
	double xlowd = 0;
	double xupd = 100;

	TH1D *hDiffB = new TH1D("Diff B","Differences; difference between spline and linear interpolation; Number of Events", nbinsd, xlowd, xupd); 
	hDiffB->SetStats(0);
	TH1D *hDiffC = new TH1D("Diff C","Differences; difference between spline and linear interpolation; Number of Events", nbinsd, xlowd, xupd); 
	hDiffC->SetLineColor(kRed);
	hDiffC->SetStats(0);

	vector<double> diffB = diff(yBSplineValues, yLinearInterpValues);
	vector<double> diffC = diff(yCSplineValues, yLinearInterpValues);	

	for(int i=0; i<(int)diffB.size(); i++) 
	{
		hDiffB->Fill(diffB.at(i));
		hDiffC->Fill(diffC.at(i));	
	}

//Test Graphs for splines_____________________________________________________________________

	TGraph *GCspline = new TGraph(xCSplineValues[iEventLook].size(), &xCSplineValues[iEventLook][0], &yCSplineValues[iEventLook][0]);
	GCspline->SetLineColor(kRed);
	TGraph *GBspline = new TGraph(xBSplineValues[iEventLook].size(), &xBSplineValues[iEventLook][0], &yBSplineValues[iEventLook][0]);
	TGraph *Gdata = new TGraph(xEvents[0].size(), &xEvents[iEventLook][0], &yEvents[iEventLook][0]);
	Gdata->SetMarkerStyle(20);
	Gdata->GetHistogram()->SetMinimum(-5);
	Gdata->GetHistogram()->SetMaximum(25);
	TGraph *lin = new TGraph((int)xLinearInterpValues[iEventLook].size(), &xLinearInterpValues[iEventLook][0], &yLinearInterpValues[iEventLook][0]);
	lin->SetLineColor(kBlue);
	lin->SetMarkerStyle(2);

//Draws______________________________________________________________________________________
	//Interpolation 
	TLegend *legInterp = new TLegend(0.9,0.70,0.75,0.85);
	legInterp->SetLineColor(kWhite); 
	legInterp->SetFillColor(kWhite);
	legInterp->SetMargin(0.3); 
	legInterp->AddEntry(hInterpB,"b-spline","l");
	legInterp->AddEntry(hInterpC,"c-spline","l");
	legInterp->SetTextSize(0.05);

	TCanvas *c1 = new TCanvas("c1", "Interpolation distance");
	c1->cd();
	hInterpB->Draw("");
	hInterpC->Draw("same");
	legInterp->Draw();

	//Time
	TLegend *legTime = new TLegend(0.9,0.70,0.75,0.85);
	legTime->SetLineColor(kWhite); 
	legTime->SetFillColor(kWhite);
	legTime->SetMargin(0.3); 
	legTime->AddEntry(hTimeB,"b-spline","l");
	legTime->AddEntry(hTimeC,"c-spline","l");
	legTime->SetTextSize(0.05);

	TCanvas *c2 = new TCanvas("c2", "Computation time");
	c2->cd();
	hTimeB->Draw("");
	hTimeC->Draw("same");
	legTime-> Draw();

	//Differences
	TLegend *legDiff = new TLegend(0.9,0.70,0.75,0.85);
	legDiff->SetLineColor(kWhite); 
	legDiff->SetFillColor(kWhite);
	legDiff->SetMargin(0.3); 
	legDiff->AddEntry(hDiffB,"b-spline","l");
	legDiff->AddEntry(hDiffC,"c-spline","l");
	legDiff->SetTextSize(0.05);

	TCanvas *c3 = new TCanvas("c3", "Differences");
	c3->cd()->SetLogx();
	c3->cd()->SetLogy();
	hDiffC->Draw("");
	hDiffB->Draw("same");
	legDiff-> Draw();

	TCanvas *c4 = new TCanvas("c4", "Test splines");
	c4->cd();
	Gdata->Draw("ap");
	GCspline->Draw("samel");
	GBspline->Draw("samel");
	lin->Draw("samel");

	//Free the memory used
	gsl_spline_free (spline_GLOB); 
 	gsl_interp_accel_free (acc_GLOB);
	gsl_bspline_free(bw);
	gsl_vector_free(xControl);
	gsl_vector_free(yControl);
	gsl_vector_free(w);
	gsl_vector_free(B);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_multifit_linear_free(mw);
	gsl_matrix_free(cov);
}

