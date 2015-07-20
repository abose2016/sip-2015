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
#include "TH1.h"
#include "TLegend.h"
#include <sstream>

using namespace std;

//Data definition
gsl_interp_accel *acc_GLOB;
gsl_spline *spline_GLOB;
gsl_bspline_workspace *bw_GLOB;

vector<double> xData_GLOB, yData_GLOB, yErrorData_GLOB;
//_______________________________________________________________________________
TH1D *histograms( vector <double> data, int nbins = 8, int xlow = 0, int xup = 8) 
{
	//make histogram
	TH1D *hello1 = new TH1D("title","title; xAxis; yAxis", nbins, xlow, xup); 
	for ( int i=0; i<nbins; i++) hello1->Fill(i,data.at(i));

	return hello1;
} 

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
 TGraph *cSpline(int nPoints, int npar, vector <double> xData, vector <double> yData, vector <double> yErrorData, double stepSpline =.01, double start = 10., double step = 0.01)
{
	//Populate the global variables
	xData_GLOB = xData;
	yData_GLOB = yData;
	yErrorData_GLOB = yErrorData;

	//Initialize Minuit
	vector<double> vstart, vstep;
	for(int i=0; i<npar; i++) 	//set starting values and step sizes for parameters
	{
		vstart.push_back(start);
		vstep.push_back(step);
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
	gsl_spline_init (spline_GLOB, &xData[0], &bestFitParams[0], xData.size()); //initialize the spline

	int nPointsSpline = int ((xData[nPoints-1]-xData[0])/stepSpline); //calculate the number of points of the spline
	vector< double > xSpline, ySpline;
	for (int i= 0; i < nPointsSpline; i++){
		xSpline.push_back(xData[0]+i*stepSpline);
		ySpline.push_back(gsl_spline_eval (spline_GLOB, xSpline.back(), acc_GLOB));
	}

	
	// Graph the spline
	TGraph *grSplineC = new TGraph(nPointsSpline, &xSpline[0], &ySpline[0]);
	grSplineC->SetTitle("");
	grSplineC->GetXaxis()->SetTitle("X-axis  [A.U.]");
	grSplineC->GetYaxis()->SetTitle("Y-axis  [A.U.]");
	
	return grSplineC;

}
//_________________________________________________________________________
TGraph *bSpline(int nControl, int npar, vector <double> xDataB, vector <double> yDataB, vector <double> yErrorDataB, double stepSpline = 0.01, double xmin = 0, double xmax = 9)
{
	//Declare and allocate memory to compose data set of control points
	gsl_vector *xControl = gsl_vector_alloc(nControl);
	gsl_vector *yControl = gsl_vector_alloc(nControl);
	gsl_vector *w = gsl_vector_alloc(nControl);

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
	gsl_vector *B = gsl_vector_alloc(npar);
	gsl_matrix *X = gsl_matrix_alloc(nControl, npar);
	for (int i = 0; i < nControl; ++i)
	{
		gsl_bspline_eval(gsl_vector_get(xControl, i), B, bw_GLOB);//Compute B_j(xi) for all j 
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
	gsl_vector_free(B);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_multifit_linear_free(mw);
	gsl_matrix_free(cov);


/*	//histograms
	//interpolation
	int nbinsI = 8, xlowI = 0, xupI = 8;
	vector <double> interp;
	for(int i=1; i< nControl; i++) interp.push_back(yDataB.at(i)/gsl_vector_get(yControl,i));
	TH1D *helloInterp = new TH1D("Interp","Interp; xAxis; yAxis", nbinsI, xlowI, xupI); 
	for ( int i=0; i<nControl-1; i++) helloInterp->Fill(interp.at(i),1);
//	TCanvas *canInterp = new TCanvas("c5", "Histogram", 400, 20, 1000, 1000);
//	canInterp->cd();
//	helloInterp->Draw("");
 
	//linearity
	int nbinsL = 8, xlowL = 0, xupL = 8;
	vector <double> linear;
	for(int i=1; i< nControl; i++) linear.push_back(splineError.at(i)/yErrorDataB.at(i));
	TH1D *helloInterp = new TH1D("Interp","Interp; xAxis; yAxis", nbinsL, xlowL, xupL); 
	for ( int i=0; i<nControl-1; i++) helloInterp->Fill(linear.at(i),1);
//	TCanvas *canInterp = new TCanvas("c5", "Histogram", 400, 20, 1000, 1000);
//	canInterp->cd();
//	helloInterp->Draw("");
*/


	for(int i=1; i< nControl; i++)
	{
		cout<<yDataB.at(i)/gsl_vector_get(yControl,i)<<"   "<<yDataB.at(i)<<"   "<<gsl_vector_get(yControl,i)<< endl;
	}

	for(int i=1; i< nControl; i++)
	{
		cout<<splineError.at(i)<<"   "<<yErrorDataB.at(i)<<"   "<<splineError.at(i)/yErrorDataB.at(i)<< endl;
	}

	TGraph *grSplineB = new TGraph(xValues.size(), &xValues[0], &yValues[0]);
	return grSplineB;
} 

//_________________________________________________________________________________
void integratedSplinesV2a(double seed = 231) 
{
	//Load the data
	int nPoints = 9;
	vector <double> xData, yData, yErrorData; 
	FillRandVectors(nPoints, xData, yData, yErrorData, seed); //create random vectors for y values and y error vector

	//Intialization of the variables
	const int npar = nPoints;
	const int orderSpline = 4;
	const int nbreak = npar+2-orderSpline;
	double stepSpline = 0.01;
	double xminBSplineWorkspace = 0;
	double xmaxBSplineWorkspace = 9;

	acc_GLOB = gsl_interp_accel_alloc ();
	spline_GLOB = gsl_spline_alloc (gsl_interp_cspline, nPoints);	
	bw_GLOB = gsl_bspline_alloc(orderSpline, nbreak);
	
	//B-spline
	clock_t tbstart = clock();
	TGraph *bGraph = bSpline(nPoints, npar, xData, yData, yErrorData, stepSpline, xminBSplineWorkspace, xmaxBSplineWorkspace);
	bGraph->SetTitle("");
	bGraph->GetXaxis()->SetTitle("X-axis Arbitrary Units");
	bGraph->GetYaxis()->SetTitle("Y-axis Arbirary Units");
	clock_t tbstop = clock();


	//C-spline
	clock_t tcstart = clock();
	TGraph *cGraph = cSpline(nPoints, npar, xData, yData, yErrorData, stepSpline);
	clock_t tcstop = clock();
	cGraph->SetLineColor(kRed);


	//Control points
	TGraphErrors *pGraph = new TGraphErrors(nPoints, &xData[0], &yData[0], 0 ,&yErrorData[0]);
	pGraph-> SetMarkerStyle(20);
	pGraph->SetMarkerColor(kBlue);

	//Free the memory for the spline
	gsl_spline_free (spline_GLOB); //frees the memory used by the spline
 	gsl_interp_accel_free (acc_GLOB);
	gsl_bspline_free(bw_GLOB);

	// time histogram
	int nbins = 100;
	double xlow = 0;
	double xup = 0.001;

	TH1D *hello2 = new TH1D("TimeB","TimeB; time; number of runs", nbins, xlow, xup); 
	double timeb = ((float)tbstop-(float)tbstart)/ CLOCKS_PER_SEC;
	hello2->Fill(timeb);

	TH1D *hello3 = new TH1D("TimeC","TimeC; time; number of runs", nbins, xlow, xup); 
	hello3->SetLineColor(kRed);
	double timec = ((float)tcstop-(float)tcstart)/ CLOCKS_PER_SEC;
	hello3->Fill(timec);

	//Legends
	TLegend *leg = new TLegend(0.15,0.70,0.4,0.85);
	leg->SetLineColor(kWhite); 
	leg->SetFillColor(kWhite);
	leg->SetMargin(0.3); 
	leg->AddEntry(pGraph,"Simulated data","lp");
	leg->AddEntry(bGraph,"b-spline","l");
	leg->AddEntry(cGraph,"c-spline","l");


	//Draw to canvas
	TCanvas *can1 = new TCanvas("c2", "Timing");
	can1->cd();
	hello2->Draw("");
	hello3->Draw("same");

	TCanvas *c1 = new TCanvas("c1", "Spline Comparison");
 	c1->cd();
	bGraph->Draw("al");
	cGraph->Draw("same l");
	pGraph->Draw("same pz");
	leg->Draw();
}

