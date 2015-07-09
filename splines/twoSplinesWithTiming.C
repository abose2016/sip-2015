#include <gsl/gsl_spline.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TRandom3.h"
#include <time.h>

//Basis spline
TGraph *bSpline( double stepSpline, std::vector <double> xControlBasis, std::vector <double> yControlBasis)
{
	int nControl = xControlBasis.size();
	int nParams = nControl;//can be tuned with nParams<=nControl
	int orderSpline = 4;
	int nBreak = nParams+2-orderSpline;

	//Declare and allocate memory to compose data set of control points
	gsl_vector *xControl = gsl_vector_alloc(nControl);
	gsl_vector *yControl = gsl_vector_alloc(nControl);
	gsl_vector *w = gsl_vector_alloc(nControl);

	//Populate data set with monotonically increasing, uniform x values and randomized y values
	for (int i = 0; i < nControl; ++i)
	{
		double xi = xControlBasis.at(i);
		double yi = yControlBasis.at(i);
		gsl_vector_set(xControl, i, xi);
		gsl_vector_set(yControl, i, yi);
		gsl_vector_set(w, i, 1.);
	}

	//Create a b spline workspace and allocate its memory
	gsl_bspline_workspace *bw = gsl_bspline_alloc(orderSpline, nBreak);
	gsl_bspline_knots_uniform(xControlBasis.front(), xControlBasis.back(), bw);//uniform breakpoints

	//Set up the variables for the fit matrix
	gsl_vector *B = gsl_vector_alloc(nParams);
	gsl_matrix *X = gsl_matrix_alloc(nControl, nParams);//fit matrix X
	for (int i = 0; i < nControl; ++i)
	{
		gsl_bspline_eval(gsl_vector_get(xControl, i), B, bw);//Compute B_j(xi) for all j 
		for (int j = 0; j < nParams; ++j) gsl_matrix_set(X, i, j, gsl_vector_get(B, j));//Fill in row i of X
	}

	//Declare variables for the fit and allocate their memory
	gsl_vector *c = gsl_vector_alloc(nParams);
	gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(nControl, nParams);

	//Do the fit
	double chisq;
	gsl_matrix *cov = gsl_matrix_alloc(nParams, nParams);

	gsl_multifit_wlinear(X, w, yControl, c, cov, &chisq, mw);

	//Output the curve and store the values of the spline in two vectors
	vector<double> xValues, yValues;
	int nValues = 1+int((xControlBasis.back() - xControlBasis.front())/stepSpline);
	for (int i = 0; i < nValues; i++)
	{
		double xi = xControlBasis.front()+i*stepSpline;
		gsl_bspline_eval(xi, B, bw);
		double yerr, yi;
		gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
		xValues.push_back(xi);
		yValues.push_back(yi);
	}

	//Free the memory used
	gsl_vector_free(xControl);
	gsl_vector_free(yControl);
	gsl_vector_free(w);
	gsl_bspline_free(bw);
	gsl_vector_free(B);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_multifit_linear_free(mw);
	gsl_matrix_free(cov);

	//Load graphs
	TGraph *grSpline = new TGraph(xValues.size(), &xValues[0], &yValues[0]);
	return grSpline;
}

//Cubic
TGraph *naturalCubic(double stepSpline, vector <double> xControlCubic, vector <double> yControlCubic)
{
		int nControl = xControlCubic.size();

    //initialize the spline
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, nControl);

    //compute the spline
    gsl_spline_init (spline, &xControlCubic[0], &yControlCubic[0], nControl);

	  int nSpline = 1+int((xControlCubic.back()-xControlCubic.front())/stepSpline);
    vector <double> xSpline, ySpline;
	  //evaluate the spline
	  for (int i= 0; i < nSpline; i++)
	  {
	      xSpline.push_back(xControlCubic.front()+i*stepSpline);
	      ySpline.push_back(gsl_spline_eval (spline, xSpline.at(i), acc));
	  }

    //clear the spline
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    //load the graphs
    TGraph *grSpline = new TGraph(nSpline,&xSpline[0],&ySpline[0]);

    return grSpline;
}

//main of the program
void twoSplines(int seed = 7898, double stepSpline =.01)
{
	int nControl = 15;

	//Initialize variables
	TRandom3 *jrand= new TRandom3(seed);
	vector <double> xControl, yControl;
	for (int i = 0; i < nControl; ++i)
	{
		double xi = ((double)nControl / (nControl - 1)) * i;
	 	double yi = jrand->Uniform(20);
		xControl.push_back(xi);
		yControl.push_back(yi);
	}


	//B-spline
	clock_t tbstart = clock();
	TGraph *bGraph = bSpline(stepSpline, xControl, yControl);
	clock_t tbstop = clock();
	bGraph->SetTitle("");
  bGraph->GetXaxis()->SetTitle("X-axis  [A.U.]");
  bGraph->GetYaxis()->SetTitle("Y-axis  [A.U.]");

	//C-spline
	clock_t tcstart = clock();
	TGraph *cGraph = naturalCubic(stepSpline, xControl, yControl);
	clock_t tcstop = clock();
	cGraph->SetLineColor(kRed);

	//Control points
	TGraph *pGraph = new TGraph (nControl, &xControl[0], &yControl[0]);
	pGraph-> SetMarkerStyle(20);
	pGraph->SetMarkerColor(kBlue);

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
