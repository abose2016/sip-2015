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

void bSplineGSLDemoV3 (int seed = 7898, double stepSpline = 0.01)
{
	//Initialize variables
	const int n = 15;
	const int ncoeffs = 12;
	const int nbreak = ncoeffs-2;

	//Declare and allocate memory to compose data set of control points
	gsl_vector *xControl, *yControl;
	vector<double> xOrigin, yOrigin;
	xControl = gsl_vector_alloc(n);
	yControl = gsl_vector_alloc(n);
	TRandom3 *jrand = new TRandom3(seed);

	//Populate data set with monotonically increasing, uniform x values and randomized y values
	for (int i = 0; i < n; ++i)
		{
			double sigma;
			double xi = (15.0 / (n - 1)) * i;
			double yi = jrand->Uniform(20);
			sigma = 0.1 * yi;
			gsl_vector_set(xControl, i, xi);
			xOrigin.push_back(xi);
			gsl_vector_set(yControl, i, yi);
			yOrigin.push_back(yi);
			std::cout << xi << "   " << yi << std::endl;
		 }

	//Create a b spline workspace and allocate its memory
	gsl_bspline_workspace *bw;
	bw = gsl_bspline_alloc(4, nbreak);

	//Use uniform breakpoints on [0, 15]
	gsl_bspline_knots_uniform(0.0, 15.0, bw);

	//Set up the variables for the fit matrix
	gsl_vector *B;
	B = gsl_vector_alloc(ncoeffs);
	
	//Construct the fit matrix X
	gsl_matrix *X, *cov;
	X = gsl_matrix_alloc(n, ncoeffs);
	cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
	for (int i = 0; i < n; ++i)
	{
		double xi = gsl_vector_get(xControl, i);

		//Compute B_j(xi) for all j 
		gsl_bspline_eval(xi, B, bw);

		//Fill in row i of X
		for (int j = 0; j < ncoeffs; ++j)
		{
			double Bj = gsl_vector_get(B, j);
			gsl_matrix_set(X, i, j, Bj);
		}
	}

	//Declare variables for the fit and allocate their memory
	double chisq;	
	gsl_vector *c;
	gsl_multifit_linear_workspace *mw;
	c = gsl_vector_alloc(ncoeffs);
	mw = gsl_multifit_linear_alloc(n, ncoeffs);

	//Do the fit
	gsl_multifit_linear(X, yControl, c, cov, &chisq, mw);
	

	//Output the curve and store the values of the spline in two vectors
	double xi, yi, yerr;
	vector<double> xValues, yValues;
	int index = 0;
	for (xi = 0.0; xi < 15.0; xi += stepSpline)
	{
		gsl_bspline_eval(xi, B, bw);
		gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
		xValues.push_back(xi);
//		yi = gsl_vector_get(B, index);
		yValues.push_back(yi);

		std::cout<< xi<< "   " << yi << std::endl;
		index++;
	}

	//Free the memory used
	gsl_bspline_free(bw);
	gsl_vector_free(B);
	gsl_vector_free(xControl);
	gsl_vector_free(yControl);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(mw);

	//Load graphs
	int numSplinePoints = xValues.size(), numControlPoints = xOrigin.size();
	TGraph *grControlPoints = new TGraph(numControlPoints, &xOrigin[0], &yOrigin[0]);
 	grControlPoints->SetMarkerStyle(20);

	TGraph *grSpline = new TGraph(numSplinePoints, &xValues[0], &yValues[0]);
	grSpline->SetTitle("");
	grSpline->GetXaxis()->SetTitle("X-axis  [A.U.]");
	grSpline->GetYaxis()->SetTitle("Y-axis  [A.U.]");

	//Draw to canvas
	TCanvas *c1 = new TCanvas("c1", "Graph", 200, 10, 700, 500);
 	c1->cd();
	grSpline->Draw("al");
	grControlPoints->Draw("SAME p");
} 
