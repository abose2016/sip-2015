#include <math.h>
#include <vector>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TRandom3.h"

void bSplineGSLDemoV4 (int seed = 7898, double stepSpline = 0.01)
{
	//Initialize variables
	const int nControl = 15;
	const int ncoeffs =  nControl;
	const int orderSpline = 4;
	const int nbreak = ncoeffs+2-orderSpline;

	double xmin = 0;
	double xmax = 15;

	//Declare and allocate memory to compose data set of control points
	gsl_vector *xControl = gsl_vector_alloc(nControl);
	gsl_vector *yControl = gsl_vector_alloc(nControl);

	//Populate data set with monotonically increasing, uniform x values and randomized y values
	TRandom3 *jrand = new TRandom3(seed);
	vector<double> xPlot, yPlot;
	for (int i = 0; i < nControl; ++i)
	{
		double xi = xmin + i*(xmax-xmin)/(nControl-1.);
		double yi = jrand->Uniform(20);
		gsl_vector_set(xControl, i, xi);
		gsl_vector_set(yControl, i, yi);
		xPlot.push_back(xi);
		yPlot.push_back(yi);
	}

	//Create a b spline workspace and allocate its memory
	gsl_bspline_workspace *bw = gsl_bspline_alloc(4, nbreak);
	gsl_bspline_knots_uniform(xmin, xmax, bw);

	//Set up the variables for the fit matrix
	gsl_vector *B = gsl_vector_alloc(ncoeffs);
	gsl_matrix *X = gsl_matrix_alloc(nControl, ncoeffs);
	for (int i = 0; i < nControl; ++i)
	{
		gsl_bspline_eval(gsl_vector_get(xControl, i), B, bw);//Compute B_j(xi) for all j 
		for (int j = 0; j < ncoeffs; ++j)	gsl_matrix_set(X, i, j, gsl_vector_get(B, j));//Fill in row i of X
	}

	//Declare variables for the fit and allocate their memory
	double chisq;	
	gsl_vector *c = gsl_vector_alloc(ncoeffs);
	gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(nControl, ncoeffs);

	//Do the fit
	gsl_matrix *cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
	gsl_multifit_linear(X, yControl, c, cov, &chisq, mw);

	//Output the curve and store the values of the spline in two vectors
	int nValues = 1+int((xPlot.back() - xPlot.front())/stepSpline);
	vector<double> xValues, yValues;
	for (int i = 0; i < nValues; i++)
	{
		double xi = xPlot.front()+i*stepSpline;
		gsl_bspline_eval(xi, B, bw);
		double max;
		for(int j = 0; j < (B->size - 1.0); j++) 
		{
			double next = gsl_vector_get(B, j+1);
			double curr = gsl_vector_get(B, j);
			if(next > curr) 
				max = next;
			else
				max = curr;
		}
		xValues.push_back(xi);
		yValues.push_back(max);
	}

	//Free the memory used
	gsl_vector_free(xControl);
	gsl_vector_free(yControl);
	gsl_bspline_free(bw);
	gsl_vector_free(B);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_multifit_linear_free(mw);
	gsl_matrix_free(cov);

	//Load graphs
	TGraph *grControlPoints = new TGraph(xPlot.size(), &xPlot[0], &yPlot[0]);
 	grControlPoints->SetMarkerStyle(20);

	TGraph *grSpline = new TGraph(xValues.size(), &xValues[0], &yValues[0]);
	grSpline->SetTitle("");
	grSpline->GetXaxis()->SetTitle("X-axis  [A.U.]");
	grSpline->GetYaxis()->SetTitle("Y-axis  [A.U.]");

	//Draw to canvas
	TCanvas *c1 = new TCanvas("c1", "Graph", 200, 10, 700, 500);
 	c1->cd();
	grSpline->Draw("al");
	grControlPoints->Draw("same p");
} 
