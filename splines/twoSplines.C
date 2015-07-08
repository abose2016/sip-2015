#include <gsl/gsl_spline.h>
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

const int nControl = 15;
const int ncoeffs = 12;
const int nbreak = ncoeffs-2;

TGraph *bSpline( double stepSpline, std::vector <double> xControlBasis, std::vector <double> yControlBasis)
{

	//Declare and allocate memory to compose data set of control points
	gsl_vector *xControl, *yControl, *w;
	vector<double> xOrigin, yOrigin;
	xControl = gsl_vector_alloc(nControl);
	yControl = gsl_vector_alloc(nControl);
	w = gsl_vector_alloc(nControl);
	double sigma, xi, yi;

	//Populate data set with monotonically increasing, uniform x values and randomized y values
	for (int i = 0; i < nControl; ++i)
		{

			xi = xControlBasis.at(i);
			yi = yControlBasis.at(i);
			sigma = 0.1 * yi;
			gsl_vector_set(xControl, i, xi);
			gsl_vector_set(yControl, i, yi);
			gsl_vector_set(w, i, 1.0 / (sigma * sigma));
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
	X = gsl_matrix_alloc(nControl, ncoeffs);
	cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
	for (int i = 0; i < nControl; ++i)
	{
		xi = gsl_vector_get(xControl, i);

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
	double chisq, Rsq, dof, tss;
	gsl_vector *c;
	gsl_multifit_linear_workspace *mw;
	c = gsl_vector_alloc(ncoeffs);
	mw = gsl_multifit_linear_alloc(nControl, ncoeffs);

	//Do the fit
	gsl_multifit_wlinear(X, w, yControl, c, cov, &chisq, mw);
	
	dof = nControl - ncoeffs;
	tss = gsl_stats_wtss(w->data, 1, yControl->data, 1, yControl->size);
	Rsq = 1.0 - chisq / tss;

	std::cout<<stderr << "chisq/dof "<< chisq/dof << ", Rsq = " << Rsq << std::endl;

	//Output the curve and store the values of the spline in two vectors
	double yerr;
	vector<double> xValues, yValues;
	for (xi = 0.0; xi < 15.0; xi += stepSpline)
	{
		gsl_bspline_eval(xi, B, bw);
		gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
		xValues.push_back(xi);
		yValues.push_back(yi);

		std::cout<< xi<< "   " << yi << std::endl;
	}

	//Free the memory used
	gsl_bspline_free(bw);
	gsl_vector_free(B);
	gsl_vector_free(xControl);
	gsl_vector_free(yControl);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_vector_free(w);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(mw);

	//Load graphs
	int numSplinePoints = xValues.size();
	TGraph *grSpline = new TGraph(numSplinePoints, &xValues[0], &yValues[0]);
	return grSpline;
}

TGraph *naturalCubic( double stepSpline, vector <double> xControlCubic, vector <double> yControlCubic)
{
    //initialize the spline
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, nControl);

    //compute the spline
    gsl_spline_init (spline, &xControlCubic[0], &yControlCubic[0], nControl);

    //evaluate the spline
    int nSpline = int((xControlCubic.at(nControl-1)-xControlCubic.at(0))/stepSpline);
    vector <double> xSpline, ySpline;
    for (int i= 0; i < nSpline; i++)
    {
        xSpline.push_back(xControlCubic.at(0)+i*stepSpline);
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

	//Initialize variables
	vector <double> xControl, yControl;
	double xi, yi;
	TRandom3 *jrand= new TRandom3(seed);

	//Populate data set with monotonically increasing, uniform x values and randomized y values
	for (int i = 0; i < nControl; ++i)
		{
			xi = (15.0 / (nControl - 1)) * i;
		 	yi = jrand->Uniform(20);
			xControl.push_back(xi);
			yControl.push_back(yi);
			std::cout<< xControl.at(i) << "   " << yControl.at(i)<< std::endl;
		 }

	TGraph *bGraph = bSpline(stepSpline, xControl, yControl);
	bGraph->SetTitle("");
  bGraph->GetXaxis()->SetTitle("X-axis  [A.U.]");
  bGraph->GetYaxis()->SetTitle("Y-axis  [A.U.]");
	TGraph *cGraph = naturalCubic(stepSpline, xControl, yControl);
	cGraph->SetLineColor(kRed);
	TGraph *pGraph = new TGraph (nControl, &xControl[0], &yControl[0]);
	pGraph-> SetMarkerStyle(20);
	pGraph->SetMarkerColor(kBlue);

//Draw to canvas
	TCanvas *c1 = new TCanvas("c1", "Graph", 200, 10, 700, 500);
 	c1->cd();
	bGraph->Draw("al");
	cGraph->Draw("SAME l");
	pGraph->Draw("SAME p");
}
