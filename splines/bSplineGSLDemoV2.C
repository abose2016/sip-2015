#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TAxis.h"
#include <vector>
#include "TRandom3.h"

void bSplineGSLDemoV2 ()
{
	const size_t nControl = 15;
	int nControlInt = 15;
	const size_t ncoeffs = 12;
	const size_t nbreak = ncoeffs - 2;
	size_t i, j;
	gsl_bspline_workspace *bw;
	gsl_vector *B;
	gsl_vector *c, *w;
	gsl_vector *xControl, *yControl;
	gsl_matrix *X, *cov;
	gsl_multifit_linear_workspace *mw;
	double chisq, Rsq, dof, tss;
	double splineStep = 0.1;
	int nValues = nControlInt/splineStep;
	double xOrigin[nControlInt], yOrigin[nControlInt], xValues[nValues], yValues[nValues];
 
	/* allocate a cubic bspline workspace (k = 4) */
	bw = gsl_bspline_alloc(4, nbreak);
	B = gsl_vector_alloc(ncoeffs);

	xControl = gsl_vector_alloc(nControl);
	yControl = gsl_vector_alloc(nControl);
	X = gsl_matrix_alloc(nControl, ncoeffs);
	c = gsl_vector_alloc(ncoeffs);
	w = gsl_vector_alloc(nControl);
	cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
	mw = gsl_multifit_linear_alloc(nControl, ncoeffs);
	TRandom3 jrand;

	for (i = 0; i < nControl; ++i)
			 {
				double sigma;
				double xi = (nControlInt / (nControlInt - 1)) * i;
				double yi = jrand.Uniform(20);
				sigma = 0.1 * yi;
				gsl_vector_set(xControl, i, xi);
				xOrigin[i] = xi;
				gsl_vector_set(yControl, i, yi);
				yOrigin[i] = yi;
				gsl_vector_set(w, i, 1.0 / (sigma * sigma));

				 printf("%f %f\n", xi, yi);
			 }

	
	/* use uniform breakpoints on [0, 15] */
	gsl_bspline_knots_uniform(0.0, nControlInt, bw);

	/* construct the fit matrix X */
	for (i = 0; i < nControl; ++i)
	{
		double xi = gsl_vector_get(xControl, i);

		/* compute B_j(xi) for all j */
		gsl_bspline_eval(xi, B, bw);

		/* fill in row i of X */
		for (j = 0; j < ncoeffs; ++j)
		{
			double Bj = gsl_vector_get(B, j);
			gsl_matrix_set(X, i, j, Bj);
		}
	}

	/* do the fit */
	gsl_multifit_wlinear(X, w, yControl, c, cov, &chisq, mw);

	dof = nControl - ncoeffs;
	tss = gsl_stats_wtss(w->data, 1, yControl->data, 1, yControl->size);
	Rsq = 1.0 - chisq / tss;

	fprintf(stderr, "chisq/dof = %e, Rsq = %f\n", 
		                chisq / dof, Rsq);

	/* output the smoothed curve */
	
		double xi, yi, yerr;
		int nSpline = int((xOrigin[nControlInt-1]-xOrigin[0])/splineStep);
		printf("BREAK\n\n\n\n\n");	
		
		for (xi = 0.0; xi < nControl; xi += splineStep)
		{
			gsl_bspline_eval(xi, B, bw);
			gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
			int index = (int)xi*10;
			xValues[index] = xi;
			yValues[index] = yi;

			printf("%f %f\n", xi, yi);	
		}
	
	gsl_bspline_free(bw);
	gsl_vector_free(B);
	gsl_vector_free(xControl);
	gsl_vector_free(yControl);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_vector_free(w);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(mw);

    TGraph *grControlPoints = new TGraph(nControlInt, xOrigin, yOrigin);
    grControlPoints->SetMarkerStyle(20);

    TGraph *grSpline = new TGraph(nSpline, xValues, yValues);
    grSpline->SetTitle("");
    grSpline->GetXaxis()->SetTitle("X-axis  [A.U.]");
    grSpline->GetYaxis()->SetTitle("Y-axis  [A.U.]");

    //plot
    TCanvas *c1 = new TCanvas("c1", "Graph", 200, 10, 700, 500);
    c1->cd();
    grSpline->Draw("al");
    grControlPoints->Draw("same p");
} 
