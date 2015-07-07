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

/* number of data points to fit */
#define N 15

/* number of fit coefficients */
#define NCOEFFS 12

/* nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 */
#define NBREAK (NCOEFFS - 2)

// make graph from vectors
TGraph *LoadGraphFromVectors(std::vector< double > xVector, std::vector< double > yVector)
{
    int n = xVector.size();

    if((xVector.size()== yVector.size()))
    {
        //Create a graph
        TGraph *gr = new TGraph(n, &xVector[0], &yVector[0]);
        gr->SetTitle("");
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(.2);
        gr->SetLineWidth(100);
        gr->GetXaxis()->SetTitle("X axis [Arbitrary Units]");
        gr->GetXaxis()->CenterTitle();
        gr->GetYaxis()->SetTitle("Y axis [Arbitrary Units]");
        gr->GetYaxis()->CenterTitle();
        return gr;
    }
    else
    {
        TGraph *gr0 = new TGraph();
        return gr0;
    }
}

void bSplineGSLDemo ()
{
	const size_t n = N;
	const size_t ncoeffs = NCOEFFS;
	const size_t nbreak = NBREAK;
	size_t i, j;
	gsl_bspline_workspace *bw;
	gsl_vector *B;
	gsl_rng *r; 
	gsl_vector *c, *w;
	gsl_vector *xControl, *yControl;
	gsl_matrix *X, *cov;
	gsl_multifit_linear_workspace *mw;
	double chisq, Rsq, dof, tss;
	vector<double> xValues, yValues, xOrigin, yOrigin;
 
	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);

	/* allocate a cubic bspline workspace (k = 4) */
	bw = gsl_bspline_alloc(4, nbreak);
	B = gsl_vector_alloc(ncoeffs);

	xControl = gsl_vector_alloc(n);
	yControl = gsl_vector_alloc(n);
	X = gsl_matrix_alloc(n, ncoeffs);
	c = gsl_vector_alloc(ncoeffs);
	w = gsl_vector_alloc(n);
	cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
	mw = gsl_multifit_linear_alloc(n, ncoeffs);
	TRandom3 jrand;

	printf("#m=0,S=0\n");
	for (i = 0; i < n; ++i)
			 {
				double sigma;
				double xi = (15.0 / (N - 1)) * i;
				double yi = jrand.Uniform(20);
				sigma = 0.1 * yi;
				 gsl_vector_set(xControl, i, xi);
				xOrigin.push_back(xi);
				 gsl_vector_set(yControl, i, yi);
				yOrigin.push_back(yi);
				gsl_vector_set(w, i, 1.0 / (sigma * sigma));

				 printf("%f %f\n", xi, yi);
			 }

	
	/* use uniform breakpoints on [0, 15] */
	gsl_bspline_knots_uniform(0.0, 15.0, bw);

	/* construct the fit matrix X */
	for (i = 0; i < n; ++i)
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

	dof = n - ncoeffs;
	tss = gsl_stats_wtss(w->data, 1, yControl->data, 1, yControl->size);
	Rsq = 1.0 - chisq / tss;

	fprintf(stderr, "chisq/dof = %e, Rsq = %f\n", 
		                chisq / dof, Rsq);

	/* output the smoothed curve */
	{
		double xi, yi, yerr;
		printf("BREAK\n\n\n\n\n");	
		
		printf("#m=1,S=0\n");
		for (xi = 0.0; xi < 15.0; xi += 0.1)
		{
			gsl_bspline_eval(xi, B, bw);
			gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
			xValues.push_back(xi);
			yValues.push_back(yi);

			printf("%f %f\n", xi, yi);	
		}
	}

	gsl_rng_free(r);
	gsl_bspline_free(bw);
	gsl_vector_free(B);
	gsl_vector_free(xControl);
	gsl_vector_free(yControl);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_vector_free(w);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(mw);

	TGraph *gr = LoadGraphFromVectors(xValues, yValues);
	TGraph *gr1 = LoadGraphFromVectors(xOrigin, yOrigin);
	TCanvas *c1 = new TCanvas("c1", "Graph", 200, 10, 700, 500);
	gr->Draw("al");
	gr1->Draw("SAME p");

	c1->Update();
} 
