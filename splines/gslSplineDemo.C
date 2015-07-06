#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TAxis.h"

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

void gslSplineDemo ()
{
	int i;
	double xi, yi, x[10], y[10];
	vector<double> xValues, yValues;

	printf ("#m=0,S=2\n");
	for (i = 0; i < 10; i++)
	{
		x[i] = i + 0.5 * sin (i);
		y[i] = i + cos (i * i);
		printf ("%g %g\n", x[i], y[i]);
	}
		printf ("#m=1,S=0\n");
	{
		gsl_interp_accel *acc = gsl_interp_accel_alloc ();
		gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, 10);
		gsl_spline_init (spline, x, y, 10);
		for (xi = x[0]; xi < x[9]; xi += 0.01)
		{
			yi = gsl_spline_eval (spline, xi, acc);
			xValues.push_back(xi);
			yValues.push_back(yi);
		}
	
		gsl_spline_free (spline);
		gsl_interp_accel_free (acc);
	}

	TGraph *gr = LoadGraphFromVectors(xValues, yValues);
	TCanvas *c1 = new TCanvas("c1", "Graph", 200, 10, 700, 500);
	gr->Draw("apz");

	c1->Update();
}
