#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TRandom3.h"

// make graph from vectors
TGraph *LoadGraphFromVectors(std::vector<double> xVector, std::vector<double> yVector)
{

	int n= xVector.size();
	if(xVector.size()== yVector.size())
	{
		//Create a graph
		TGraph *gr = new TGraph(n, &xVector[0], &yVector[0]);
		gr->SetTitle("");
		gr->SetMarkerStyle(20);
		gr->SetMarkerSize(1);
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


//main of the program
void gslSplineDemoV2(double nstep =.01)
{

	//iterator
	int i;
	//number of control points
	const int n=10;
	//initialize data arrays
	double xi, yi, x[n]= {1,2,3,4,5,6,7,8,9,10}, y[n];
	//initialize data vectors
	std::vector<double> xValues, yValues, xControlPoints, yControlPoints;
	// make a random array
	TRandom3 jrand;
	jrand.RndmArray(n,y);
	
	//Fill control point vectors for graphing
	for ( i=0; i<n ; i++)
	{
		xControlPoints.push_back(x[i]);
		yControlPoints.push_back(y[i]);
	}

	// this is like initializing a vector, except we're initializing a spline
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);

	// here we make the spline, this is the line where all the important stuff happens
	gsl_spline_init (spline, x, y, n);

	// evaluate the spline  for each x value and fill x and y vectors
	for (xi = x[0]; xi < x[9]; xi += nstep)
	{
		yi = gsl_spline_eval (spline, xi, acc);
		xValues.push_back(xi);
		yValues.push_back(yi);
	}
	
	//freeing the memory space, equivalent to deleting variables
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);

	//plot data
	TGraph *gr1 = LoadGraphFromVectors(xValues, yValues);
	TGraph *gr2 = LoadGraphFromVectors(xControlPoints,yControlPoints);
	TCanvas *c1 = new TCanvas("c1", "Graph", 200, 10, 700, 500);
	gr1->Draw("al");
	gr2->Draw("SAME p");
	c1->Update();
}
