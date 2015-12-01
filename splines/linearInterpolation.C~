#include <vector>
#include <iostream>
#include <iomanip>
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include <string>
#include "TRandom3.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TF1.h"

using namespace std;

//_____________________________________________________________________________
// Fill y and error vectors with random points from TRandom3
void FillRandVectors(int nPoints, vector<double> &xVector, vector< double > &yVector, int seed = 250, double lowerBound= 10, double upperBound= 20)
{
	//Call TRandom3
	TRandom3 *jrand = new TRandom3(seed);
	xVector.clear();
	yVector.clear();

	for (int i = 0; i < nPoints; i++)
	{
		xVector.push_back(i);
		yVector.push_back(jrand->Uniform(lowerBound, upperBound));
	}
	delete jrand;
}
//_____________________________________________________________________________

vector <vector<double> > linearInterpolation( vector<double> xData, vector<double> yData, double stepSpline = 0.1)
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
void grapher()
{
	vector<double> xData, yData;
	int nPoints = 8;
	FillRandVectors(nPoints, xData, yData);
	vector <vector<double> > interpPoints = linearInterpolation(xData, yData);

	TGraph *lin = new TGraph(nPoints, &interpPoints[0][0], &interpPoints[1][0]);
	lin->SetMarkerColor(kBlue);
	lin->SetMarkerStyle(20);

	TGraph *data = new TGraph(nPoints, &xData[0], &yData[0]);
	data->SetMarkerColor(kRed);
	data->SetMarkerStyle(20);

	TCanvas *c1 = new TCanvas("c1", "Linear");
	c1->cd();
//	lin->Draw("ap");
	data->Draw("ap"); 
}

	
	
