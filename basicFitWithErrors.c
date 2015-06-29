#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include <vector>
#include <string>
#include <iostream>
#include "TRandom3.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TF1.h"
#include "TStyle.h"

using namespace std;



//fill x,y and error vectors with random points from TRandom#
void FillRandVectors(vector<double> &xVector, vector< double > &yVector, vector< double > &xErrorVector, vector< double > &yErrorVector, int n, int seed = 250)
{
	//Call TRandom3
	TRandom3 *jrand = new TRandom3(seed);
	for (int i = 0; i<n; i = i + 1)
	{
		xVector.push_back(jrand->Uniform(20.0));
		yVector.push_back(jrand->Uniform(5.0, 20.0));
		xErrorVector.push_back(jrand->Uniform(0.5, 5.0));
		yErrorVector.push_back(jrand->Uniform(0.5, 5.0));
	}
	delete jrand;
}

//Construct a graph from vectors and error vectors
TGraphErrors *LoadGraphFromVectorsWithError(std::vector< double > xVector, std::vector< double > yVector, std::vector< double > xErrorVector, std::vector< double > yErrorVector, string xTitle, string yTitle)
{
	int n = xVector.size();

	if ((xVector.size() == yVector.size()) &&
		(yVector.size() == xErrorVector.size()) &&
		(xErrorVector.size() == yErrorVector.size()))
	{
		//Create a graph
		TGraphErrors *gr = new TGraphErrors(n, &xVector[0], &yVector[0], &xErrorVector[0], &yErrorVector[0]);
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

// construct a graph without errors
TGraph *LoadGraphFromVectors(std::vector< double > xVector, std::vector< double > yVector, string xTitle, string yTitle )
{
	int n = xVector.size();

	if ((xVector.size() == yVector.size()))
	{
		//Create a graph
		TGraph *gr = new TGraph(n, &xVector[0], &yVector[0]);
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
		TGraph *gr0 = new TGraph();
		return gr0;
		delete gr0;
	}
}


//main of the program
void basicFitWithErrors() {

	// setting variables
	const int n = 7;
	vector<string> polynomials;
	vector<double> xVector, yVector, xErrorVector, yErrorVector, newPlotY, 	 	newPlotX, tempX, tempY;

	
	FillRandVectors(xVector, yVector, xErrorVector, yErrorVector, n, 807340);
	
	polynomials.push_back("pol0");
	polynomials.push_back("pol1");
	polynomials.push_back("pol2");
	polynomials.push_back("pol3");
	polynomials.push_back("pol4");
	polynomials.push_back("pol5");
	polynomials.push_back("pol6");
	polynomials.push_back("pol7");
	polynomials.push_back("pol8");

	
	TGraphErrors *g1 = LoadGraphFromVectorsWithError(xVector, yVector, 	 	xErrorVector, yErrorVector, "X Axis (arbitrary units)", "Y Axis  	 	(arbitrary units)");
	gStyle->SetOptFit(1111);
	TCanvas *c1 = new TCanvas("c1", "interpolation", 0, 0, 1000, 800); 




	for (int i = 0; i < (int)polynomials.size(); i++) {
		string curr = polynomials.at(i);

		TF1 *fa1 = new TF1("fa1", curr.c_str(), 0, 10);

		g1->Draw("ap");
		g1->Fit(fa1);
		double chi2 = fa1->GetChisquare();
		int nParInt = fa1->GetNpar();
		double nPar = (double)nParInt;

		newPlotX.push_back(nPar);
		newPlotY.push_back(chi2);

		c1->Update();
		gSystem->Sleep(1000);
	}

	TGraph *g2 = LoadGraphFromVectors(newPlotX, newPlotY, "Number of Parameters", "Chi Squared");
	TCanvas *c2 = new TCanvas("c2", "Chi square", 0, 0, 1000, 800);
	g2->Draw("ap");
	c2-> Update();
}

