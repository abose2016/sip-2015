#include "TCanvas.h"
#include "TGraphErrors.h"
#include <vector>
#include <iostream>
#include "TRandom.h"
using namespace std;

vector<TString> polynomials;
vector<double> xVector, yVector, xErrorVector, yErrorVector, newPlotY, newPlotX, tempX, tempY;
const int n = 7;

void basicFitWithErrors() {
	TCanvas *c1 = new TCanvas("c1", "interpolation", 0, 0, 1000, 800); 

	FillRandVectors(xVector, yVector, xErrorVector, yErrorVector, n);
	TGraphErrors *g1 = LoadGraphFromVectors(xVector, yVector, xErrorVector, yErrorVector);

	polynomials.push_back("pol0");
	polynomials.push_back("pol1");
	polynomials.push_back("pol2");
	polynomials.push_back("pol3");
	polynomials.push_back("pol4");
	polynomials.push_back("pol5");
	polynomials.push_back("pol6");
	polynomials.push_back("pol7");
	polynomials.push_back("pol8");

	gStyle->SetOptFit(1111);

	for (int i = 0; i < polynomials.size(); i++) {
		TString curr = polynomials.at(i);

		TF1 *fa1 = new TF1("fa1", curr, 0, 10);

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

	TCanvas *c2 = new TCanvas("c2", "Chi square", 0, 0, 1000, 800);
	TGraph *g2 = LoadGraphFromVectorsWithUnits(newPlotX, newPlotY);
	g2->Draw("ap");
	c2-> Update();
}

void FillRandVectors(vector<double> &xVector, vector< double > &yVector, vector< double > &xErrorVector, vector< double > &yErrorVector, int n)
	{
		//Call TRandom3

		int seed = 250;
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

TGraphErrors *LoadGraphFromVectors(std::vector< double > xVector, std::vector< double > yVector, std::vector< double > xErrorVector, std::vector< double > yErrorVector)
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
		gr->GetXaxis()->SetTitle("X axis [Arbitrary Units]");
		gr->GetXaxis()->CenterTitle();
		gr->GetYaxis()->SetTitle("Y axis [Arbitrary Units]");
		gr->GetYaxis()->CenterTitle();
		return gr;
	}
	else
	{
		TGraphErrors *gr0 = new TGraphErrors();
		return gr0;
	}
}

TGraph *LoadGraphFromVectorsWithUnits(std::vector< double > xVector, std::vector< double > yVector)
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
		gr->GetXaxis()->SetTitle("Number of parameters");
		gr->GetXaxis()->CenterTitle();
		gr->GetYaxis()->SetTitle("Chi square");
		gr->GetYaxis()->CenterTitle();
		return gr;
	}
	else
	{
		TGraph *gr0 = new TGraph();
		return gr0;
	}
}