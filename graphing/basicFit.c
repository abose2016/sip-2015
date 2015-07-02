#include "TCanvas.h"
#include "TGraphErrors.h"
#include <vector>
#include <iostream>
#include "TRandom.h"
using namespace std;

vector<TString> polynomials;
vector<double> xVector, yVector, newPlotY, newPlotX;
const int n = 7;

void basicFit() {
	TCanvas *c1 = new TCanvas("c1", "interpolation", 0, 0, 1000, 800);

	FillRandVectors(xVector, yVector, n);
	TGraph *g1 = LoadGraphFromVectors(xVector, yVector);

	polynomials.push_back("pol0");
	polynomials.push_back("pol1");
	polynomials.push_back("pol2");
	polynomials.push_back("pol3");
	polynomials.push_back("pol4");
	polynomials.push_back("pol5");
	polynomials.push_back("pol6");

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
	c2->Update();
}

void FillRandVectors(vector<double> &xVector, vector< double > &yVector, int n)
{
	//Call TRandom3

	int seed = 7898;
	TRandom3 *jrand = new TRandom3(seed);
	for (int i = 0; i<n; i = i + 1)
	{
		xVector.push_back(jrand->Uniform(20.0));
		yVector.push_back(jrand->Uniform(5.0, 20.0));
	}
	TGraphErrors *gr = LoadGraphFromVectors(xVector, yVector);

	delete jrand;
}

TGraph *LoadGraphFromVectors(std::vector< double > xVector, std::vector< double > yVector)
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