#include "TFrame.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TAxis.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include "TSpline.h"

using namespace std;

//put random numbers into some vectors
void FillRandVectors(vector< double > &xVector,
                     vector< double > &yVector,
					 int n)
{
    //Call TRandom3

    int seed = 23382;
    TRandom3 *jrand = new TRandom3(seed);
    for(int i =0; i<n; i=i+1)
    {
        xVector.push_back(jrand->Uniform(20.0));
        yVector.push_back(jrand->Uniform(5.0, 20.0));
    }
	TGraph *gr = LoadGraphFromVectors(xVector, yVector);

    delete jrand;
}

// make graph from vectors
TGraph *LoadGraphFromVectors(std::vector< double > xVector, std::vector< double > yVector)
{
    int size = xVector.size();

    if(xVector.size()== yVector.size())
    {
        //Create a graph
		TGraph *gr = new TGraph(size, &xVector[0], &yVector[0]);
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

		return null;
    }
}

//Main of the program
void randSpline()
{
    // making vectors
    vector <double> xVector, yVector, xErrorVector, yErrorVector;
    const int n = 6;

    // filling vectors with random numbers
    FillRandVectors(xVector, yVector, n);
    // making graph
	TGraph *gr = LoadGraphFromVectors(xVector, yVector);

    //Plot
    TCanvas *c1 = new TCanvas("c1","My Awesome Test Graph!!",200,10,700,500);
    gr->Draw("apz");

// Draw the Cubic Spline to the same canvas
   TSpline3 *sp = new TSpline3("Cubic Spline", gr,"b2e2", 0, 0);
   
   sp->SetLineColor(kRed);
   sp->Draw("lsame");

    c1->Update();
}
