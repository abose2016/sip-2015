#include "TCanvas.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TAxis.h"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

//put random numbers into some vectors
void FillRandVectors(vector< double > &xVector,
                     vector< double > &yVector,
                     vector< double > &xErrorVector, 	 			     vector< double > &yErrorVector, int n)
{
    //Call TRandom3

    int seed = 546563;
    TRandom3 *jrand = new TRandom3(seed);
    for(int i =0; i<n; i=i+1)
    {
        xVector.push_back(jrand->Uniform(20.0));
        yVector.push_back(jrand->Uniform(5.0, 20.0));
        xErrorVector.push_back(jrand->Uniform(5.0));
        yErrorVector.push_back(jrand->Uniform(5.0));
    }
    TGraphErrors *gr = LoadGraphFromVectors(xVector,yVector,xErrorVector,yErrorVector);

    delete jrand;
}

// make graph from vectors
TGraphErrors *LoadGraphFromVectors(std::vector< double > xVector, std::vector< double > yVector, std::vector< double > xErrorVector, std::vector< double > yErrorVector)
{
    int n = xVector.size();

    if((xVector.size()== yVector.size())&&
            (yVector.size()== xErrorVector.size())&&
            (xErrorVector.size()== yErrorVector.size()))
    {
        //Create a graph
        TGraphErrors *gr = new TGraphErrors(n,&xVector[0],
                                            &yVector[0], &xErrorVector[0],
                                            &yErrorVector[0]);
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

//Main of the program
void TwoPlots()
{

    // making vectors
    vector <double> xVector, yVector, xErrorVector,
           yErrorVector;
    const int n = 10;

    // filling vectors with random numbers
    FillRandVectors(xVector, yVector, xErrorVector,
                    yErrorVector, n);
    // making graph
    TGraphErrors *gr = LoadGraphFromVectors(xVector, yVector, xErrorVector,
                                            yErrorVector);
   
    //Plot
    TCanvas *c1 = new TCanvas("c1","My Awesome Test Graph!!",200,10,700,500);
    gr->Draw("apz");

    // Draw a function on top
    TF1 * jFun= new TF1("new function","10+100*sin(x)/x", 0.0,20.0);
    jFun-> Draw("SAME");

    c1->Update();
}




