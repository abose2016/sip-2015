#include "TCanvas.h"
#include "TGraphErrors.h"
#include <vector>
#include <iostream>
#include <algorithm> //I think I need this header to use sort because sort is a subfunction of algorithm
// I need to include something for the TRandom class if I want it to complie, but I'm not sure what to include.

using namespace std;

void TestGraphRandBetter() {// named with the title so it will run automatically

   TCanvas *c1 = new TCanvas("c1","My Awesome Test Graph!!",200,10,700,500);

   c1->SetFillColor(70);
  c1->SetGrid();
  c1->GetFrame()->SetFillColor(100);
  c1->GetFrame()->SetBorderSize(1000);

   const int n = 10;

TRandom3 jrand;

vector<double> xVector;
vector<double> yVector;
vector<double> xErrorVector;
vector<double> yErrorVector;

for(int i =0; i<n; i=i+1){
xVector.push_back(jrand.Integer(20));
yVector.push_back(jrand.Integer(20));
xErrorVector.push_back(jrand.Integer(5));
yErrorVector.push_back(jrand.Integer(5));
}// I changed my parameters so that the errors would be out of a smaller number than the x and y values. 

sort(xVector.begin(), xVector.end());
sort(yVector.begin(), yVector.end());
// This puts my vectors in order so that when I graph it, the points will be plotted from least to greatest. 

   TGraphErrors *gr = new TGraphErrors(n,&xVector[0],&yVector[0],&xErrorVector[0],&yErrorVector[0]);
 
   gr->SetTitle("Test Graph With Random Data and Error Bars");
   gr->SetMarkerColor(2);
   gr->SetMarkerStyle(21);
   gr->Draw("ALP");


   c1->Update();

}
