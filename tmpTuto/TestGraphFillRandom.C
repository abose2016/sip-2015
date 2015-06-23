#include "TCanvas.h"
#include "TGraphErrors.h"
#include <vector>
#include <iostream>
// I need to include something for the TRandom class

using namespace std;

void TestGraphFillRandom() {

   TCanvas *c1 = new TCanvas("c1","My Awesome Test Graph!!",200,10,700,500);

   c1->SetFillColor(70);
  c1->SetGrid();
  c1->GetFrame()->SetFillColor(100);
  c1->GetFrame()->SetBorderSize(1000);

   const int n = 10;
double xarray[n], yarray[n], exarray[n], eyarray[n];// empty arrays of size n. 

TRandom3 jrand;// TRandom3 is the best so I'm making an object of the type TRandom3. 

jrand.RndmArray(n, xarray);// RndmArray(size of array, name of array) fills the arrays with random numbers.
jrand.RndmArray(n, yarray);
jrand.RndmArray(n, exarray);
jrand.RndmArray(n, eyarray);


vector<double> xVector;
vector<double> yVector;
vector<double> xErrorVector;
vector<double> yErrorVector;

for(int i =0; i<n; i=i+1){
xVector.push_back(xarray[i]);
yVector.push_back(yarray[i]);
xErrorVector.push_back(exarray[i]);
yErrorVector.push_back(eyarray[i]);
}// this for loop just fills the vectors with the values in the arrays




   TGraphErrors *gr = new TGraphErrors(n,&xVector[0],&yVector[0],&xErrorVector[0],&yErrorVector[0]);
 
   gr->SetTitle("Test Graph With Random Data and Error Bars");
   gr->SetMarkerColor(2);
   gr->SetMarkerStyle(21);
   gr->Draw("ALP");


   c1->Update();

}
 
