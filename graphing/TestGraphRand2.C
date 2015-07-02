#include "TCanvas.h"
#include "TGraphErrors.h"
#include <vector>
#include <iostream>

// I need to include something for the TRandom class if I want it to complie, but I'm not sure what to include.


using namespace std;

void TestGraphRand2() {// named with the title so it will run automatically

   TCanvas *c1 = new TCanvas("c1","My Awesome Test Graph!!",200,10,700,500);

   c1->SetFillColor(70);
  c1->SetGrid();
  c1->GetFrame()->SetFillColor(100);
  c1->GetFrame()->SetBorderSize(1000);

   const int n = 10;

TRandom3 jrand;// make an object jrand of the type TRandom3, this is the same syntax as saying int jrand, execpt I don't want an int. 

vector<double> xVector;
vector<double> yVector;
vector<double> xErrorVector;
vector<double> yErrorVector;

for(int i =0; i<n; i=i+1){
xVector.push_back(jrand.Integer(10));
yVector.push_back(jrand.Integer(10));
xErrorVector.push_back(jrand.Integer(10));
yErrorVector.push_back(jrand.Integer(10));
}// here I'm using my for loop to fill my vectors with random numbers instead of the parts of an array. I found a function Integer(), which in combo with a random number, makes a random number between either 1 or 0 and, in this case, 10. 



   TGraphErrors *gr = new TGraphErrors(n,&xVector[0],&yVector[0],&xErrorVector[0],&yErrorVector[0]);
 
   gr->SetTitle("Test Graph With Random Data and Error Bars");
   gr->SetMarkerColor(2);
   gr->SetMarkerStyle(21);
   gr->Draw("ALP");


   c1->Update();

}
