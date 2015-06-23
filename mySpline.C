#include "TSpline.h"


// making some arrays
//const int n = 4;// number of points
//double Point_x[n] = {0, 1, 1.5, 2};
//double Point_y[n] = {1.43, 0.368, 0.135, 0.018};

 const int n = 10;
double xarray[n], yarray[n];// empty arrays of size n. 

 int seed = 546563;
    TRandom3 *jrand = new TRandom3(seed);// TRandom3 is the best so I'm making an object of the type TRandom3. 

jrand.RndmArray(n, xarray);// RndmArray(size of array, name of array) fills the arrays with random numbers.
jrand.RndmArray(n, yarray);


// main function with the same name of this file
// it will be executed when you type ".x interpolation.C" in ROOT
void mySpline()
{
// Section 1. Draw the points on a canvas
   TCanvas *c1 = new TCanvas("c1","interpolation",0,0,1000,800);

   TGraph *g1 = new TGraph(n, xarray, yarray);// just plots points

   g1->SetMarkerStyle(20);
   g1->SetMarkerSize(2);
   
   g1->GetXaxis()->SetLimits(-1, 3);        // set real range
   g1->GetXaxis()->SetRangeUser(-0.5, 2.5); // set visible range
   g1->GetXaxis()->SetTitle("X");
   g1->GetXaxis()->CenterTitle();

   g1->GetYaxis()->SetLimits(-1, 2.0);
   g1->GetYaxis()->SetRangeUser(-0.4, 1.6);
   g1->GetYaxis()->SetTitle("Y");
   g1->GetYaxis()->CenterTitle();

   g1->Draw("ap"); // options to draw a graph are described on 
                   // http://root.cern.ch/root/html/TGraph.html#TGraph:PaintGraph

// Section 3. Draw the Cubic Spline to the same canvas
   TSpline3 *sp = new TSpline3("Cubic Spline", xarray, yarray, n, "b2e2", 0, 0);
   // refer to http://root.cern.ch/root/html/Tspline3.html for the usage of TSpline3
   // "b2e2" together with the last two "0" means that the second derivatives 
   // of the begin and end points equal to zero
   sp->SetLineColor(kRed);
   sp->Draw("lsame");


//____________________________________________________________________________
double LogChi2(double *x, double *par) 
{
   double ChiSqr = 0;
   for (int i=0; i<n; i++)
      ChiSqr = ChiSqr + (yarray[i] - x[0]*exp(-x[1]*xarray[i]))**2;
   return log(ChiSqr);
}}


