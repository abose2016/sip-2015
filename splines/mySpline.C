#include "TSpline.h"
#define ARRAY_SIZE(array) (sizeof((array))/sizeof((array[0])))

const int size = 5;
double xarray[5], yarray[5];// empty arrays of size "size". 

int seed = 546563;
TRandom3 *jrand = new TRandom3(seed);// TRandom3 is the best so I'm making an object of the type TRandom3. 

jrand->RndmArray(5, xarray);// RndmArray(size of array, name of array) fills the arrays with random numbers.
jrand->RndmArray(5, yarray);

// main function with the same name of this file
// it will be executed when you type ".x interpolation.C" in ROOT
void mySpline()
{
// Section 1. Draw the points on a canvas
	cout << "Size of x array: " << ARRAY_SIZE(xarray) << endl;  // prints 5 as expected
	cout << "Size of y array: " << ARRAY_SIZE(yarray) << endl;  // prints 5 as expected

   TCanvas *c1 = new TCanvas("c1","interpolation",0,0,1000,800);

   TGraph *g1 = new TGraph(5, xarray, yarray);// just plots points

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
   TSpline3 *sp = new TSpline3("Cubic Spline", xarray, yarray, 5, "b2e2", 0, 0);
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


