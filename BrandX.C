#include "TCanvas.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include "Math/Polynomial.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"

// How to plot an interpolator? -- the ROOT way, with a global variable.
ROOT::Math::Interpolator *gInter;

double interp_wrap(double *x, double *par) {
   double s;
   s = gInter->Eval(x[0]);
   return s;
}


// Sample code to exercise some of the spline interpolator
// functionality already present in Root.  This can act as a
// comparison against the more specialized approaches.

// Kinds of test data sets to prepare
enum tfun_t {trig, sort_rand};

// Documented in the gnu GSL manual, section 27.3

/*
  kLINEAR, kPOLYNOMIAL, kCSPLINE, kCSPLINE_PERIODIC, kAKIMA,
  kAKIMA_PERIODIC
*/

void fill_test_data(vector<double> &xv, vector<double> &yv, tfun_t t){
  switch (t)
    {
    case sort_rand:
      // Fill with random data and sort.  Gives a bumpy straight line
      {
	TRandom3 jrand(123454321);
	const int n = 8;
	for(int i =0; i<n; i=i+1) {
	  xv.push_back(jrand.Uniform(20.0));
	  yv.push_back(jrand.Uniform(6.5, 36.5));
	}
	sort(xv.begin(), xv.end());
	sort(yv.begin(), yv.end());
      }
      return;
    case trig:
      {
	for (int i = 0; i<8; i++){
	  double x = i*0.55+.3;
	  xv.push_back(x);
	  yv.push_back(sin(x));
	}
      }
      return;
    }
}

// Utility for testing.
void print_vector(vector<double> myvec) {
    cout << "[ ";
    for (vector<double>::iterator it = myvec.begin(); it != myvec.end(); it++) {
        cout << *it << ", ";   // NO endl
    }
    cout << "]" << endl;
}


void BrandX(){
  vector<double> xv;
  vector<double> yv;

  fill_test_data(xv, yv, sort_rand);

  print_vector(xv);
  print_vector(yv);

  Double_t w = 600;  // These coordinates are in pixels
  Double_t h = 600;
  TCanvas *c1 = new TCanvas("c hello", "Test data and interpolation", w, h);
  c1->SetWindowSize(w + (w - c1->GetWw()), h + (h - c1->GetWh()));

  c1->SetGrid();
  TGraph* gr = new TGraph(xv.size(),&xv[0],&yv[0]);

  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("AP");

  // Create and interpolator, and superimpose the curve on top of the
  // plotted data points.
  
  // Here are the choices.  These are documented in the Gnu GSL
  // manual, section 27.3
  /**
     enum Type { kLINEAR, kPOLYNOMIAL, kCSPLINE, kCSPLINE_PERIODIC,
         kAKIMA, kAKIMA_PERIODIC };
  */

  ROOT::Math::Interpolation::Type itype;
  itype = ROOT::Math::Interpolation::kCSPLINE;

  ROOT::Math::Interpolator *inter = 
    new ROOT::Math::Interpolator(xv, yv, itype);
  
  gInter = inter;
    
  TF1* gr2 =  new TF1("interp", interp_wrap, xv[0], *(xv.end()-1), 0);
  gr2->Draw("SAME");

}
	
