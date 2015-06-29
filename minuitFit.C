#include <vector>
#include <iostream>
#include "TMinuit.h"

using namespace std;


//Data definition
const int nbins = 5;
vector <double> x,y,xError;

//______________________________________________________________________________
double func(float x,double *par)
{
 double value=( ((par[0]*x*x)+ (par[1]*x)+ par[2]));
 return value;
}

//______________________________________________________________________________
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
   //calculate chisquare
   double chisq = 0;
   for (int i=0;i<nbins; i++) {
     double delta  = (y[i]-func(x[i],par))/xError[i];
     chisq += delta*delta;
   }
   f = chisq;
}

//______________________________________________________________________________
void minuitFit()
{
	// x values
	x.push_back(1);
	x.push_back(2.5);
	x.push_back(7.67);
	x.push_back(8);
	x.push_back(9.1);
	// y values
	y.push_back(10);
	y.push_back(3.5);
	y.push_back(8.67);
	y.push_back(25);
	y.push_back(2);
	// xError
	xError.push_back(.2);
	xError.push_back(.3);
	xError.push_back(.8);
	xError.push_back(.1);
	xError.push_back(1);


   int npar = 4;
   TMinuit *gMinuit = new TMinuit(npar);  //initialize TMinuit with a maximum of 5 params
   gMinuit->SetFCN(fcn);

   double arglist[10];
   int ierflg = 0;
   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

// Set starting values and step sizes for parameters
   static vector <double> vstart, step;
	vstart.push_back(3);
	vstart.push_back(1);
	vstart.push_back(.1);
	vstart.push_back(.01);
	
	step.push_back(.1);
	step.push_back(.1);
	step.push_back(.01);
	step.push_back(.001);

   gMinuit->mnparm(0, "a1", vstart.at(0), step.at(0), 0,0,ierflg);
   gMinuit->mnparm(1, "a2", vstart.at(1), step.at(1), 0,0,ierflg);
   gMinuit->mnparm(2, "a3", vstart.at(2), step.at(2), 0,0,ierflg);
   gMinuit->mnparm(3, "a4", vstart.at(3), step.at(3), 0,0,ierflg);

// Now ready for minimization step
   arglist[0] = 500;
   arglist[1] = 1.;
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

// Print results
   double amin,edm,errdef;
   int nvpar,nparx,icstat;
   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   //gMinuit->mnprin(3,amin);
}

