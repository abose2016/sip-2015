#include <iostream>
#include <TMinuit.h>
#include <TApplication.h>
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"

//Globally defined function
TF1 *f2 = new TF1("f2","(sin(x))",-10,10);

//The function to be minimised
void fcn(int& npar, double* deriv, double& f, double par[], int flag){
  f=f2->Eval(par[0],par[1]);
}   
     
int minuitExample(){
  //First Draw function
  f2->Draw("LEGO2");
  cout<<"Program finished.."<<endl;
  
  //Setup the minimisation procedure
  const int npar = 2;              // the number of parameters
  TMinuit minuit(npar);
  minuit.SetFCN(fcn);

  double par[npar];               // the start values
  double stepSize[npar];          // step sizes 
  double minVal[npar];            // minimum bound on parameter 
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];

  par[0] = 100.0;      // a guess
  stepSize[0] = 0.1; // take e.g. 0.1 of start value
  minVal[0] = -10;   // if min and max values = 0, parameter is unbounded.
  maxVal[0] = 10;
  parName[0] = "x";

  par[1] = 100.0;      // a guess
  stepSize[1] = 0.1; // take e.g. 0.1 of start value
  minVal[1] = -10;   // if min and max values = 0, parameter is unbounded.
  maxVal[1] = 10;
  parName[1] = "y";

  //Setup parameters
  for (int i=0; i<npar; i++){
    minuit.DefineParameter(i, parName[i].c_str(), 
			   par[i], stepSize[i], minVal[i], maxVal[i]);
  }
  
  cout<<"Running Migrad()...."<<endl;
  // Do the minimization!
  minuit.Migrad();       // Minuit's best minimization algorithm
  cout<<"Migrad() completed..."<<endl;
 
  //Get the Minuit results
  double outpar[npar], err[npar];
  for (int i=0; i<npar; i++){
      minuit.GetParameter(i,outpar[i],err[i]);
      cout<<"Fitted parameter:"<<i<<" is:"<<outpar[i]<<" +/- "<<err[i]<<endl;
  }

  cout<<"Program Finished"<<endl;

 return 0;
}