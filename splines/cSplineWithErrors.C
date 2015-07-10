#include <gsl/gsl_spline.h>
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TRandom3.h"

//main of the program
void gslSplineDemoV3(double stepSpline =.01)
{

    //initialize data arrays
    const int nControl=10;
    double xControl[nControl], yControl[nControl], yError[nControl];
		int seed = 7898;
		TRandom3 *jrand = new TRandom3(seed);
		for (int i=0; i<nControl; i++) 
		{
		xControl[i]=i; 
		yControl[i]= jrand->Uniform(10,20); 
		yError[i]= jrand->Uniform(1,2);
		}

    //initialize the spline
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, nControl);

    //compute the spline
    gsl_spline_init (spline, xControl, yControl, nControl);

    //evaluate the spline
    int nSpline = int((xControl[nControl-1]-xControl[0])/stepSpline);
    double xSpline[nSpline], ySpline[nSpline];
    for (int i= 0; i < nSpline; i++)
    {
        xSpline[i] = xControl[0]+i*stepSpline;
        ySpline[i] = gsl_spline_eval (spline, xSpline[i], acc);
    }
    
    //clear the spline
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    //load the graphs
    TGraphErrors *grControlPoints = new TGraphErrors(nControl,xControl,yControl,0,yError);
    grControlPoints->SetMarkerStyle(20);

    TGraph *grSpline = new TGraph(nSpline,xSpline,ySpline);
    grSpline->SetTitle("");
    grSpline->GetXaxis()->SetTitle("X-axis  [A.U.]");
    grSpline->GetYaxis()->SetTitle("Y-axis  [A.U.]");

    //plot
    TCanvas *c1 = new TCanvas("c1", "Graph", 200, 10, 700, 500);
    c1->cd();
    grSpline->Draw("al");
    grControlPoints->Draw("same p");
}
