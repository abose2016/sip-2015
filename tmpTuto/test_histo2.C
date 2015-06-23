#include "TH1D.h"
#include "TCanvas.h"

void DrawHisto(){

	int nbins = 2;
	double xmin = 0;
	double xmax = 1;
	TH1D *h0 = new TH1D("h0","Title of my histogram",nbins,xmin,xmax);
	h0->SetBinContent(1,5.4);
	h0->SetBinContent(2,3.1);

	TCanvas *c1 = new TCanvas("c1","Title of my canvas");
	c1->cd();
	h0->Draw();
}
