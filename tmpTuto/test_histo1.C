{
	int nbins = 2;
	double xmin = 0;
	double xmax = 1;
	TH1D *h0 = new TH1D("h0","Title of my histogram",nbins,xmin,xmax);

	h0->SetBinContent(1,5.4);
	h0->SetBinContent(2,3.1);

	h0->Draw();
}
