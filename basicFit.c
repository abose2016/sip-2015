const Long64_t n = 4;
Double_t Point_x[n] = { 0, 1, 1.5, 2 };
Double_t Point_y[n] = { 1.43, 0.368, 0.135, 0.018 };
vector<TString> polynomials;


void basicFit() {
	// Section 1. Draw the points on a canvas
	TCanvas *c1 = new TCanvas("c1", "interpolation", 0, 0, 1000, 800);

	TGraph *g1 = new TGraph(n, Point_x, Point_y);

	polynomials.push_back("pol0");
	polynomials.push_back("pol1");
	polynomials.push_back("pol2");
	polynomials.push_back("pol3");
	polynomials.push_back("pol4");
	polynomials.push_back("pol5");
	polynomials.push_back("pol6");


	for (int i = 0; i < polynomials.size(); i++) {
		TString curr = polynomials.at(i);

		TF1 *fa1 = new TF1("fa1", curr, 0, 10);

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
		g1->Fit(fa1);

		c1->Update();
		gSystem->Sleep(1000);
	}
}