{
  gSystem->Load("libTree.so");
  gSystem->Load("libHist.so");
  gSystem->Load("libGpad.so");
  gSystem->Load("/usr/lib/libgsl.so");
  gSystem->Load("/usr/lib/libgslcblas.so");
  gSystem->SetMakeSharedLib("cd $BuildDir ; g++ -c $Opt -pipe -Wall -W -Woverloaded-virtual -fPIC -O3 -g -Iinclude -pthread $IncludePath $SourceFiles ;  g++ $ObjectFiles -shared -Wl,-soname,$LibName.so -O  -g -o $SharedLib");

	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetTitleXOffset(1.2);
	gStyle->SetTitleYOffset(1.1);

	gStyle->SetPadLeftMargin(0.1);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.05);

	gStyle->SetTextSize(0.045);
	gStyle->SetTitleSize(0.045,"xyz");
	gStyle->SetLabelSize(0.045,"xyz");

	gStyle->SetTickLength(0.02,"Y");
	gStyle->SetTickLength(0.02,"X");
}



