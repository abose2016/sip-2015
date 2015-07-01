#ifndef __CINT__
#include "Riostream.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TGraphErrors.h"
#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "TF1.h"
#endif

using namespace std;

void FillRandVectors(vector<double> &xVector, int n, int seed = 250, double lowerBound= 5, double upperBound= 20)
{
	//Call TRandom3
	TRandom3 *jrand = new TRandom3(seed);
	for (int i = 0; i<n; i = i + 1)
	{
		xVector.push_back(jrand->Uniform(lowerBound, upperBound));
		
	}
	delete jrand;
}


void matrixMake(){
	Double_t det1;
  TMatrixD a(4,4);
	TMatrixD b(4,4);
	vector <double> data;
	FillRandVectors(data, 16);
	int counter=0;


	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{
			
			a[i][j]= data.at(counter);
			b[i][j]= data.at(counter);
			counter++;

		}
	} 


	for(int i=0; i<4; i++)
		{
			for(int j=0; j<4; j++)
			{
				cout<< a[i][j]<< ",";
			}
		cout<< "\n";
		} 


	a.InvertFast(&det1);
	cout << "After inversion" <<endl;

	for(int i=0; i<4; i++)
		{
			for(int j=0; j<4; j++)
			{
				std::cout<< a[i][j]<< ",";
			}
		std::cout<< "\n";
		} 

	
}
