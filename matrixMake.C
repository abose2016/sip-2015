#include <iostream>
#include <vector>
#include "TMatrixD.h"
#include "TMath.h"
#include "TRandom3.h"

//Fill x,y and error vectors with random points from TRandom#
void FillRandVectors(std::vector< double > &xVector, int n, int seed = 250, double lowerBound = 5, double upperBound = 20)
{
	//Call TRandom3
	TRandom3 *jrand = new TRandom3(seed);
	for (int i = 0; i < n; i = i + 1)
	{
		xVector.push_back(jrand->Uniform(lowerBound, upperBound));
	}
	delete jrand;
}

//Main function - creates randomly populated matrix A, inverts it, and then multiplies the two to create identity matrix
void matrixMake(){
	int counter = 0;
	int matrixSize = 3;

	//Loading random values of matrix into vector labeled data
	std::vector< double > data;
	FillRandVectors(data, matrixSize*matrixSize);

	//Creating two square matrices
	TMatrixD a(matrixSize, matrixSize);
	TMatrixD b(matrixSize, matrixSize);

	//Loading values of data vector into matrices A and B
	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = 0; j < matrixSize; j++)
		{

			a[i][j] = data.at(counter);
			b[i][j] = data.at(counter);
			counter++;
		}
	}

	//Printing the values of matrix A to console
	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = 0; j < matrixSize; j++)
		{
			cout << a[i][j] << ",       ";
		}
		cout << "\n\n";
	}

	//Inversion of matrix A and printing new values to console
	a.Invert();
	cout << "Matrix A after inversion" << endl;

	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = 0; j < matrixSize; j++)
		{
			std::cout << setprecision(2) << a[i][j] << ",       ";
		}
		std::cout << "\n\n";
	}

	//Multiplying the inverted matrix by the original (now stored as matrix B) in order to get the identity matrix
	a *= b;
	cout << "Inverted matrix A multiplied by its original (identity matrix)" << endl;

	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = 0; j < matrixSize; j++)
		{
			std::cout << setprecision(2) << a[i][j] << ",   ";
		}
		std::cout << "\n\n";
	}
}
