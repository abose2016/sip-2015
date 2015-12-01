#include <vector>

using namespace std;

double tj(double ti, double P0[], double P1[], double alpha = 0.5) 
{
		double xi = P0[0];
		double yi = P0[1];
		double xj = P1[0];
		double yj = P1[1];

    	return (( pow(( pow((xj-xi),2) + pow(pow((yj-yi),2), .5)), alpha) + ti)
	}


void CatmullRomSpline(double P0[], double P1[], double P2[], double P3[], int nPoints=100)
{
	double P0[], P1[], P2[], P3[]; 

	//Calculate t0 to t4
	double alpha = 0.5;
  
	double t0 = 0
	double t1 = tj(t0, P0, P1)
 	double t2 = tj(t1, P1, P2)
	double t3 = tj(t2, P2, P3)

	double t[];
	t[0] = t1;
  //Only calculate points between P1 and P2
  for(int i = 0; i < nPoints; i++)
	{
			t[i+1] = (t2 - t1)/nPoints + t[i];
	}
  // Reshape so that we can multiply by the points P0 to P3
  // and get a point for each value of t.
	
	vector<double> reshapedPoints;
	for(int i = 0; i < t.size(); i++)
	{
		reshapedPoints.push_back(t[i]);
	}

  A1 = (t1-reshapedPoints)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1
  A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2
  A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3

  B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2
  B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3

  C  = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2
  return C

	void CatmullRomChain(TMatrixD *P)
	{
		int sz= P.size();

		// The curve C will contain an array of (x,y) points.
		C = []
		for i in range(sz-3)
		{
		  c = CatmullRomSpline(P[i], P[i+1], P[i+2], P[i+3])
		  C.extend(c)
		}
		return C
	}
}

void main(){
	int matrixSize = 7;
	//Define a set of points for curve to go through
	TMatrixD points(matrixSize, matrixSize);
	points[0]= [0,6];
	points[1]=[2,2];
	points[2]=[3,.5];
	points[3]=[4,8];
	points[4]=[5,3];
	points[5]=[6,1];
	points[6]=[7,3.5];

	//Calculate the Catmull-Rom splines through the points
	c = CatmullRomChain(points)

	// Convert the Catmull Rom curve points into x and y arrays and plot
	c[0] = x[]
	x,y = zip(*c)
	plt.plot(x,y)

	//Plot the control points
	px, py = zip(*Points)
	plt.plot(px,py,'or')

	plt.show()
}

