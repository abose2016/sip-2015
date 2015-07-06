

void CatmullRomSpline(P0, P1, P2, P3, nPoints=100)
{
  vector <double> P0, P1, P2, P3; 
//Calculate t0 to t4
  double alpha = 0.5;
  void tj(ti, Pi, Pj)
			{
    	xi, yi = Pi
    	xj, yj = Pj
    	return (( pow(( pow((xj-xi),2) + pow(pow((yj-yi),2), .5)), alpha) + ti)
			}
  t0 = 0
  t1 = tj(t0, P0, P1)
  t2 = tj(t1, P1, P2)
  t3 = tj(t2, P2, P3)

  //Only calculate points between P1 and P2
  t = numpy.linspace(t1,t2,nPoints)

  // Reshape so that we can multiply by the points P0 to P3
  // and get a point for each value of t.
  t = t.reshape(len(t),1)

  A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1
  A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2
  A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3

  B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2
  B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3

  C  = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2
  return C

	void CatmullRomChain(P)
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

main(){
//Define a set of points for curve to go through
Points = [[0,6],[2,2],[3,.5],[4,8],[5,3],[6,1],[7,3.5]]

//Calculate the Catmull-Rom splines through the points
c = CatmullRomChain(Points)

// Convert the Catmull Rom curve points into x and y arrays and plot
x,y = zip(*c)
plt.plot(x,y)

//Plot the control points
px, py = zip(*Points)
plt.plot(px,py,'or')

plt.show()
}

