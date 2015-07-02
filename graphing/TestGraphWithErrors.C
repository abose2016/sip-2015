#include "TCanvas.h"// This lets me draw a canvas
#include "TGraphErrors.h"// this is for the class TGraphErrors
#include <vector>//This lets me use vectors
#include <iostream>//I'm not quite sure why I need this. look up?

using namespace std;// This means that I can use all of the functions from the standard library without having to say std:: in front of everything. I know this has some disadvantage because Dr. Biteau doesn't do it. look up?

void TestGraphWithErrors() { // naming a function. Remember to name it the same name as the file so it runs automatically. There is a way to load a file and then call specific functions but I haven't figured it out yet. look up?

   TCanvas *c1 = new TCanvas("c1","My Awesome Test Graph!!",200,10,700,500);// Here I'm making a new instance of the TCanvas method. 

   c1->SetFillColor(70);//Sets color of background.70 is blue-green, 80 is green, 90 is yellow, 100 is red.
  c1->SetGrid();// no idea. I tried commenting it out and nothing discernable changed.
  c1->GetFrame()->SetFillColor(100);// same with this one. 
  c1->GetFrame()->SetBorderSize(1000);// and this one


   const int n = 10;// This is the size of my arrays and vectors
   double x[n]  = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};// these are the x coordinates of each point.
   double y[n]  = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
   double ex[n] = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};// x coordinates of the errors.
   double ey[n] = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};



vector<double> xVector; // defines an empty vector
vector<double> yVector;
vector<double> xErrorVector;
vector<double> yErrorVector;

for(int i =0; i<n; i=i+1){
xVector.push_back(x[i]);
yVector.push_back(y[i]);
xErrorVector.push_back(ex[i]);
yErrorVector.push_back(ey[i]);
}// this for loop fills my vectors using the push_back() method to add values to the end of the vector one after the other until we've filled it up to n places. It has a dot between the vector and the function because vectors are objects rather than pointers. If it were an array (pointer to first thing in the array), it would use an -> instead.




   TGraphErrors *gr = new TGraphErrors(n,&xVector[0],&yVector[0],&xErrorVector[0],&yErrorVector[0]);// Here I call a new instace of the TGraphErrors thing. I'm being a little bit sneaky in my coding because the TGraphErrors doesn't have a version that takes the kind of vectors I've defined as arguments. But it does have a version that takes arrays as arguments so I used &s in front of my vectors so that I could treat them like arrays. If you have an & before an object, it gives you the address of that object. If you use * in front of a pointer, it gives you the thing in that address. Then you can treat it as if it were the actual thing. 
 
   gr->SetTitle("Test Graph With Error Bars");// Title
   gr->SetMarkerColor(2);// This sets the marker color
   gr->SetMarkerStyle(21);// This actually draws the markers, so if you don't have this, then you won't have markers.
   gr->Draw("ALP");// Not sure what this does, but if you comment it out, it won't draw you a graph. 


   c1->Update();// I think this updates the canvas you created at the beggining to include all of the stuff we just called. 

}
 
