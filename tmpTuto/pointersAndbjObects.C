#include <iostream>
#include <vector>

using namespace std;
// Simple examples of various ways of defining and using values,
// pointers and references.
int f1(int n)
// Ordinary call by value semantics
{
    n = n+1;
    return 2*n;
}


int f2(int *np) {
    // Value is passed in via a pointer.  The value pointed to may be
    // changed in the calling function.  Note the parentheses around
    // *np - without these, it would have a different effect.
    (*np)++;
    return 2*(*np);
}

// THis one uses a reference parameter.  Within the code, (and look at
// the caller), it looks like pass by value, but the behavior is as if
// you passed a pointer.  The term for this sort of thing is
// "syntactic sugar".

int f3(int &n) {
    n = n+1;
    return 2*n;
}

// Pass by reference - thus is usual with vectors
void load_fib(vector<double> &dv, int n) {
    dv.push_back(0.0);
    dv.push_back(1.0);
    for (int i =0; i<n; i++) {
        dv.push_back(*(dv.end()-1) + *(dv.end()-2));
    }
}

void print_vector(vector<double> &myvec) {
    cout << "[ ";
    for (vector<double>::iterator it = myvec.begin(); it != myvec.end(); it++) {
        cout << *it << ", ";   // NO endl
    }
    cout << "]" << endl;

}


int main() {

    cout << "Hello world!  Lets call some functions!" << endl;


    int i=7;
    int j = f1(i);
    cout << "output from f1 with argument 7 is " << j <<
      " and argument variable is now " << i << endl;

    // Note that we have to pass the address of i, rather than i itself.
    i = 7;
    j = f2(&i);
    cout << "output from f2 with argument 7 is " << j <<
      " and argument variable is now " << i << endl;

    i=7;
    j = f3(i);  // Note - NOT &i in the caller.
    cout << "output from f3 with argument 7 is " << j <<
      " and argument variable is now " << i << endl;

    // Same idea, but now with vectors

    cout << "Fibonacci numbers with a vector" << endl;
    vector<double> av;
    load_fib(av, 16);

    print_vector(av);

    return 0;
}
