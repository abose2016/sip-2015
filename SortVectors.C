//Look for the std library sort in <algorithm>


using namespace std;
#include <vector>
#include <algorithm>


void SortVectors(){

    std::vector<double> vx;
    vx.push_back(3.5);
    vx.push_back(3);
    vx.push_back(3.141592653589793);

    cout << "original "  << vx[0] << " " << vx[1] << " " << vx[2] << " " << endl;
    std::sort(vx.begin(), vx.end());
    cout << "sorted "  << vx[0] << " " << vx[1] << " " << vx[2] << " " << endl;
};
