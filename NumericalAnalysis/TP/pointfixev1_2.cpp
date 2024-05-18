//  WRITTEN BY ABDENNOUR HARCHE
// 17/11/2022

#include <iostream>
#include <cmath>
#include <vector>
using std::vector;

// Define the fixed point functions vector G(X)
vector<double> G(vector<double> const& x){
    return {
        sin(x[0]+x[1]),
        cos(x[0]+x[1])
    };
}

// Quadratic Norm
double norm3(vector<double> v){
    return sqrt(pow(v[0],2)+pow(v[1],2));
}

int main(){
    double err = 1;
    // Accepted error
    double tol = 1e-6;
    int niter=0;
    int nmax=50;
    vector<double> x ={0.,0.};
    // Store a copy of the current variables' value.
    vector<double> y = x;
    while (err>tol && niter<nmax)
    {
        // Apply the method
        x = G(x);
        ++niter;
        vector<double> v(3); 
        for(int i(0);i<3;++i) {v[i]=x[i]-y[i];}
        // calculate the error.
        err = norm3(v);
        printf("%2i \t %.7f \t %.7f  \t %.8g \n",niter,x[0],x[1],err);
        y=x;
    }
    std::cout << "END\n";
}