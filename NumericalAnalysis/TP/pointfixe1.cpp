//  WRITTEN BY ABDENNOUR HARCHE
// 17/11/2022

#include <iostream>
#include <cmath>
#include <vector>
using std::vector;

vector<double> G(vector<double> const& x){
    return {
        1/12.*(pow(x[1],2)+2*exp(-x[2])),
        1/6.*(1-x[0]+sin(x[2])),
        1/6.*(pow(x[0],2)+pow(x[1],2)+pow(x[2],2))
    };
}
double norm3(vector<double> v){
    return sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
}

int main(){
    double err = 1;
    double tol = 1e-9;
    int niter=0;
    int nmax=20;
    vector<double> x ={0.,0.,0.};
    vector<double> y = x;
    while (err>tol && niter<nmax)
    {
        x = G(x);
        ++niter;
        vector<double> v(3); 
        for(int i(0);i<3;++i) {v[i]=x[i]-y[i];}
        err = norm3(v);
        printf("%2i \t %.10f \t %.10f \t %.10f \t %.10g \n",niter,x[0],x[1],x[2],err);
        y=x;
    }
    std::cout << "END\n";
}