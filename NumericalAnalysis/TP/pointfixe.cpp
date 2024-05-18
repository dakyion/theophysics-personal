//  WRITTEN BY ABDENNOUR HARCHE
// 17/11/2022

#include <iostream>
#include <cmath>

double g1(double x1,double x2,double x3){
    return 1/12.*(pow(x2,2)+2*exp(-x3));
}
double g2(double x1,double x2,double x3){
    return 1/6.*(1-x1+sin(x3));
}
double g3(double x1,double x2,double x3){
    return 1/6.*(pow(x1,2)+pow(x2,2)+pow(x3,2));
}
double norm3(double* v){
    return sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
}
int main(){
    double err = 1;
    double tol = 1e-9;
    int niter=0;
    int nmax=20;
    double x[3]={0.,0.,0.};
    
    while (err>tol && niter<nmax)
    {
        double y[3];
        for (int i = 0; i < 3; ++i){y[i]=x[i];}
        x[0] = g1(y[0],y[1],y[2]);
        x[1] = g2(y[0],y[1],y[2]);
        x[2] = g3(y[0],y[1],y[2]);
        ++niter;
        double v[3] ;
        for (int i(0); i < 3; ++i){v[i] = x[i] - y[i];}
        err = norm3(v);
        printf("%2i \t %.10f \t %.10f \t %.10f \t %.10g \n",niter,x[0],x[1],x[2],err);
    }
    std::cout << "END\n";
}