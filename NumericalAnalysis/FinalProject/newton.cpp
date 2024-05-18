//  WRITTEN BY ABDENNOUR HARCHE
// 24/12/2022

#include <iostream>
#include <cmath>
#include <vector>
using std::vector ;
// Add two matrices
vector<vector<double>> addm(vector<vector<double>> const& A, vector<vector<double>> const& B, int const& m,int const& n){
    vector<vector<double>> x(m);
    for (int i = 0; i < m; ++i) x[i].resize(n);

    for(int i(0);i<m;++i){
        for (int j = 0; j < n; ++j) x[i][j]= A[i][j] + B[i][j];
    }

    return x;
}
// Add two vectors
vector<double> addv(vector<double> const& v, vector<double> const& u,int const& d){
    vector<double> x(d);
    for (int j = 0; j < d; ++j) x[j]= u[j] + v[j];
    return x;
}
// Multiply a vector by a scalar
vector<double> multsv(double s, vector<double> const& v, int dim){
    vector<double> V = v;
    for (int i = 0; i < dim; i++) V[i]*=s;
    return V;
}
// Multiply a matrix by a vector
vector<double> multmv(vector<vector<double>> const& A, vector<double> const& v, int const& m,int const& n){
    vector<double> x(n);
    
    for (int i = 0; i < m; ++i)
    {
        double s = 0;
        for (int j = 0; j < n; ++j) s+=A[i][j]*v[j] ;
        x[i]= s;
    }
    return x;
}
// Multiply a matrix by a matrix
vector<vector<double>> multmm(vector<vector<double>> const& A, vector<vector<double>> const& B, int m,int n,int p){
    vector<vector<double>> x(m);
    for (int i = 0; i < m; ++i) x[i].resize(p);
    for (int i = 0; i < m; ++i)
    {
        for (int r = 0; r < p; ++r)
        {
            double s = 0;
            for (int j = 0; j < n; ++j) s += A[i][j] * B[j][r];
            x[i][r] = s;
        }
    }
    return x;
}
// Transpose a matrix
vector<vector<double>> trans(vector<vector<double>> const& A, int p,int q){
    vector<vector<double>> T = A;
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < q; j++)
        {
            if (i>=j) continue;
            
            double _ = T[i][j];
            T[i][j] = T[j][i];
            T[j][i] = _;
        }
    }
    return T;
}
// Generate Gauss-Jordan elimination matrices
vector<vector<double>> Pn(vector<vector<double>> const& A,int const& dimension,int const& q){
    vector<vector<double>> P(dimension);
    for (int i = 0; i < dimension; i++) P[i].resize(dimension);
    
    for (int i = 0; i < dimension; ++i)
    {
        for (int j = 0; j < dimension; ++j)
        {
            if (i == j && i != q) P[i][i] = 1.;
            else if (i == j && i == q) P[i][i] = 1. / A[i][i];
            else if (j == q) P[i][j] = -A[i][j] / A[j][j];
            else P[i][j] = 0.;
        }
    }
    return P;
}
// Solve a linear system using Gauss-Jordan Method
vector<double> gjsolve(vector<vector<double>> A, vector<double> const& B, int dimension){
    
    vector<vector<double>> P(dimension);
    for (int i = 0; i < dimension; i++) P[i].resize(3);
    vector<double> Pb = B;
    for (int q = 0; q < dimension; ++q)
    {
        P = Pn(A,dimension,q);
        Pb = multmv(P,Pb,dimension,dimension);
        A =  multmm(P,A,dimension,dimension,dimension);
    }

    return Pb;
}

// Get the minimum of a vector

// Display a vector into the console
void displayV(vector<double> const& V,int dim){
    for (int j = 0; j < dim; j++) std::cout << V[j] << "\t";
    std::cout << "\n\n";
}
// Display a matrix into the console
void displayM(vector<vector<double>> const& M, int n, int p){
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < p; j++) std::cout << M[i][j] << "\t";
        std::cout << "\n";
    }
    std::cout << "\n";
}
// Fill a vector
vector<double> fillV(int dim){
    vector<double> v(dim);
    std::cout << "Enter Vector: \n";
    for (int i = 0; i < dim; i++) std::cin >> v[i];
    return v;
}

vector<double> F(vector<double> const& x){
    return {
        1/12.*(pow(x[1],2)+2*exp(-x[2])),
        1/6.*(1-x[0]+sin(x[2])),
        1/6.*(pow(x[0],2)+pow(x[1],2)+pow(x[2],2))
    };
}
vector<vector<double>> J(vector<double> const& x){
    return {
        {1/12.*(pow(x[1],2)+2*exp(-x[2]))},
        {1/6.*(1-x[0]+sin(x[2]))},
        {1/6.*(pow(x[0],2)+pow(x[1],2)+pow(x[2],2))}
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
    vector<double> x = fillV(3);
    vector<double> y = x;
    while (err>tol && niter<nmax)
    {
        vector<double> H = gjsolve(J(x),multsv(-1,F(x),3),3);
        x= addv(y,H,3);
        ++niter;
        vector<double> v(3); 
        for(int i(0);i<3;++i) {v[i]=x[i]-y[i];}
        err = norm3(v);
        printf("%2i \t %.10f \t %.10f \t %.10f \t %.10g \n",niter,x[0],x[1],x[2],err);
        y=x;
    }
    std::cout << "END\n";
}