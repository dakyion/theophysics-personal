//  WRITTEN BY ABDENNOUR HARCHE
// 16/12/2022

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
vector<double> gjsolve(vector<vector<double>>& A, vector<double> const& B, int dimension){
    
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
double minV(vector<double> v,int dim){
    double t= abs(v[0]);
    for (int i = 1; i < dim; i++)
    {
        if(v[i] !=0 && t > abs(v[i])) t = v[i];
    }
    return t;
}
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
// Fill a matrix
vector<vector<double>> fillM(int p,int q){
    vector<vector<double>> m(q);
    for(int i(0);i<p;i++) m[i].reserve(3);
    std::cout << "Enter The matrix: \n";
    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < p; j++) std::cin >> m[i][j];
    }
    return m;
}
// Fill a vector
vector<double> fillV(int dim){
    vector<double> v(dim);
    std::cout << "Enter Vector: \n";
    for (int i = 0; i < dim; i++) std::cin >> v[i];
    return v;
}

int main(){
    
    // vector<vector<double>> A={{5,-4,6},{3,1,3},{2,4,1}};
    vector<vector<double>> A=fillM(3,3);
    
    vector<double> y0 = {0,1,1};
    vector<double> y1 = multmv(A,y0,3,3);
    vector<double> y2 = multmv(A,y1,3,3);
    vector<double> y3 = multmv(A,y2,3,3);
    y3 = multsv(-1,y3,3);
    
    vector<vector<double>> Y = {y2,y1,y0};
    Y = trans(Y,3,3);
    displayV(y3,3);
    displayM(Y,3,3);

    vector<double> x = gjsolve(Y,y3,3);
    std::cout << "The coefficients of the Caracteristic Polynomial of the Matrix \n";
    displayV(x,3);
    std::cout << "Use them to determine the eigenvalues and give them to the waiting prompt: \n";
    vector<double> l = fillV(3);
    // double s[3] = {0 , 0, 0};
    // double p[3] = {1,1,1};
    // vector<vector<double>> U(3);
    // for(int i = 0;i<3;i++) U[i].resize(3);
    // for (int i = 0; i < 3; i++)
    // {
    //     for (int j = 0; j < 3; j++) if(i != j) {s[i]+=l[j]; p[i]*=l[j];}
    //     vector<double> U[i] = addv(addv(y2, multsv(-s[i], y1, 3), 3), multsv(p[i], y0, 3), 3);
    //     displayV(multsv(1 / minV(U[i], 3), U[i], 3), 3);
    // }
    std::cout << "Eigenvectors are : \n";
    double s1 = l[1]+l[2];
    double p1 = l[1]*l[2];
    vector<double> u1 = addv(addv(y2,multsv(-s1,y1,3),3),multsv(p1,y0,3),3);
    displayV(multsv(1/minV(u1,3),u1,3),3);
    double s2 = l[0]+l[2];
    double p2 = l[0]*l[2];
    vector<double> u2 = addv(addv(y2,multsv(-s2,y1,3),3),multsv(p2,y0,3),3);
    displayV(multsv(1/minV(u2,3),u2,3),3);
    double s3 = l[0]+l[1];
    double p3 = l[0]*l[1];
    vector<double> u3 = addv(addv(y2,multsv(-s3,y1,3),3),multsv(p3,y0,3),3);
    displayV(multsv(1/minV(u3,3),u3,3),3);
    return 0;
}