//  WRITTEN BY ABDENNOUR HARCHE
// 17/11/2023

#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include<fstream>
#include<string>

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

double dotVV(vector<double> const& a,vector<double> const& b, int dim){
    double sum=0.;
    for (int i = 0; i < dim; i++)
    {
        sum += a[i]*b[i];
    }
    return sum;
}

class Vortex{
    public:
        double charge;
        vector<double> position = {0,0};
        double angle;
        double velocity = 0.;
        Vortex(){}
        Vortex(double c,double x,double y){
            charge = c;
            position = {x,y};
            angle = atan2(position[1],position[0]);
        }
};

vector<vector<double>> updateField(vector<Vortex> const& v, int const& n, double const& k){
    vector<vector<double>> ef (2*n);
    for (int i = 0; i < 2*n; i++)
    {
        ef[i] = {0.,0.};
    }
    for (int i = 0; i < 2*n-1; i++)
    {
        for (int j = i+1; j < 2*n; j++)
        {
            vector<double> r = addv(v[i].position,multsv(-1,v[j].position,2),2);
            double nrm = 1./dotVV(r,r,2);
            vector<double> def1 = multsv(k*v[j].charge*nrm,r,2);
            vector<double> def2 = multsv(-k*v[i].charge*nrm,r,2);
            ef[i] = addv(ef[i],def1,2);
            ef[j] = addv(ef[j],def2,2);
        }   
    }
    return ef;
}

void updatePosition(vector<Vortex>& v,int const& n, vector<vector<double>> const& ef, double const& dt){
    
    for (int i = 0; i < 2*n; i++)
    {  
        double dx = ef[i][1]*dt;
        double dy = -ef[i][0]*dt;
        vector<double> dr = {dx,dy};
        v[i].position = addv(v[i].position,dr,2);
    }
}

int main(){
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 10.0);

    double q = 1.602e-19;
    double e0 = 8.854e-12;
    double dt = 1e-12;
    int n;
    std::cout << "Enter n \t" ;
    std::cin >> n;
    double scale = 1e-8;
    double ldensity = q * 1000;
    double k = 1./(2*e0*M_PI);
    // vector<vector<double>> electricField (2*n);
    // for (int i = 0; i < n; i++)
    // {
    //     electricField[i].reserve(2);
    //     electricField[i][0] = 0.;
    //     electricField[i][1] = 0.;
    // }
    
    vector<Vortex> vortices (2*n);
    for (int i = 0; i < 2*n; i++)
    {
        if (i<n)
        {
            double x = dis(gen)*scale;
            double y = dis(gen)*scale;
            vortices[i] = Vortex(ldensity,x,y);
        }else{
            double x = dis(gen)*scale;
            double y = dis(gen)*scale;
            vortices[i] = Vortex(-ldensity,x,y);
        }
    }
    
    for (int step = 0; step < 1000; step++)
    {
        vector<vector<double>> electricField = updateField(vortices,n,k);
        updatePosition(vortices,n,electricField,dt);
        std::string filename = "data/step" + std::to_string(step) + ".dat";
        std::ofstream outfile(filename);
        for(Vortex vort : vortices){
            outfile << vort.position[0] << "\t" << vort.position[1] << "\n";
        }
        // for (int i = 0; i < 2*n; i++)
        // {
        //     // vector<double> p = vortices[i].position;
        //     outfile << vortices[i].position[0] << "\t" << vortices[i].position[1] << "\n";
        // }
        
        outfile.close();
    }
    
    return 0;
}
