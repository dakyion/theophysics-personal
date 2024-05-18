#include<iostream>
#include<fstream>
// sudo apt install gnuplot
int main(){
    const int size = 100; // grid size
    const double d = size- 1; // (delta)^-1
    const int num_iter = 100;
    double pot[size][size];

    std::ofstream outfile("ex2.dat");
    for (int i = 0; i < size; i++)
    {
        pot[i][0] = 0.0; // 
        for (int j = 1; j < size; j++)
        {
            pot[i][j] = 0; // rest of grid is zero as an initial Guess!
        }
    }
    for (int iter = 0; iter < num_iter; iter++)
    {
        for (int i = 1; i < size-1; i++)
        {
            for (int j = 1; j < size-1; j++)
            {
                pot[i][j] = 0.25*(pot[i+1][j]+pot[i-1][j]+pot[i][j+1]+pot[i][j-1]-1/(d*d)*i/d*j/d);
            }
            int j =size -1;
            pot[i][j] = 0.25*(pot[i+1][j]+pot[i-1][j]+2*i*i*1/d*1/d*1/d+2*pot[i][j-1]-1/(d*d)*i/d*j/d);
        }
    }
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            outfile << pot[i][j] << "\n";
        }
        outfile << "\n";
    }
    std::cout << "data stored \n";
    outfile.close();
    return 0;
}