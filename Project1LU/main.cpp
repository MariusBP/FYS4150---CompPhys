//Project 1
#include <iostream>
#include <armadillo>
#include <vector>
#include <ctime>
#include <cstdio>

using namespace std;
using namespace arma;

int main()
  {
    int n;
    cout << "Give me n valuables! *Points gun*" << endl;
    cin >> n;
    //initialisation
    int start = clock();
    mat A(n+1,n+1);
    vec dt(n+1);
    vec ft(n+1);
    vec x;
    x = linspace(0,1,n+1);
    vec f;
    double h = (x[n]-x[0])/n;
    f = h*h*100*exp(-10*x);
    vec u = 1 - (1 - exp(-10))*x - exp(-10*x);


    A(0,0) = 2;
    A(0,1) = -1;
    A(n,n) = 2;
    A(n,n-1) = -1;
    dt[0] = A(0,0);
    ft[0] = f[0];

    //Setting up the rest of A
    for(int i=1; i<n; i++){
        A(i,i) = 2;
        A(i,i-1) = -1;
        A(i,i+1) = -1;
    }

    vec v = solve(A,f);

    ofstream myfile;
    string name = "values" + to_string(n) + ".txt";
    myfile.open (name);
    for(int i = 0; i<n+1; i++){
        myfile << x[i] <<"\t" << v[i] << "\t" <<u[i] <<"\n";
    }
    myfile.close();

    vec err = log10(abs((v-u)/u));
    err(0) = 0;
    err(n) = 0;
    cout << "The error is: "<< err.max()<<endl;
    cout << "Program run time: "<<(clock()-start)/float(CLOCKS_PER_SEC)<< " seconds" << endl;

}
