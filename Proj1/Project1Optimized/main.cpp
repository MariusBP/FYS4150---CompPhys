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
    vec dt(n+1);
    vec ft(n+1);
    vec v(n+1);
    v(0) = 0;
    v[n] = 0;
    vec x;
    x = linspace(0,1,n+1);
    vec f;
    double h = (x[n]-x[0])/n;
    f = h*h*100*exp(-10*x);
    vec u = 1 - (1 - exp(-10))*x - exp(-10*x);

    dt[0] = 2;
    ft[0] = f[0];

    //forwards substitution
    for(double i = 1; i<n+1; i++){
        dt[i] = (i+1)/i;
        ft[i] = f[i] + ft[i-1]/dt[i-1];
    }

    v[n-1] = ft[n-1]/dt[n-1];

    //backwards substitution
    for(int i=n; i>0; i--){
        v[i] = (ft[i] + v[i+1])/dt[i-1];
    }
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
