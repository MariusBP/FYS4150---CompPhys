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
    int start = clock();
    int n;
    cout << "Give me n valuables! *Points gun*" << endl;
    cin >> n;
    //initialisation
    vec d = randu(n+1);
    vec e = randu(n+1);
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
    vec sol = 1 - (1 - exp(-10))*x - exp(-10*x);

    dt[0] = d[0];
    ft[0] = f[0];

    //forwards substitution
    for(int i = 1; i<n+1; i++){
        dt[i] = d[i] - (e[i-1]*e[i-1])/dt[i-1];
        ft[i] = f[i] - e[i-1]*ft[i-1]/dt[i-1];
    }

    v[n-1] = ft[n-1]/dt[n-1];

    //backwards solution
    for(int i=n-1; i>0; i--){
        v[i] = (ft[i] - e[i]*v[i+1])/dt[i];
    }
    ofstream myfile;
    string name = "values" + to_string(n) + ".txt";
    myfile.open (name);
    for(int i = 0; i<n+1; i++){
        myfile << v[i] << "\t" <<sol[i] <<"\n";
    }
    myfile.close();
    cout << "Program run time: "<<(clock()-start)/float(CLOCKS_PER_SEC)<< " seconds" << endl;
}
