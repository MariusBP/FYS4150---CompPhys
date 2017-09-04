//Project 1
#include <iostream>
#include <armadillo>
#include <vector>
#include <ctime>
using namespace std;
using namespace arma;

int main()
  {
    clock_t start;
    int n;
    cout << "Give me your valuables! *Points gun*" << endl;
    cin >> n;
    //initialisation
    vec d(n+1);
    vec e(n);
    vec dt(n+1);
    vec ft(n+1);
    vec u(n+1);
    u(0) = 0;
    u[n] = 0;
    vec x;
    x = linspace(0,1,n+1);
    vec f;
    double h = (x[n]-x[0])/n;
    f = h*h*100*exp(-10*x);

    d.fill(2.0);
    e.fill(-1.0);
    dt[0] = d[0];
    ft[0] = f[0];

    for(int i = 1; i<n+1; i++){
        dt[i] = d[i] - (e[i-1]*e[i-1])/dt[i-1];
        ft[i] = f[i] - e[i-1]*ft[i-1]/dt[i-1];
    }

    u[n-1] = ft[n-1]/dt[n-1];

    for(int i=n-1; i>0; i--){
        u[i] = (ft[i] - e[i]*u[i+1])/dt[i];
    }
    cout << (clock()-start)/float(CLOCKS_PER_SEC) << endl;
    for(int i = 0; i<n+1; i++)
    cout << u[i] << "\t" << 1 - (1 - exp(-10))*x[i] - exp(-10*x[i])<< endl;
}
