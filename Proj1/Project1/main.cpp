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
    int start = clock();//initialising clock

    vec d(n+1);//creating the d vector and filling it with 2s
    d.fill(2);    
    vec e(n+1);//creating the e-vector and filling it with -1s
    e.fill(-1);
    vec dt(n+1); //d-tilde
    vec ft(n+1);//f-tilde
    vec v(n+1); //numerical solution with initial values
    v(0) = 0;
    v[n] = 0;
    vec x; //x from 0 to 1
    x = linspace(0,1,n+1);
    vec f;
    double h = (x[n]-x[0])/n; //step length
    f = h*h*100*exp(-10*x); //in task b-tilde
    vec u = 1 - (1 - exp(-10))*x - exp(-10*x); //exact solution

    dt[0] = d[0];//initialising d-tilde and f-tilde
    ft[0] = f[0];

    //forwards substitution
    for(int i = 1; i<n+1; i++){
        dt[i] = d[i] - (e[i-1]*e[i-1])/dt[i-1];
        ft[i] = f[i] - e[i-1]*ft[i-1]/dt[i-1];
    }

    v[n-1] = ft[n-1]/dt[n-1];//initial value for numerical solution

    //backwards substitution
    for(int i=n; i>0; i--){
        v[i] = (ft[i] - e[i]*v[i+1])/dt[i];
    }
    /*
    ofstream myfile;//writing to file
    string name = "values" + to_string(n) + ".txt";
    myfile.open (name);
    for(int i = 0; i<n+1; i++){
        myfile << x[i] <<"\t" << v[i] << "\t" <<u[i] <<"\n";
    }
    myfile.close();
    */
    vec err = log10(abs((v-u)/u));//error
    err(0) = 0;//setting the endpoints to 0 since because of initial values they are inf
    err(n) = 0;
    cout << "The error is: "<< err.max()<<endl;//printing the max error.
    cout << "Program run time: "<<(clock()-start)/float(CLOCKS_PER_SEC)<< " seconds" << endl;//printing time taken

}
