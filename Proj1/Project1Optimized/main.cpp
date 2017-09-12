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
    int n;//taking n as input from command line
    cout << "Give me n valuables! *Points gun*" << endl;
    cin >> n;
    //initialisation
    int start = clock(); //Tracking time
    vec dt(n+1); //d-tilde
    vec ft(n+1); //f-tilde
    vec v(n+1); //our numerical approximation to u
    v(0) = 0;//initial values
    v[n] = 0;
    vec x; //x in (0,1)
    x = linspace(0,1,n+1);
    vec f;
    double h = (x[n]-x[0])/n; //our step length
    f = h*h*100*exp(-10*x);
    vec u = 1 - (1 - exp(-10))*x - exp(-10*x);//the analytical solution

    dt[0] = 2;//initial values for d-tilde and f-tilde
    ft[0] = f[0];

    //forwards substitution optimized code
    for(double i = 1; i<n+1; i++){
        dt[i] = (i+1)/i;
        ft[i] = f[i] + ft[i-1]/dt[i-1];
    }

    v[n-1] = ft[n-1]/dt[n-1]; //initial value for v

    //backwards substitution
    for(int i=n; i>0; i--){
        v[i] = (ft[i] + v[i+1])/dt[i];
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
    vec err = log10(abs((v-u)/u)); //Calculating error and ignoring the infinite terms at the end-points
    err(0) = 0;
    err(n) = 0;
    cout << "The error is: "<< err.max()<<endl; //printing the results
    cout << "Program run time: "<<(clock()-start)/float(CLOCKS_PER_SEC)<< " seconds" << endl;

}
