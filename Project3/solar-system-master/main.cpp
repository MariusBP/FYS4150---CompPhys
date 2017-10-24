#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "solarsystem.h"
#include "verlet.h"
using namespace std;

int main(int numArguments, char **arguments)
{
    int numTimesteps = 100000;
    if(numArguments >= 2) numTimesteps = atoi(arguments[1]);

    SolarSystem solarSystem;
    // We create new bodies like this. Note that the createCelestialBody function returns a reference to the newly created body
    // This can then be used to modify properties or print properties of the body if desired
    // Use with: solarSystem.createCelestialBody( position, velocity, mass );
    // We don't need to store the reference, but just call the function without a left hand side

    double M_sun = 2e30;//kg
    /*
    //solarSystem.createCelestialBody(vec3(-5.2e-3,0,0), vec3(0,-2.79e-3,0), 1.0 , "the Sun"); //This is for 3-body Sun-Earth-Jupiter system
    solarSystem.createCelestialBody(vec3(-1.02e-2,0,0), vec3(0,-3.51e-3,0), 1.0 , "the Sun");
    solarSystem.createCelestialBody(vec3(0.39, 0, 0), vec3(0, 2*M_PI*0.39/0.24, 0), 3.3e23/M_sun, "Mercury");
    solarSystem.createCelestialBody(vec3(0.72, 0, 0), vec3(0, 2*M_PI*0.72/0.615, 0), 4.9e24/M_sun, "Venus");
    solarSystem.createCelestialBody(vec3(1, 0, 0), vec3(0, 2*M_PI, 0), 3e-6 , "Earth");
    solarSystem.createCelestialBody(vec3(1.52, 0, 0), vec3(0, 2*M_PI*1.52/1.88, 0), 6.6e23/M_sun, "Mars");
    solarSystem.createCelestialBody(vec3(5.2, 0, 0), vec3(0, 2*M_PI*5.2/11.8, 0), 1e-3, "Jupiter");
    solarSystem.createCelestialBody(vec3(9.54, 0, 0), vec3(0, 2*M_PI*9.54/29.44, 0), 5.5e26/M_sun, "Saturn");
    solarSystem.createCelestialBody(vec3(19.19, 0, 0), vec3(0, 2*M_PI*19.19/84.01, 0), 8.8e25/M_sun, "Uranus");
    solarSystem.createCelestialBody(vec3(30.06, 0, 0), vec3(0, 2*M_PI*30.06/164.79, 0), 1.03e26/M_sun, "Neptune");
    solarSystem.createCelestialBody(vec3(39.53, 0, 0), vec3(0, 2*M_PI*39.53/247.92, 0), 1.31e22/M_sun, "not-a-planet");
    */
    solarSystem.createCelestialBody(vec3(2.187003065211543E-03, 5.768166559108312E-03,-1.294147734354897E-04), vec3(-5.280315745929900E-06, 5.460827268770582E-06, 1.244590551885596E-07)*365.242199, 1.0 , "the Sun");
    solarSystem.createCelestialBody(vec3(-2.139370590565288E-01,-4.028814669327753E-01,-1.369419923866817E-02), vec3(1.919605110412534E-02,-1.178939724774313E-02,-2.725194971976164E-03)*365.242199, 3.3e23/M_sun, "Mercury");
    solarSystem.createCelestialBody(vec3(-6.915411411024813E-01, 1.907707656342085E-01, 4.244080591200550E-02), vec3(-5.318149817917223E-03,-1.962946252497859E-02,3.737564099986648E-05)*365.242199, 4.9e24/M_sun, "Venus");
    solarSystem.createCelestialBody(vec3(8.589987108796383E-01, 5.110680605545673E-01,-1.568623415833688E-04), vec3(-9.031139666444341E-03, 1.476253462344845E-02,-1.943130198126555E-07)*365.242199, 3e-6 , "Earth");
    solarSystem.createCelestialBody(vec3(-1.590488403238053E+00, 4.879788693373922E-01, 4.906264799695509E-02), vec3(-3.536379998688235E-03,-1.219225340287500E-02,-1.688251041034333E-04)*365.242199, 6.6e23/M_sun, "Mars");
    solarSystem.createCelestialBody(vec3(-4.556745348155565E+00, -2.963008457339381E+00, 1.142108603087190E-01), vec3(4.025290127025506E-03, -5.968165948463064E-03, -6.526566072536134E-05)*365.242199, 1e-3, "Jupiter");
    solarSystem.createCelestialBody(vec3(-3.158467085324504E-01,-1.005065028034512E+01, 1.873222298678130E-01), vec3(5.270276974251006E-03,-1.927338311574419E-04,-2.066633331055191E-04)*365.242199, 5.5e26/M_sun, "Saturn");
    solarSystem.createCelestialBody(vec3(1.784724616993690E+01, 8.833225342557650E+00,-1.984072076144117E-01), vec3(-1.773481248035921E-03, 3.341641213902314E-03, 3.532859459458915E-05)*365.242199, 8.8e25/M_sun, "Uranus");
    solarSystem.createCelestialBody(vec3(2.862016046630078E+01,-8.800209679340895E+00,-4.783572794559496E-01), vec3(9.014315108950128E-04, 3.018736862613127E-03,-8.311114083455642E-05)*365.242199, 1.03e26/M_sun, "Neptune");

    // To get a list (a reference, not copy) of all the bodies in the solar system, we use the .bodies() function
    vector<CelestialBody> &bodies = solarSystem.bodies();


    for(int i = 0; i < (int)bodies.size(); i++) {

        CelestialBody &body = bodies[i]; // Reference to this body
        solarSystem.calculateForcesAndEnergy();
        cout << "The position of " << body.objectName << " is " << body.position << " with velocity " << body.velocity << " and angular momentum " << body.angularMomentum << endl;
    }

    solarSystem.calculateForcesAndEnergy();
    cout << "Total energy before: " << solarSystem.totalEnergy() << endl;
    cout << "Total angular momentum before: " << solarSystem.angularMomentum() << endl;

    double dt = 0.001;
    Verlet integrator(dt);
    for(int timestep=0; timestep<numTimesteps; timestep++) {
        integrator.integrateOneStep(solarSystem);
        solarSystem.writeToFile("positions_f.xyz");
        //cout << "Total energy for iteration "<< timestep << " : " << solarSystem.totalEnergy() << endl;
        //cout << "Total angular momentum: " << solarSystem.angularMomentum() << endl;

    }

    solarSystem.calculateForcesAndEnergy();
    cout << "Total energy after: " << solarSystem.totalEnergy() << endl;
    cout << "Total angular momentum after: " << solarSystem.angularMomentum() << endl;




    cout << "I just created a solar system that has " << solarSystem.bodies().size() << " objects." << endl;
    //cout << "The potential energy of the system is "<< solarSystem.potentialEnergy()<<endl;
    cout << "The program simulates " << numTimesteps*dt << " years." << endl;
    return 0;
}
