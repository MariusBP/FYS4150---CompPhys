#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "solarsystem.h"
#include "verlet.h"
using namespace std;

int main(int numArguments, char **arguments)
{
    int numTimesteps = 10000;
    if(numArguments >= 2) numTimesteps = atoi(arguments[1]);

    SolarSystem solarSystem;
    // We create new bodies like this. Note that the createCelestialBody function returns a reference to the newly created body
    // This can then be used to modify properties or print properties of the body if desired
    // Use with: solarSystem.createCelestialBody( position, velocity, mass );

    solarSystem.createCelestialBody( vec3(0,0,0), vec3(0,0,0), 1.0 );

    // We don't need to store the reference, but just call the function without a left hand side
    solarSystem.createCelestialBody( vec3(1, 0, 0), vec3(0, 2*M_PI, 0), 3e-6 );//Earth
    solarSystem.createCelestialBody(vec3(5.2, 0, 0), vec3(0,2*M_PI*5.2/11.8,0), 1e-3*1000);//Jupiter


    // To get a list (a reference, not copy) of all the bodies in the solar system, we use the .bodies() function
    vector<CelestialBody> &bodies = solarSystem.bodies();

    for(int i = 0; i<bodies.size(); i++) {
        CelestialBody &body = bodies[i]; // Reference to this body
        cout << "The position of this object is " << body.position << " with velocity " << body.velocity << endl;
    }

    double dt = 0.001;
    Verlet integrator(dt);
    for(int timestep=0; timestep<numTimesteps; timestep++) {
        integrator.integrateOneStep(solarSystem);
        solarSystem.writeToFile("positions.xyz");
    }

    cout << "I just created a solar system that has " << solarSystem.bodies().size() << " objects." << endl;
    cout << "The potential energy of the system is "<< solarSystem.potentialEnergy()<<endl;
    cout << "The program ran for " << dt*numTimesteps<<" yrs."<<endl;
    return 0;
}
