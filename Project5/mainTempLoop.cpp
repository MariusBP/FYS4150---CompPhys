#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main(int numberOfArguments, char **argumentList)
{
    int numberOfUnitCells = 5;
    double initialTemperature = 0;
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms
    ofstream m_file;

    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));

    double dt = UnitConverter::timeFromSI(1e-14); // Measured in seconds.
    int totalTimeSteps = 100;

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;

    cout << setw(20) << "Timestep" <<
            setw(20) << "Time" <<
            setw(20) << "Temperature" <<
            setw(20) << "KineticEnergy" <<
            setw(20) << "PotentialEnergy" <<
            setw(20) << "TotalEnergy" <<
            setw(20) << "Density" <<
            setw(20) << "Diffusion Constant" << endl;
    m_file.open("diffusion.txt", ofstream::out);
    for(double temperature =0; temperature<1000; temperature+=20){
        System system;
        initialTemperature = UnitConverter::temperatureFromSI(temperature); // measured in Kelvin
        cout << "Initial Temperature = " << initialTemperature << endl;
        system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
        system.potential().setEpsilon(UnitConverter::temperatureFromSI(119.8));
        system.potential().setSigma(3.405);

        system.removeTotalMomentum();
        cout << "Total momentum is " << system.m_momentum << setprecision(8)<<endl;



    StatisticsSampler statisticsSampler;
//    IO movie("movie.xyz"); // To write the state to file

        for(int timestep=0; timestep<totalTimeSteps; timestep++) {
            system.step(dt); //Update positions of the particles
            if(timestep == 1){
    statisticsSampler.sample(system);

    m_file << setw(20) << UnitConverter::temperatureToSI(statisticsSampler.temperature());
    cout << setw(20) << system.steps() <<
            setw(20) << UnitConverter::timeToSI( system.time()) <<
            setw(20) << UnitConverter::temperatureToSI(statisticsSampler.temperature()) <<
            setw(20) << UnitConverter::energyToEv(statisticsSampler.kineticEnergy()) <<
            setw(20) << UnitConverter::energyToEv(statisticsSampler.potentialEnergy()) <<
            setw(20) << UnitConverter::energyToEv(statisticsSampler.totalEnergy()) <<
            setw(20) << UnitConverter::massToSI(statisticsSampler.density())/1e-30 <<
            setw(20) << UnitConverter::diffusionToSI(statisticsSampler.diffusionConstant()) << endl;

            }
/*            statisticsSampler.sample(system);
                    if( timestep % 100 == 0 ) {
            // Print the timestep every 100 timesteps
            cout << setw(20) << system.steps() <<
                    setw(20) << UnitConverter::timeToSI( system.time()) <<
                    setw(20) << UnitConverter::temperatureToSI(statisticsSampler.temperature()) <<
                    setw(20) << UnitConverter::energyToEv(statisticsSampler.kineticEnergy()) <<
                    setw(20) << UnitConverter::energyToEv(statisticsSampler.potentialEnergy()) <<
                    setw(20) << UnitConverter::energyToEv(statisticsSampler.totalEnergy()) <<
                    setw(20) << UnitConverter::massToSI(statisticsSampler.density())/1e-30 <<
                    setw(20) << UnitConverter::diffusionToSI(statisticsSampler.diffusionConstant()) << endl;
*/        }

        //        movie.saveState(system);
        //        statisticsSampler.saveToFile(system);

    statisticsSampler.sample(system);
    m_file << setw(20) << UnitConverter::temperatureToSI(statisticsSampler.temperature())<<
              setw(20) << UnitConverter::diffusionToSI(statisticsSampler.diffusionConstant()) <<
              setw(20) << UnitConverter::timeToSI( system.time()) << endl;

    cout << setw(20) << system.steps() <<
            setw(20) << UnitConverter::timeToSI( system.time()) <<
            setw(20) << UnitConverter::temperatureToSI(statisticsSampler.temperature()) <<
            setw(20) << UnitConverter::energyToEv(statisticsSampler.kineticEnergy()) <<
            setw(20) << UnitConverter::energyToEv(statisticsSampler.potentialEnergy()) <<
            setw(20) << UnitConverter::energyToEv(statisticsSampler.totalEnergy()) <<
            setw(20) << UnitConverter::massToSI(statisticsSampler.density())/1e-30 <<
            setw(20) << UnitConverter::diffusionToSI(statisticsSampler.diffusionConstant()) << endl;

}
//    movie.close();
    m_file.close();
    return 0;
}
