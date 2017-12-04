#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using std::ofstream; using std::cout; using std::endl; using std::setw;

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system)
{
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already
/*    if(!m_file.good()) {
        m_file.open("statistics.txt", ofstream::out);
        // If it's still not open, something bad happened...
        if(!m_file.good()) {
            cout << "Error, could not open statistics.txt" << endl;
            exit(1);
        }*/
        m_file.open("statistics.txt", ofstream::out);
        m_file << setw(20) << system.steps() <<
            setw(20) << system.time() <<
            setw(20) << UnitConverter::temperatureToSI(temperature()) << " "<<
            setw(20) << UnitConverter::energyToEv(kineticEnergy()) << " "<<
            setw(20) << UnitConverter::energyToEv(potentialEnergy()) << " "<<
            setw(20) << UnitConverter::energyToEv(totalEnergy()) << "/n";

    // Print out values here
}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleDiffusion(system);
    saveToFile(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {
        m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
    }
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_totalPotentialEnergy = system.potential().totalPotentialEnergy();
}

void StatisticsSampler::sampleTemperature(System &system)
{
    double atomNumber = system.atoms().size();
    m_temperature = 2.0/3.0*(m_kineticEnergy/atomNumber);
}

void StatisticsSampler::sampleDensity(System &system)
{
    m_density = system.m_mass/system.volume();
}

void StatisticsSampler::sampleDiffusion(System &system)
{
    m_diffusionConstant = 0;
    for(Atom *atom: system.atoms())
    {
        vec3 r = atom->position - atom->initialPosition;
        m_diffusionConstant += r.lengthSquared();
    }

    m_diffusionConstant = (m_diffusionConstant/(double) system.atoms().size())/(6*system.time());
}
