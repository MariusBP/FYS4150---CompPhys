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

void StatisticsSampler::open(const char *filename)
{
    m_file.open(filename);
}

void StatisticsSampler::saveToFile(System &system, double initialTemperature)
{
    // Save the statistical properties for each timestep for plotting etc.
        m_file << setw(20) << system.steps() <<
            setw(20) << UnitConverter::timeToSI(system.time()) <<
            setw(20) << UnitConverter::temperatureToSI(initialTemperature) <<
            setw(20) << UnitConverter::temperatureToSI(temperature()) << " "<<
            setw(20) << UnitConverter::energyToEv(kineticEnergy()) << " "<<
            setw(20) << UnitConverter::energyToEv(potentialEnergy()) << " "<<
            setw(20) << UnitConverter::energyToEv(totalEnergy()) << endl;

    // Print out values here
}

void StatisticsSampler::sample(System &system, double initialTemperature)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleDiffusion(system);
    saveToFile(system, initialTemperature);
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
    double meanR = 0;
    for(Atom *atom: system.atoms())
    {
        vec3 r = atom->realPosition - atom->initialPosition;
        meanR += r.lengthSquared();
    }

    m_diffusionConstant = (meanR/(double) system.atoms().size())/(6*system.time());
}
