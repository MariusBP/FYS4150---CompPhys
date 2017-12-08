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
    //for opening a file
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

}

void StatisticsSampler::sample(System &system, double initialTemperature)
{
    // different statistical properties written to the file
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleDiffusion(system);
    saveToFile(system, initialTemperature);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {
        m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
    }
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_totalPotentialEnergy = system.potential().totalPotentialEnergy();//steal the potential energy from the system class
}

void StatisticsSampler::sampleTemperature(System &system)
{
    double atomNumber = system.atoms().size();
    m_temperature = 2.0/3.0*(m_kineticEnergy/atomNumber); //compute temperature fro kinetic energy
}

void StatisticsSampler::sampleDensity(System &system)
{
    m_density = system.m_mass/system.volume();//The dansity of the system
}

void StatisticsSampler::sampleDiffusion(System &system)
{
    //for calculating the diffusion constant
    m_diffusionConstant = 0;
    double meanR = 0;
    for(Atom *atom: system.atoms())
    {
        vec3 r = atom->realPosition - atom->initialPosition;
        meanR += r.lengthSquared();//adding the length squared between all the particles and their respective initial position
    }

    m_diffusionConstant = (meanR/(double) system.atoms().size())/(6*system.time());//taking the mean length and dividing by 6*the time
}
