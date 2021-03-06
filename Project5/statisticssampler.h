#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H
#include <fstream>

class System; // Promise the compiler that this is a class even though we haven't included system.h here

class StatisticsSampler
{
private:
    std::ofstream m_file;
    double m_kineticEnergy = 0;
    double m_totalPotentialEnergy = 0;
    double m_temperature = 0;
    double m_density = 0;
    double m_diffusionConstant = 0;
public:
    StatisticsSampler();
    void open(const char *filename);
    void saveToFile(System &system, double initialTemperature);
    void sample(System &system, double initialTemperature);
    void sampleKineticEnergy(System &system);
    void samplePotentialEnergy(System &system);
    void sampleTemperature(System &system);
    void sampleDensity(System &system);
    void sampleDiffusion(System &system);
    double kineticEnergy() { return m_kineticEnergy; }
    double potentialEnergy() { return m_totalPotentialEnergy; }
    double totalEnergy() { return m_kineticEnergy+m_totalPotentialEnergy; }
    double temperature() { return m_temperature; }
    double density() { return m_density; }
    double diffusionConstant() {return m_diffusionConstant; }
};
#endif
