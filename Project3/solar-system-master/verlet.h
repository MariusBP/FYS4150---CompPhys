#ifndef VERLET_H
#define VERLET_H

class Verlet
{
public:
    double m_dt;
    Verlet(double dt);
    void integrateOneStep(class SolarSystem &system);
};

#endif // VERLET_H
