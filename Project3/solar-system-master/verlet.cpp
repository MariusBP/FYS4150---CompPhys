#include "verlet.h"
#include "solarsystem.h"

Verlet::Verlet(double dt):
m_dt(dt)
{};

void Verlet::integrateOneStep(SolarSystem &system)
{
    system.calculateForcesAndEnergy();

    for(CelestialBody &body : system.bodies()) {
        vec3 a = body.force / body.mass;
        body.velocity += a* m_dt;
        body.position += body.velocity*m_dt;
    }
};
