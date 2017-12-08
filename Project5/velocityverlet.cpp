#include "velocityverlet.h"
#include "system.h"
#include "atom.h"

void VelocityVerlet::integrate(System &system, double dt)
{
    if(m_firstStep) {
        system.calculateForces();
        m_firstStep = false;
    }

    for(Atom *atom : system.atoms()) {
        atom->velocity += atom->force*dt/atom->mass();
        vec3 change = atom->velocity*dt;
        atom->position += change;
        atom->realPosition += change; //Updating the "actual" position of the atom
    }
    system.applyPeriodicBoundaryConditions();
    system.calculateForces(); // New positions, recompute forces
}
