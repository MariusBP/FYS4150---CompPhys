#include "lennardjones.h"
#include "system.h"
#include <cmath>

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::totalPotentialEnergy() const
{
    return m_totalPotentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculatePotentialEnergy(Atom &atom1,Atom &atom2)
{
    vec3 r12           = atom2.position - atom1.position;
    double rl = r12.length();
    double R12 = sigma()/rl;
    double R6 = sigma()/rl;
    for(int i=0;i<11; i++){
        R12 *= sigma()/rl;
    }
    for(int i=0; i<5; i++){
        R6 *=sigma()/rl;
    }
    m_potentialEnergy  = 4.0*epsilon()*(R12 - R6); //calculating the potential U between two atoms
    m_totalPotentialEnergy += m_potentialEnergy; // Add the E_p between atom1 and atom2 to the total E_p
}

void LennardJones::calculateForce(System &system, Atom &atom1,Atom &atom2)
{
    vec3 Force12(0,0,0);
/*    double x1x2 = atom2.position.x() - atom1.position.x();
    double y1y2 = atom2.position.y() - atom1.position.y();
    double z1z2 = atom2.position.z() - atom1.position.z();
    if(abs(x1x2) > system.systemSize().x()*0.5){
        if(x1x2 < 0){
            x1x2 += system.systemSize().x()*0.5;
        }
        else{
            x1x2 -= system.systemSize().x()*0.5;
        }
    }
    if(abs(y1y2) > system.systemSize().y()*0.5){
        if(y1y2 < 0){
            y1y2 += system.systemSize().y()*0.5;
        }
        else{
            y1y2 -= system.systemSize().y()*0.5;
        }
    }
    if(abs(z1z2) > system.systemSize().z()*0.5){
        if(z1z2 < 0){
            z1z2 += system.systemSize().z()*0.5;
        }
        else{
            z1z2 -= system.systemSize().z()*0.5;
        }
    }
*/
    vec3 r12 = atom2.position - atom1.position;//(x1x2,y1y2, z1z2);
    double rl = r12.length();
    double R12 = sigma()/rl;
    double R6 = sigma()/rl;
    for(int i=0;i<11; i++){
        R12 *= sigma()/rl;
    }
    for(int i=0; i<5; i++){
        R6 *=sigma()/rl;
    }

    m_potentialEnergy  = 4.0*epsilon()*(R12 - R6); //calculating the potential U between two atoms
    m_totalPotentialEnergy += m_potentialEnergy; // Add the E_p between atom1 and atom2 to the total E_p
    Force12 = 24*epsilon()*(2*R12 - R6)*(r12/r12.lengthSquared());
    atom1.force -= Force12;
    atom2.force += Force12;
}
