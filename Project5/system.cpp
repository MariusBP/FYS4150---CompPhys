#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    //Sizes of the system in x,y and z direction
    double size_x = m_systemSize.x();
    double size_y = m_systemSize.y();
    double size_z = m_systemSize.z();

    for (Atom *atom: m_atoms){
        //Looping the x-position
        if (atom->position.x() <  0){
            atom->initialPosition.setX(atom->initialPosition.x() + size_x);
            atom->position.setX(atom->position.x() + size_x);
         }
        if (atom->position.x() >=  size_x){
            atom->initialPosition.setX(atom->initialPosition.x() - size_x);
            atom->position.setX( atom->position.x() - size_x);
         }
        //Looping y-position
        if (atom->position.y() <  0){
            atom->initialPosition.setY(atom->initialPosition.y() + size_y);
            atom->position.setY(atom->position.y() + size_y);
            }
        if (atom->position.y() >=  size_y){
            atom->initialPosition.setY(atom->initialPosition.y() - size_y);
            atom->position.setY( atom->position.y() - size_y);
            }
        //Looping z-position
        if (atom->position.z() <  0){
            atom->initialPosition.setZ(atom->initialPosition.z()+ size_z);
            atom->position.setZ(atom->position.z() + size_z);
         }
        if (atom->position.z() >=  size_z){
            atom->initialPosition.setZ(atom->initialPosition.z()-size_z);
            atom->position.setZ( atom->position.z() - size_z);
            }
        };
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
}

void System::removeTotalMomentum() {
    m_momentum.zeros();
    for (Atom *atom: m_atoms){
        m_momentum += atom->mass()*atom->velocity;
    }
    double atomNumber = m_atoms.size();
    vec3 scaling_factor = m_momentum/atomNumber;

    for (Atom *atom: m_atoms){
        atom->velocity -= scaling_factor/atom->mass();
    }
    m_momentum.zeros();
    for (Atom *atom: m_atoms){
        m_momentum += atom->mass()*atom->velocity;
    }
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You implemented this function properly.
    setSystemSize(vec3(0, 0, 0) + numberOfUnitCellsEachDimension * latticeConstant); // Remember to set the correct system size!
/*
    Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26)); //Argon weight converted to ~40 m_u

    atom1->position.set(0,0,0);
    atom1->velocity = vec3(0,0,0);
    m_atoms.push_back(atom1);

    Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26)); //Argon weight converted to ~40 m_u
    atom2->position.set(1,0,0);
    atom2->velocity = vec3(0,0,0);
    m_atoms.push_back(atom2);

*/
    //different positions of the atoms in a unit cell
    vec3 a(0,0,0);
    vec3 b(latticeConstant/2,latticeConstant/2,0);
    vec3 c(0,latticeConstant/2,latticeConstant/2);
    vec3 d(latticeConstant/2,0,latticeConstant/2);

    vec3 latticePositions[4] = {a,b,c,d};

    for(int i=0; i < numberOfUnitCellsEachDimension; i++){
        for(int j=0; j < numberOfUnitCellsEachDimension; j++){
            for(int k=0; k < numberOfUnitCellsEachDimension; k++){
                vec3 R(i*latticeConstant,j*latticeConstant, k*latticeConstant);//Loop
                for(vec3 r: latticePositions){
                    Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26)); //Argon weight converted to ~40 m_u
                    atom->position = R + r; //position of the atom we are going to place
                    atom->initialPosition = atom->position;
                    atom->resetVelocityMaxwellian(temperature); //set initial velocity due to a gaussian distribution
                    m_atoms.push_back(atom);
                    m_mass += atom->mass(); //add atom mass to total mass
                }
            }
        }
    }
}
    /*
    for(int i=0; i<100; i++) {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26)); //Argon weight converted to ~40 m_u
        double x = Random::nextDouble(0, 10); // random number in the interval [0,10]
        double y = Random::nextDouble(0, 10);
        double z = Random::nextDouble(0, 10);
        atom->position.set(x,y,z);
        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);
    }
    setSystemSize(vec3(10, 10, 10)*/

void System::calculateForces()
{
    m_potential.m_totalPotentialEnergy = 0; //setting the total potential energy to 0
    for(Atom *atom : m_atoms) {
        atom->resetForce();
        }
    for(int i=0; i< (int) m_atoms.size(); i++)
    {
        Atom *atom1 = m_atoms[i];
        for(int j=i+1; j< (int) m_atoms.size(); j++)
        {
            Atom *atom2 = m_atoms[j];
            m_potential.calculateForce(*this, *atom1, *atom2); // Calculate and assign the force from atom1 to atom2 and vise versa
        }
    }
}
/*
 * void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}
*/

void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}
