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

// Account for periodic boundary conditions of the system of atoms
void System::applyPeriodicBoundaryConditions()
{
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    double  x_size  = m_systemSize.x();
    double  y_size  = m_systemSize.y();
    double  z_size  = m_systemSize.z();
    for(Atom *atom : m_atoms)
    {
        // Update x,y, and z dimensions
        atom->position.setX(atom->position.x() - floor(atom->position.x() / x_size) * x_size);
        atom->position.setY(atom->position.y() - floor(atom->position.y() / y_size) * y_size);
        atom->position.setZ(atom->position.z() - floor(atom->position.z() / z_size) * z_size);
    }
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
    setSystemSize(vec3(0, 0, 0) + numberOfUnitCellsEachDimension * latticeConstant); //Setting system size N_cells*b

    //different positions of the atoms in a given unit cell:
    vec3 a(0,0,0);
    vec3 b(latticeConstant/2,latticeConstant/2,0);
    vec3 c(0,latticeConstant/2,latticeConstant/2);
    vec3 d(latticeConstant/2,0,latticeConstant/2);

    vec3 latticePositions[4] = {a,b,c,d};

    Random::randomSeed();
    Random::nextGaussian(1.0, 0.5);

    for(int i=0; i < numberOfUnitCellsEachDimension; i++){
        for(int j=0; j < numberOfUnitCellsEachDimension; j++){
            for(int k=0; k < numberOfUnitCellsEachDimension; k++){
                vec3 R(i*latticeConstant,j*latticeConstant, k*latticeConstant);//Loop
                for(vec3 r: latticePositions){
                    Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26)); //Argon weight converted to ~40 m_u
                    atom->position = R + r; //Position of the atom we are going to place
                    atom->realPosition = R + r; //Position of the atom with no boundaries applied
                    atom->initialPosition = R + r; //setting initial position for use in computing diffusionConstant
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
    for(int i=0; i< (int) m_atoms.size(); i++) //For all atoms
    {
        Atom *atom1 = m_atoms[i];
        for(int j=i+1; j< (int) m_atoms.size(); j++) //For all atoms other than the ones already looped thorugh
        {
            Atom *atom2 = m_atoms[j];
            m_potential.calculateForce(*this, *atom1, *atom2); // Calculate and assign the force from atom1 to atom2 and vise versa
        }
    }
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}
