#ifndef LENNARDJONES_H
#define LENNARDJONES_H

class LennardJones
{
private:
    double m_sigma = 1.0;
    double m_epsilon = 1.0;
    double m_potentialEnergy = 0;

public:
    double m_totalPotentialEnergy = 0;
    LennardJones() { }
    double totalPotentialEnergy() const;
    double potentialEnergy() const;
    double sigma() const;
    void setSigma(double sigma);
    double epsilon() const;
    void setEpsilon(double epsilon);
    void calculateForce(class System &system, class Atom &atom1, class Atom &atom2);
    void calculatePotentialEnergy(class Atom &atom1, class Atom &atom2);
};
#endif
