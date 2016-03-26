#ifndef PBFSOLVER_H
#define PBFSOLVER_H
#include "fluidSolver.hpp"
#include <map>
# define M_PI          3.141592653589793238462643383279502884f /* pi */
#define EPSILON 0.001f;
#define del_q 0.1
#define PBF_H 0.5f
class PBFSolver : public FluidSolver
{
public:
    PBFSolver();

    void ApplyForces(Particle &p, float del_t);
    void initializeGrid();
    float CalculateLambda(Particle p);
    vec3 CalculateDeltaP(Particle p, int i);
protected:
    float CalculateSCorr(Particle p);
    float CalculateDensity(Particle p);
    float CalculateGradientConstraint(Particle p, Particle k);
    void NeighborSearch(Particle *p, std::vector<Particle> &neighbors);
    float rDensity;
    std::map<int, Particle*> iNeighbors;
private:
    int numParticles;
};
float CalculateSpiky(float r, float h);
float CalculatePoly6(float r, float h);
#endif // PBFSOLVER_H
