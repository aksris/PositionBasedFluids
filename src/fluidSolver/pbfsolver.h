#ifndef PBFSOLVER_H
#define PBFSOLVER_H
#include "fluidSolver.hpp"
#include <map>
# define M_PI          3.141592653589793238462643383279502884f /* pi */
#define EPSILON 0.001f;
#define del_q 0.1
#define PBF_H 0.5f
#define NEIGHBOR_RADIUS 1 //find neighbor radius

class Cell
{
public:
    Cell() : count(0) {}
    std::vector<Particle *> particles;
    int count;
};

class Grid
{
public:
    Grid() {}
    Grid(vec3 min_bounds, vec3 max_bounds, float cell_size) :
        minBounds(min_bounds), maxBounds(max_bounds), cellSize(cell_size) {
        dimensions = vec3((maxBounds[0]-minBounds[0])/cellSize,
                (maxBounds[1]-minBounds[1])/cellSize, (maxBounds[2]-minBounds[2])/cellSize);
        cells = std::vector<Cell *>(dimensions[0]*dimensions[1]*dimensions[2], new Cell());
    };
    vec3 minBounds;
    vec3 maxBounds;
    vec3 dimensions;
    float cellSize;
    std::vector<Cell *> cells;

    Cell* operator() (int i, int j, int k);

    vec3 getIndices(vec3 pos);
    void update(const std::vector<Particle> &particles);
    void clear();
    int flatIndex(int i, int j, int k);
};

class PBFSolver : public FluidSolver
{
public:
    PBFSolver(vec3 minBounds, vec3 maxBounds);

    void ApplyForces(Particle &p, float del_t);
    float CalculateDensity(Particle *p, float h);
    std::vector <Particle*> FindNeighbors(Particle *p);
protected:
    float rDensity;
private:
    int numParticles;
    Grid uGrid;
};
float CalculateGradSpiky(float r, float h);
float CalculatePoly6(float r, float h);
#endif // PBFSOLVER_H
