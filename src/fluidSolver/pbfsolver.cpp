#include "pbfsolver.h"

PBFSolver::PBFSolver(vec3 minBounds, vec3 maxBounds)
{
    rDensity = 1000.f;
    numParticles = ParticlesContainer.size();
    uGrid = Grid(minBounds, maxBounds, 1);
}

void PBFSolver::ApplyForces(Particle& p, float del_t){
    //gravity force
    p.speed += glm::vec3(0.f, -9.81f, 0.f) * del_t;
    p.pos += p.speed * del_t;
}

float CalculatePoly6(vec3 r, float h){
    float r2 = pow(length(r),2);
    return (315.f * pow((h*h - r2),3) / (64.f * M_PI * pow(h, 9)));
}

vec3 CalculateGradSpiky(vec3 r, float h){
    float mag_r = (length(r));
    return -45.f * pow((h - mag_r), 2) / (M_PI * pow(h, 6)) * (r / mag_r);
}

// Given position in world space, return the index of the corresponding cell in grid space.
vec3 Grid::getIndices(const vec3 &pos) {
    int i = (int)((pos[0] - minBounds[0]) / cellSize);
    int j = (int)((pos[1] - minBounds[1]) / cellSize);
    int k = (int)((pos[2] - minBounds[2]) / cellSize);
    return vec3(i, j, k);
}

int Grid::flatIndex(int i, int j, int k) {
    return i + dimensions[1] * (j + dimensions[2] * k);
}

Cell* Grid::operator() (int i, int j, int k) {
    return cells[flatIndex(i, j, k)];
}

void Grid::update(const std::vector<Particle> &particles) {
    for (Particle p : particles) {
        vec3 indices = getIndices(p.pos);
        Cell * cell = cells[flatIndex(indices[0], indices[1], indices[2])];
        cell->particles.push_back(&p);
        cell->count++;
    }
}

void Grid::clear() {
    cells = std::vector<Cell *>(dimensions[0]*dimensions[1]*dimensions[2], new Cell());
}

//find all the neighbors of this particle. NOTE: this includes the current particle in consideration
void PBFSolver::FindNeighbors(Particle *p){
    p->neighbors.clear();
    vec3 indices = uGrid.getIndices(p->pos);
    for (int i = indices[0] - NEIGHBOR_RADIUS; i <= indices[0] + NEIGHBOR_RADIUS; ++i){
        for (int j = indices[1] - NEIGHBOR_RADIUS; j <= indices[1] + NEIGHBOR_RADIUS; ++j){
            for (int k = indices[2] - NEIGHBOR_RADIUS; k <= indices[2] + NEIGHBOR_RADIUS; ++k){
                if (i > 0 && j > 0 && k > 0){
                    if (i < uGrid.dimensions[0] && j < uGrid.dimensions[1] && k < uGrid.dimensions[2]){
                        std::vector<Particle *> tmp_cell = uGrid(i, j, k)->particles;
                        p->neighbors.insert( p->neighbors.end(), tmp_cell.begin(), tmp_cell.end() );
                    }
                }
            }
        }
    }
}

float PBFSolver::CalculateDensity(Particle *p, float h){
    float totalDensity = 0.f;

    for(Particle *n : p->neighbors){
        vec3 r = ((n->pos - p->pos));
        totalDensity += CalculatePoly6(r, h);
    }
    // this is because we have the current particle in the neighbor list, so we subtract the contribution to density
    return totalDensity - CalculatePoly6(vec3(0.f), h);
}

vec3 PBFSolver::GradientAtN(Particle* n, Particle* p){
    vec3 r = (p->pos - n->pos);
    float h = 1.f;
    return -1.f / rDensity * CalculateGradSpiky(r, h);
}

vec3 PBFSolver::GradientAtP(Particle* p){
    vec3 totalGradientConstraint(0.f);
    float h = 1.f;
    for(Particle *n : p->neighbors){
        vec3 r = ((n->pos - p->pos));
        totalGradientConstraint += CalculateGradSpiky(r, h);
    }
    // this is because we have the current particle in the neighbor list, so we subtract the contribution to density
    return 1.f / rDensity * (totalGradientConstraint - CalculateGradSpiky(vec3(0.f), h));
}

float PBFSolver::CalculateLagrangeMultiplier(Particle* p){
    const float eps = 1.0e-6f;

    float h = 1.f, lambda;
    float density = CalculateDensity(p, h);
    // Evaluate constraint function, numerator for the lagrangian multiplier
    const float C = std::max(density / rDensity - 1.0f, 0.0f);

    if(C != 0.f){
        float sum_grad_C2 = 0.0;
        vec3 gradC_i(0.0f, 0.0f, 0.0f);

        FindNeighbors(p);

        for(Particle* n : p->neighbors){
            //1. if k = i, k is the current particle
                //summation over the neighbors
            if(length(p->pos - n->pos) < EPSILON){
                continue;
            }
            else{
                sum_grad_C2 += length2(GradientAtN(n, p));
            }
            //2. if k = j, k is a neighboring particle
                //negative grad
        }
        sum_grad_C2 += length2(GradientAtP(p));

        // Compute lambda
        lambda = -C / (sum_grad_C2 + eps);
    }
    else{
        lambda = 0.f;
    }
    return lambda;
}

vec3 PBFSolver::CalculateDeltaP(Particle *p, int particleIndex, int neighborIndex){
    vec3 del_p(0.f);
    FindNeighbors(p);

    const unsigned int numberOfParticles = p->neighbors.size();
    for (unsigned int j = 0; j < numberOfParticles; j++)
    {
        if (neighborIndex < numberOfParticles)		// Test if fluid particle
        {
            const vec3 gradC_j = -1.f / rDensity * CalculateGradSpiky(p->pos - p->neighbors.at(j)->pos, 1.f);
            del_p -= (lambda[particleIndex] + lambda[neighborIndex]) * gradC_j;
        }
    }
    return del_p;
}

