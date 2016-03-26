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

float CalculatePoly6(float r, float h){
    return (315.f * pow((h*h - r*r),3) / (64.f * M_PI * pow(h, 9)));
}

float CalculateGradSpiky(float r, float h){
    return -45.f * pow((h - r), 2) / (M_PI * pow(h, 6));
}

// Given position in world space, return the index of the corresponding cell in grid space.
vec3 Grid::getIndices(vec3 pos) {
    int i = (int)(pos[0] / cellSize);
    int j = (int)(pos[1] / cellSize);
    int k = (int)(pos[2] / cellSize);
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

std::vector<Particle *> PBFSolver::FindNeighbors(Particle *p){
    p->neighbors.clear();
    vec3 indices = uGrid.getIndices(*p);
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
        float r = glm::distance(n->pos, p->pos);
        totalDensity += CalculatePoly6(r, h);
    }
    return totalDensity;
}

