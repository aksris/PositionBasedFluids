#include "pbfsolver.h"

#define MU_VISCOSITY 0.02f
#define MAX_ITER 5
#define CELL_SIZE 0.5f

PBFSolver::PBFSolver(vec3 minBounds, vec3 maxBounds)
{
    rDensity = 0.9f;
    numParticles = ParticlesContainer.size();
    uGrid = Grid(minBounds, maxBounds, CELL_SIZE);
}

void PBFSolver::step(){
    //clearing acceleration and set gravity
    clearAccel();
    // Advection; semi implicit euler
    for(Particle &p: ParticlesContainer){
        p.pos_star = vec3(0.f);
        p.old_pos = p.pos;
        p.speed += p.accel * t_stp;
        p.pos += p.speed * t_stp;
    }
    // Clear grid of particles.
    uGrid.clear();
    // Update grid cells.
    uGrid.update(ParticlesContainer);
    //find all particle neighbors
    for(Particle &p: ParticlesContainer){
        FindNeighbors(&p);
    }
    //constraint projection
    ConstraintProjection();
    for(Particle &p: ParticlesContainer){
        p.speed = (1 / t_stp) * (p.pos - p.old_pos);
    }
    //add viscosity
    //boundary
    handleBoundary();
}

void PBFSolver::handleBoundary() {
            // Bounds clamping
    float x_dim = uGrid.dimensions[0] / 2 * CELL_SIZE;
    float y_dim = uGrid.dimensions[1] / 2 * CELL_SIZE;
    float z_dim = uGrid.dimensions[2] / 2 * CELL_SIZE;
    float offset = 0.001;
    
    for(Particle &p : ParticlesContainer){
        if (p.pos[0] < -x_dim || p.pos[0] > x_dim) {
            p.pos[0] = std::min(x_dim - offset, std::max(-x_dim + offset, p.pos[0]));
            p.speed[0] *= -0.7f;
        }
        if (p.pos[1] < -y_dim || p.pos[1] > y_dim) {
            p.pos[1] = std::min(y_dim - offset, std::max(-y_dim + offset, p.pos[1]));
            p.speed[1] *= -0.7f;
        }
        if (p.pos[2] < -z_dim || p.pos[2] > z_dim) {
            p.pos[2] = std::min(z_dim - offset, std::max(-z_dim + offset, p.pos[2]));
            p.speed[2] *= -0.7f;
        }
    }

}

void PBFSolver::ConstraintProjection(){
    int iter = 0;
    while(iter < MAX_ITER){
        for(Particle &p : ParticlesContainer){
            p.density = CalculateDensity(&p, uGrid.cellSize);
            CalculateLagrangeMultiplier(&p);
        }
        for(Particle &p: ParticlesContainer){
            vec3 corr = SolveDensityConstraint(&p);
            p.pos_star = corr;
        }
        for(Particle &p : ParticlesContainer){
            p.pos += p.pos_star;
        }
        ++iter;
    }
}

vec3 PBFSolver::SolveDensityConstraint(Particle *p){
    vec3 corr(0.f);
    for(Particle *n : p->neighbors){
        vec3 grad_w = GradW(p->pos - n->pos);
        vec3 grad_Cj = grad_w * -n->mass / rDensity;
        corr -= (p->lambda + n->lambda) * grad_Cj;
    }
    return corr;
}

float PBFSolver::W(vec3 r){

    float res = 0.f;
    const float r1 = sqrt(dot(r, r));
//            normalize(r);
    const float q = r1 / uGrid.cellSize;

    if(q <= 1.0001f){
        if(q <= 0.5f){
            const float q2 = q * q;
            const float q3 = q * q * q;
            res = (8.f / M_PI * pow(uGrid.cellSize, 3)) * (6.f * q3 - 6.f * q2 + 1.f);
        }
        else{
            res = (8.f / M_PI * pow(uGrid.cellSize, 3)) * (2.f * pow(1.f - q, 3));
        }
    }
    return res;
}

vec3 PBFSolver::GradW(vec3 r){
    vec3 res(0.f);

    const float r1 = sqrt(dot(r, r));
//            normalize(r);
    const float q = r1 / uGrid.cellSize;

    if(q <= 1.0001f){
        if(r1 > 1.0e-6){
            const vec3 gradQ = r * ((float) 1.0 / (r1 * uGrid.cellSize));
            if(q <= 0.5f){
                res = (float)(48.f / (M_PI * pow(uGrid.cellSize, 3))) * q * ((3.0f * q - 2.0f)) * gradQ;
            }
            else{
                const float factor = 1.f - q;
                res = (float)(48.f / (M_PI * pow(uGrid.cellSize, 3))) * (-factor * factor) * gradQ;
            }
        }
    }

    return res;
}

void PBFSolver::clearAccel(){
    for(Particle &p : ParticlesContainer){
        p.accel = vec3(0.f, -9.81f, 0.f);
    }
}

void PBFSolver::ApplyForces(Particle& p, float del_t){
    //gravity force
    p.speed += glm::vec3(0.f, -9.81f, 0.f) * del_t;
    p.pos_star = p.pos + p.speed * del_t;
}

float CalculatePoly6(vec3 r, float h){
    float r2 = pow(length(r),2);
    return (315.f * pow((h*h - r2),3) / (64.f * M_PI * pow(h, 9)));
}

vec3 CalculateGradSpiky(vec3 r, float h){
    float mag_r = (length(r));
    if(mag_r < 0.00001f)
        return vec3(0.f);
    return -45.f * (float) (pow((h - mag_r), 2) / (M_PI * pow(h, 6))) * (r / mag_r);
}

// Given position in world space, return the index of the corresponding cell in grid space.
ivec3 Grid::getIndices(const vec3 &pos) {
    int i = std::min(std::max((int)((pos[0] - minBounds[0]) / cellSize), 0), (int)(dimensions[0] - 1.f));
    int j = std::min(std::max((int)((pos[1] - minBounds[1]) / cellSize), 0), (int)(dimensions[1] - 1.f));
    int k = std::min(std::max((int)((pos[2] - minBounds[2]) / cellSize), 0), (int)(dimensions[2] - 1.f));
    return ivec3(i, j, k);
}

int Grid::flatIndex(int i, int j, int k) {
    int a = i + dimensions[1] * (j + dimensions[2] * k);
    return a;
}

Cell* Grid::operator() (int i, int j, int k) {
    return cells[flatIndex(i, j, k)];
}

// Grid cells keep a list of the particles that they contain.
// This function populates cells with the correct particles.
void Grid::update( std::vector<Particle> &particles) {
    for (Particle &p : particles) {
        ivec3 indices = getIndices(p.pos);
        Cell * cell = cells[flatIndex(indices[0], indices[1], indices[2])];
        cell->particles.push_back(&p);
        cell->count++;
    }
}

// Clear all the grid cells
void Grid::clear() {
    for (Cell * c : cells) {
        c->particles.clear();
        c->count = 0;
    }
}

// Find all the neighbors of this particle.
void PBFSolver::FindNeighbors(Particle *p){
    p->neighbors.clear();
    ivec3 indices = uGrid.getIndices(p->pos);
    for (int i = indices[0] - NEIGHBOR_RADIUS; i <= indices[0] + NEIGHBOR_RADIUS; ++i){
        for (int j = indices[1] - NEIGHBOR_RADIUS; j <= indices[1] + NEIGHBOR_RADIUS; ++j){
            for (int k = indices[2] - NEIGHBOR_RADIUS; k <= indices[2] + NEIGHBOR_RADIUS; ++k){
                if (i > 0 && j > 0 && k > 0){
                    if (i < uGrid.dimensions[0] && j < uGrid.dimensions[1] && k < uGrid.dimensions[2]){
                        std::vector<Particle *> &tmp_cell = uGrid(i, j, k)->particles;
                        //p->neighbors.insert( p->neighbors.end(), tmp_cell.begin(), tmp_cell.end() );
                        for (unsigned int q = 0; q < tmp_cell.size(); ++q) {
                            if (tmp_cell.at(q) != p) {
                                p->neighbors.push_back(tmp_cell.at(q));
                            }
                        }
                    }
                }
            }
        }
    }
}

float CalculatePoly8(float r, float h){
    return (315.f * pow((h*h - r*r),3) / (64.f * M_PI * pow(h, 9)));
}

float PBFSolver::CalculateSCorr(vec3 r, float h, float q, int n){
    float ret = 1.f;

    ret = 0.1f * pow((CalculatePoly6(r, h)/CalculatePoly8(q*h, h)),n);

    return ret;
}

float PBFSolver::CalculateDensity(Particle *p, float h){
    float totalDensity = 0.f;
    totalDensity = p->mass * W(vec3(0.f));
    for(Particle *n : p->neighbors){
        totalDensity += n->mass * W(p->pos - n->pos);
    }
    // this is because we have the current particle in the neighbor list, so we subtract the contribution to density
    return totalDensity;
}

vec3 PBFSolver::GradientAtN(Particle* n, Particle* p){
    vec3 r = (p->pos - n->pos);
    float h = 1.000001f;
    return -1.f / rDensity * CalculateGradSpiky(r, h);
}

vec3 PBFSolver::GradientAtP(Particle* p){
    vec3 totalGradientConstraint(0.f);
    float h = 1.000001f;
    for(Particle *n : p->neighbors){
        vec3 r = ((n->pos - p->pos));
        totalGradientConstraint += CalculateGradSpiky(r, h);
    }
    // this is because we have the current particle in the neighbor list, so we subtract the contribution to density
    return 1.f / rDensity * (totalGradientConstraint );
}

void PBFSolver::CalculateLagrangeMultiplier(Particle* p){
    const float eps = 1.0e-6f;

    float h = 1.000001f, lambda;
//    float density = CalculateDensity(p, h);
    // Evaluate constraint function, numerator for the lagrangian multiplier
    const float C = std::max(p->density / rDensity - 1.0f, 0.0f);

    if(C != 0.f){
        float sum_grad_C2 = 0.0;
        vec3 gradC_i(0.0f);

        for(Particle* n : p->neighbors){
            //1. if k = i, k is the current particle
                //summation over the neighbors
            vec3 grad_Cj = -n->mass / rDensity * GradW(p->pos - n->pos);
            sum_grad_C2 += dot(grad_Cj, grad_Cj);
                //grad_Cj.squaredNorm()
            gradC_i -= grad_Cj;
//            if(length(p->pos - n->pos) < 0.001f){
//                continue;
//            }
//            else{
//                sum_grad_C2 += length2(GradientAtN(n, p));
//            }
//            //2. if k = j, k is a neighboring particle
//                //negative grad
//        }
        }
        sum_grad_C2 += dot(gradC_i, gradC_i);
            //grad_Ci.squaredNorm();

        // Compute lambda
        lambda = -C / (sum_grad_C2 + eps);
    }
    else{
        lambda = 0.f;
    }
    p->lambda = lambda;
}

vec3 PBFSolver::CalculateDeltaP(Particle *p){
    vec3 del_p(0.f);
    float h = 1.f;
    const unsigned int numberOfParticles = p->neighbors.size();
    for (unsigned int j = 0; j < numberOfParticles; j++)
    {
            vec3 r = ((p->neighbors[j]->pos - p->pos));
            const vec3 gradC_j = -1.f / rDensity * CalculateGradSpiky(p->pos - p->neighbors.at(j)->pos, 1.f);
            del_p -= (p->lambda + p->neighbors[j]->lambda + 0.f*CalculateSCorr(r, h, 0.1f, 4)) * gradC_j;
    }
    return del_p;
}


float PBFSolver::ViscousLaplacian(const vec3 diff, float h)
{
    float rad = length(diff);
    return (45.f / M_PI) * (h - rad) / pow(h, 6);
}


void PBFSolver::CalculateViscosityForce(Particle& p, float del_t)
{
    float KERNEL_RAD = 1.001f;
    vec3 force(0.f);
    for(Particle* n : p.neighbors)
    {
        vec3 tmp = (n->speed - p.speed) * p.mass * ViscousLaplacian(p.pos - n->pos, KERNEL_RAD);
        tmp *= (MU_VISCOSITY / p.density) * p.mass;
        force += tmp;

    }
    p.speed += force * del_t;
    p.pos_star = p.pos + p.speed * del_t;
}

