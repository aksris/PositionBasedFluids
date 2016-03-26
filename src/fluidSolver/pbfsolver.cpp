#include "pbfsolver.h"

PBFSolver::PBFSolver()
{
    rDensity = 1000.f;
    numParticles = ParticlesContainer.size();
}

void PBFSolver::ApplyForces(Particle& p, float del_t){
    //gravity force
    p.speed += glm::vec3(0.f, -9.81f, 0.f) * del_t;
    p.pos += p.speed * del_t;
    int index = Grid.P.getCellIndex(p.pos.x, p.pos.y, p.pos.z);
    p.gridIdx = index;
}

void PBFSolver::initializeGrid(){
    for (Particle p : FluidSolver::ParticlesContainer){
        iNeighbors[p.gridIdx] = &p;
    }
}
float CalculatePoly6(float r, float h){
    return (315.f * pow((h*h - r*r),3) / (64.f * M_PI * pow(h, 9)));
}
void PBFSolver:: NeighborSearch(Particle *p, std::vector<Particle> &neighbors) {
    numParticles = ParticlesContainer.size();
    for (int i = 0; i < numParticles; i++) {
        float dist = glm::length(p->pos - ParticlesContainer[i].pos);
        if (dist < PBF_H) {
            neighbors.push_back(ParticlesContainer[i]);
        }
    }
}
float PBFSolver::CalculateDensity(Particle p){
    //poly6 kernel
    float h = 10.f; //cell width
    float r= 0.f; //distance between particles i and j
    float density = 0.f;

    std::vector<Particle> neighbors;
    NeighborSearch(&p, neighbors);

    //calculate r for all the neighbors and sum
    for(int i = 0; i < neighbors.size(); ++i){
        r = glm::distance(p.pos, neighbors[i].pos);
        density += neighbors.at(i).mass * CalculatePoly6(r, PBF_H);
    }
    //    std::cout << "Density " << density << std::endl;
    return density;
}

float PBFSolver::CalculateGradientConstraint(Particle p, Particle k){
    float h = 10.f;
    std::vector<Particle> neighbors;
    NeighborSearch(&p, neighbors);
    if(k.pos == p.pos){
        float sum = 0.f;


        for(int i = 0; i < neighbors.size(); ++i){
            float r = glm::distance(p.pos, neighbors[i].pos);
            //            std::cout << "distance " << r << std::endl;
            sum += CalculateSpiky(r, PBF_H);
        }
        return sum / rDensity;
    }
    else{
        float r = glm::distance(p.pos, k.pos);
        float val = -CalculateSpiky(r, PBF_H);
        return val / rDensity;
    }
    return h;
}

float PBFSolver::CalculateLambda(Particle p){
    float numerator;
    float denominator;
    std::vector<Particle> neighbors;
    NeighborSearch(&p, neighbors);
    numerator = (CalculateDensity(p) / rDensity) - 1;
    float sum =0.f;
    for(int i = 0; i < neighbors.size(); ++i){
        sum += pow(CalculateGradientConstraint(p, neighbors[i]), 2);
    }
    denominator = sum + EPSILON;
    return numerator / denominator;
}

float PBFSolver::CalculateSCorr(Particle p){
    float sum_numerator = 0.f, h= 10.f;
    std::vector<Particle> neighbors;
    NeighborSearch(&p, neighbors);
    for(int i = 0; i < neighbors.size(); ++i){
        float r = glm::distance(p.pos, neighbors[i].pos);
        sum_numerator += 15.f * pow((PBF_H - r), 3) / (M_PI * pow(PBF_H, 6));
    }

    float denominator = 15.f * pow((PBF_H - del_q * PBF_H), 3) / (M_PI * pow(PBF_H, 6));
    return pow((sum_numerator/denominator), 4);
}

float CalculateSpiky(float r, float h){
    return -45.f * pow((h - r), 2) / (M_PI * pow(h, 6));
}

vec3 PBFSolver::CalculateDeltaP(Particle p, int i){
    vec3 sum(0.f), r(0.f);
    float h = 10.f;
    std::vector<Particle> neighbors;
    NeighborSearch(&p, neighbors);
    for(int j = 0; j < neighbors.size(); ++j){
        r = (p.pos - neighbors[j].pos) + p.size;
        sum += (lambda[i] + lambda[j]) * -45.f * powf((PBF_H - glm::length(normalize(r))), 2.f) / (M_PI * powf(PBF_H, 6.f));
    }
    //    std::cout << "sum " << sum.x << ", " << sum.y << ", " << sum.z << std::endl;
    return sum / rDensity;
}
