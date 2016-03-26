//
//  fluidSolver.hpp
//  Thanda

#ifndef fluidSolver_hpp
#define fluidSolver_hpp
#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <algorithm>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
#include "macgriddata.h"
#include <iostream>
using namespace glm;

enum geomtype {AIR = 0, FLUID = 1, SOLID = 2};
class Particle{

public:
    Particle();
    int gridIdx;
    glm::vec3 pos, speed;
    unsigned char r,g,b,a; // Color
    float size, angle, mass, density;
    float life; // Remaining life of the particle. if <0 : dead and unused.
    float cameradistance; // *Squared* distance to the camera. if dead : -1.0f


    bool operator<(const Particle& that) const {
        // Sort in reverse order : far particles drawn first.
        return this->cameradistance > that.cameradistance;
    }
};

class MACGrid{
public:
    MACGrid();
    void initialize();
    MACGrid& operator=(const MACGrid& val);
    MACGridDataX vel_U;
    MACGridDataY vel_V;
    MACGridDataZ vel_W;
    MACGridData P;
    MACGridData mNeighbors;
protected:
};

class FluidSolver{
public:
    FluidSolver();
    glm::vec3 containerBounds;

    int LastUsedParticle; int MaxParticles;
    std::vector<Particle> ParticlesContainer;
    std::vector<vec3> ParticlePos;

    std::vector<float> lambda;
    std::vector<vec3> del_p;

    std::vector<Particle> particle_save; //to save particle velocity
    MACGrid Grid;
    vec3 integratePos(const vec3 pos, const vec3 speed, float time_step, bool RK2);

    int findUnusedParticles();
    void sortParticles();
    void particlesInit();
    void genParticles(float particle_separation, float boundx, float boundy, float boundz);

};
#endif /* fluidSolver_hpp */