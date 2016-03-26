//
//  fluidSolver.cpp
//  Thanda


#include "fluidSolver.hpp"
Particle::Particle(){
    pos = glm::vec3(0.f, 0.f, 0.f);
    speed = glm::vec3(0.f, 0.f, 0.f);
    r = 0;
    g = 0;
    b = 220;
    a = 230;
    size = 0.1f;
    angle = 45.f;
    mass = 1.f;
    life = 1.f;
    cameradistance = 10.f;
    density = 1.f;
}

MACGrid::MACGrid(){
}

void MACGrid::initialize(){
    vel_U.MACGridDataInitialize();
    vel_V.MACGridDataInitialize();
    vel_W.MACGridDataInitialize();
    P.MACGridDataInitialize();
}

MACGrid& MACGrid::operator =(const MACGrid& val){
    if (&val == this)
    {
        return *this;
    }
    vel_U = val.vel_U;
    vel_V = val.vel_V;
    vel_W = val.vel_W;
    P = val.P;
    return *this;
}

FluidSolver::FluidSolver(){
    LastUsedParticle = 0;
    MaxParticles = 200000;
}

int FluidSolver::findUnusedParticles(){

    for(int i=LastUsedParticle; i<MaxParticles; i++){
        if (ParticlesContainer[i].life < 0){
            LastUsedParticle = i;
            return i;
        }
    }

    for(int i=0; i<LastUsedParticle; i++){
        if (ParticlesContainer[i].life < 0){
            LastUsedParticle = i;
            return i;
        }
    }

    return 0; // All particles are taken, override the first one
}

void FluidSolver::sortParticles(){
    std::sort(&ParticlesContainer[0], &ParticlesContainer[MaxParticles]);
}

void FluidSolver::particlesInit(){
    for(int i=0; i<MaxParticles; i++){
        //        ParticlesContainer[i].life = -1.0f;
        ParticlesContainer[i].cameradistance = -1.0f;
        ParticlesContainer[i].size = 0.1f;
    }
}

void FluidSolver::genParticles(float particle_separation, float boundx, float boundy, float boundz){
    Particle p;
    for(float i = 0; i < boundx; i+= particle_separation){
        for(float j = 0; j < boundy ; j+= particle_separation){
            for(float k = 0; k <boundz ; k+= particle_separation){
                p.pos = glm::vec3(i, j, k);
                ParticlesContainer.push_back(p);
            }
        }
    }
    ParticlePos.resize(ParticlesContainer.size());
    lambda.resize(ParticlesContainer.size());
    del_p.resize(ParticlesContainer.size());
}

vec3 FluidSolver::integratePos(const vec3 pos, const vec3 speed, float time_step, bool RK2){
    vec3 new_pos(0.f);
    if(RK2){
        //RK2 integration
        vec3 k1 = speed * (pos) * time_step / 2.f;
        vec3 k2 = speed * (pos + k1) * time_step;
        new_pos = pos + k2;
    }
    else{
        //RK4
        vec3 k1 = speed * (pos) * time_step / 2.f;
        vec3 k2 = speed * (pos + k1) * time_step;
        vec3 k3 = speed * (pos + k2) * time_step;
        vec3 k4 = speed * (pos + k3);
        new_pos = pos + (0.666666f) * (k1 + 2.f * k2 + 2.f * k3 + k4);
    }
    return new_pos;
}
