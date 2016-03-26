//
//  geom.hpp
//  Thanda

#ifndef geom_hpp
#define geom_hpp
#include "GL/glew.h"
#include <vector>
#include "../fluidSolver/fluidSolver.hpp"
class Geometry{
public:
    Geometry();
    std::vector<GLfloat> g_cube_vertex_buffer_data;
    virtual bool collisionDetect(Particle *p, float dt, glm::vec3 &coll_Pos, glm::vec3 &coll_Nor);
};

#endif /* geom_hpp */
