#ifndef CUBE_H
#define CUBE_H
#include "geom.hpp"

class Cube: public Geometry
{
public:
    Cube();
    std::vector<GLfloat> g_cube_vertex_buffer_data;
    std::vector<GLuint> cub_idx;
    void push_vert_data();
    void push_indices_data();
    glm::mat4 modelMatrix;
    virtual bool collisionDetect(Particle *p, float dt, glm::vec3 &coll_Pos, glm::vec3 &coll_Nor);
};
// helpers for cube intersection
bool inUnitCube(glm::vec3 point);
float rayPlaneISX(glm::vec3 pos, glm::vec3 dir, glm::vec3 planePos, glm::vec3 norm);
bool nearlyEqual(float a, float b, float epsilon);
void checkSlab(glm::vec3 pos, glm::vec3 dir, float &nearMax, glm::vec3 &nearNorm, float &farMin, glm::vec3 &farNorm, glm::vec3 norm, glm::vec3 pos1, glm::vec3 pos2);
#endif // CUBE_H
