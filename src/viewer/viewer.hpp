//
//  viewer.hpp
//  Thanda

#ifndef viewer_hpp
#define viewer_hpp
#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

//#include "particles.h"
//#include "common/shader.hpp"
//#include "common/controls.hpp"
//#include "common/vboindexer.hpp"
//#include "common/texture.hpp"
#include "../scene/scene.hpp"
#include <stdio.h>
#include <string>
#include <vector>
#include "../fluidSolver/fluidSolver.hpp"
#include "../fluidSolver/pbfsolver.h"
#include "../camera/camera.hpp"

#include <fstream>
#include <algorithm>

using namespace std;
#include "../geom/cube.h"
#include <stdlib.h>
#include <string.h>
class Viewer {
public:
    Viewer(int width, int height, Scene s);
    void initializeGL();
    void drawParticles(int ParticlesCount);
    void allocateParticleBuffers(GLuint billboard_vertex_buffer, GLuint particles_position_buffer, GLuint particles_color_buffer);
    void display();

    int w, h;
    Scene scene;
    Cube cube;
    GLuint CameraRight_worldspace_ID, CameraUp_worldspace_ID,  ViewProjMatrixID, programID, TextureID, programIDGeometry, geomMatrixID;
    GLFWwindow* window;
};
#endif /* viewer_hpp */
