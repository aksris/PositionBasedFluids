//
//  camera.hpp
//  Thanda
//

#ifndef camera_hpp
#define camera_hpp

#include <GLFW/glfw3.h>
// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
using namespace glm;

class Camera{
public:
    Camera();
    glm::mat4 ViewMatrix;
    glm::mat4 ProjectionMatrix;
    glm::vec3 position;
    float horizontalAngle, verticalAngle, initialFOV, speed, mouseSpeed;
    void computeMatricesFromInputs(GLFWwindow *window);
    glm::mat4 getViewMatrix();
    glm::mat4 getProjectionMatrix();
};

#endif /* camera_hpp */
