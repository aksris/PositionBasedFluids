//
//  camera.cpp
//  Thanda
//

#include "camera.hpp"

Camera::Camera(){
    initialFOV = 45.f;
    position = glm::vec3(0, 0, 5.5f / tan(glm::radians(initialFOV / 2.f)));
    horizontalAngle = 3.14f;
    verticalAngle = 0.f;
    speed = 20.f;
    mouseSpeed = 0.0005f;
}
glm::mat4 Camera::getViewMatrix(){
    return ViewMatrix;
}
glm::mat4 Camera::getProjectionMatrix(){
    return ProjectionMatrix;
}
void Camera::computeMatricesFromInputs(GLFWwindow* window){

    // glfwGetTime is called only once, the first time this function is called
    static double lastTime = glfwGetTime();

    // Compute time difference between current and last frame
    double currentTime = glfwGetTime();
    float deltaTime = float(currentTime - lastTime);

    // Get mouse position
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    // Reset mouse position for next frame
    glfwSetCursorPos(window, 1024/2, 768/2);

    // Compute new orientation
    horizontalAngle += mouseSpeed * float(1024/2 - xpos );
    verticalAngle   += mouseSpeed * float( 768/2 - ypos );

    // Direction : Spherical coordinates to Cartesian coordinates conversion
    glm::vec3 direction(
        cos(verticalAngle) * sin(horizontalAngle),
        sin(verticalAngle),
        cos(verticalAngle) * cos(horizontalAngle)
    );

    // Right vector
    glm::vec3 right = glm::vec3(
        sin(horizontalAngle - 3.14f/2.0f),
        0,
        cos(horizontalAngle - 3.14f/2.0f)
    );

    // Up vector
    glm::vec3 up = glm::cross( right, direction );

    // Move forward
    if (glfwGetKey( window, GLFW_KEY_W ) == GLFW_PRESS){
        position += direction * deltaTime * speed;
    }
    // Move backward
    if (glfwGetKey( window, GLFW_KEY_S ) == GLFW_PRESS){
        position -= direction * deltaTime * speed;
    }
    // Strafe right
    if (glfwGetKey( window, GLFW_KEY_D ) == GLFW_PRESS){
        position += right * deltaTime * speed;
    }
    // Strafe left
    if (glfwGetKey( window, GLFW_KEY_A ) == GLFW_PRESS){
        position -= right * deltaTime * speed;
    }
    //Up
    if (glfwGetKey( window, GLFW_KEY_Q ) == GLFW_PRESS){
        position -= up * deltaTime * speed;
    }
    //Down
    if (glfwGetKey( window, GLFW_KEY_E ) == GLFW_PRESS){
        position -= -up * deltaTime * speed;
    }

    float FoV = initialFOV;// - 5 * glfwGetMouseWheel();

    // Projection matrix : 45 Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
    ProjectionMatrix = glm::perspective(FoV, 4.0f / 3.0f, 0.1f, 400.0f);
    // Camera matrix
    ViewMatrix       = glm::lookAt(
                                position,           // Camera is here
                                position+direction, // and looks here : at the same position, plus "direction"
                                up                  // Head is up (set to 0,-1,0 to look upside-down)
                           );

    // For the next frame, the "last time" will be "now"
    lastTime = currentTime;
}
