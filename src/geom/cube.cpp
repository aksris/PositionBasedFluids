#include "cube.h"

Cube::Cube()
{
    this->push_vert_data();
    this->push_indices_data();
    modelMatrix = glm::mat4();
}

void Cube::push_vert_data(){
    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);

    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);

    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);

    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);

    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);

    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);

    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);

    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(0.5f);
    g_cube_vertex_buffer_data.push_back(-0.5f);
}

void Cube::push_indices_data(){
    cub_idx.push_back(0);
    cub_idx.push_back(1);
    cub_idx.push_back(1);
    cub_idx.push_back(2);
    cub_idx.push_back(2);
    cub_idx.push_back(3);
    cub_idx.push_back(3);
    cub_idx.push_back(0);

    cub_idx.push_back(0);
    cub_idx.push_back(4);
    cub_idx.push_back(1);
    cub_idx.push_back(5);
    cub_idx.push_back(2);
    cub_idx.push_back(6);
    cub_idx.push_back(3);
    cub_idx.push_back(7);

    cub_idx.push_back(4);
    cub_idx.push_back(5);
    cub_idx.push_back(5);
    cub_idx.push_back(6);
    cub_idx.push_back(6);
    cub_idx.push_back(7);
    cub_idx.push_back(7);
    cub_idx.push_back(4);

}

bool Cube::collisionDetect(Particle *p, float dt, glm::vec3 &coll_Pos, glm::vec3 &coll_Nor) {
    // The collision detection algorithm here is a modification of the
    // slab-based collision checking method, which was originally
    // designed for things entering a cube.
    // We modify the algorithm to detect if a particle is leaving
    // the cube by moving the particle forward in time
    // and "rewinding" to see if it collides with the cube from outside to inside.

    // unit cube local params
    float _radius = 0.5f;
    glm::vec3 _center = glm::vec3(0.0f, 0.0f, 0.0f);

    glm::mat4 TF = modelMatrix; // local -> world transformation
    glm::mat4 inverseTF = glm::inverse(TF); // world -> local transformation

    // transform the particle's data into local coordinate system
    glm::vec4 worldVel = glm::vec4(p->speed, 0.0f);
    glm::vec4 worldPos = glm::vec4(p->pos, 1.0f) + worldVel * dt;

    glm::vec4 localPos = inverseTF * worldPos;
    glm::vec4 localVel = inverseTF * worldVel;

    glm::vec3 localPosVec3 = glm::vec3(localPos.x, localPos.y, localPos.z);
    glm::vec3 localVelVec3 = glm::vec3(localVel.x, localVel.y, localVel.z);

    // check if the particle stays inside the cube through the timestep
    if (inUnitCube(localPosVec3))
        return false;

    // if there is an intersection, determine the position and normal.
    glm::vec3 boxMax = glm::vec3(_radius, _radius, _radius);
    glm::vec3 boxMin = glm::vec3(-_radius, -_radius, -_radius);

    glm::vec3 xPlaneVector = glm::vec3(1.0f, 0.0f, 0.0f);
    glm::vec3 xPlane1Point = glm::vec3(boxMax.x, 0.0f, 0.0f);
    glm::vec3 xPlane2Point = glm::vec3(boxMin.x, 0.0f, 0.0f);
    glm::vec3 yPlaneVector = glm::vec3(0.0f, 1.0f, 0.0f);
    glm::vec3 yPlane1Point = glm::vec3(0.0f, boxMax.y, 0.0f);
    glm::vec3 yPlane2Point = glm::vec3(0.0f, boxMin.y, 0.0f);
    glm::vec3 zPlaneVector = glm::vec3(0.0f, 0.0f, 1.0f);
    glm::vec3 zPlane1Point = glm::vec3(0.0f, 0.0f, boxMax.z);
    glm::vec3 zPlane2Point = glm::vec3(0.0f, 0.0f, boxMin.z);

    float nearMax = -HUGE_VAL;
    float farMin = HUGE_VAL;
    glm::vec3 nearNorm;
    glm::vec3 farNorm;

    glm::vec3 dir = glm::normalize(-localVelVec3); // direction the particle is being "rewound" in
    // we want to cast a ray in the direction pointing towards the bounding volume, so may necessitate "correction"
    glm::vec3 dirBack = glm::normalize(-localPosVec3); // the direction going back to the center of the box
    if (glm::dot(dirBack, dir) <= 0.0f) dir = dirBack;

    if (dir.x != 0.0f)
        checkSlab(localPosVec3, dir, nearMax, nearNorm, farMin, farNorm, xPlaneVector, xPlane1Point, xPlane2Point);

    if (dir.y != 0.0f)
        checkSlab(localPosVec3, dir, nearMax, nearNorm, farMin, farNorm, yPlaneVector, yPlane1Point, yPlane2Point);

    if (dir.z != 0.0f)
        checkSlab(localPosVec3, dir, nearMax, nearNorm, farMin, farNorm, zPlaneVector, zPlane1Point, zPlane2Point);

    if (farMin < nearMax) // no intersection
    {
        return false;
    };

    glm::vec3 isxPoint = localPosVec3 + dir * nearMax;

    // get position
    glm::vec4 worldPoint = glm::vec4(isxPoint, 1.0f);
    worldPoint = TF * worldPoint;

    // get normal
    coll_Nor = -isxPoint; // we're treating the volume as only having inward normals. clamp this to a face normal.
    float xAbs = std::abs(coll_Nor.x);
    float yAbs = std::abs(coll_Nor.y);
    float zAbs = std::abs(coll_Nor.z);

    if (xAbs > yAbs && xAbs > zAbs) coll_Nor = glm::vec3(coll_Nor.x, 0.0f, 0.0f);
    else if (yAbs > xAbs && yAbs > zAbs) coll_Nor = glm::vec3(0.0f, coll_Nor.y, 0.0f);
    else if (zAbs > xAbs && zAbs > yAbs) coll_Nor = glm::vec3(0.0f, 0.0f, coll_Nor.z);
    coll_Nor = glm::normalize(coll_Nor);

    coll_Pos = glm::vec3(worldPoint.x, worldPoint.y, worldPoint.z);
    return true;
}

bool inUnitCube(glm::vec3 point) {
    // unit cube local params
    float _radius = 0.5f;
    glm::vec3 _center = glm::vec3(0.0f, 0.0f, 0.0f);
    return (
        point.x >= -_radius + _center.x &&
        point.x <= _radius + _center.x &&
        point.y >= -_radius + _center.y &&
        point.y <= _radius + _center.y &&
        point.z >= -_radius + _center.z &&
        point.z <= _radius + _center.z);
}

float rayPlaneISX(glm::vec3 pos, glm::vec3 dir, glm::vec3 planePos, glm::vec3 norm)
{
    float tIntersect = -1.0f;
    //http://www.cs.princeton.edu/courses/archive/fall00/cs426/lectures/raycast/sld017.htm
    // equation of a plane is dot(P, N) + d = 0
    // equation of the ray is point = P0 + t * v
    // t = -(P0 * N + d) / (V * N)
    glm::vec3 P = planePos;
    glm::vec3 N = norm;
    glm::vec3 P0 = pos;
    glm::vec3 V = dir;

    // watch out for rays parallel to the plane
    if (nearlyEqual(glm::dot(V, N), 0.0f, 0.00001f)) return tIntersect;

    float d = -glm::dot(P, N);

    tIntersect = -(glm::dot(P0, N) + d) / glm::dot(V, N);
    return tIntersect;
}

bool nearlyEqual(float a, float b, float epsilon)
{
    {
        const float absA = fabs(a);
        const float absB = fabs(b);
        const float diff = fabs(a - b);
        if (a == b) { // shortcut
            return true;
        }
        else if (a * b == 0) { // a or b or both are zero
            // relative error is not meaningful here
            return diff < (epsilon * epsilon);
        }
        else { // use relative error
            return diff / (absA + absB) < epsilon;
        }
    }
}

void checkSlab(glm::vec3 pos, glm::vec3 dir, float &nearMax, glm::vec3 &nearNorm, float &farMin, glm::vec3 &farNorm, glm::vec3 norm, glm::vec3 pos1, glm::vec3 pos2) {

    float nearReturn = rayPlaneISX(pos, dir, pos1, norm);

    float farReturn = rayPlaneISX(pos, dir, pos2, norm);

    if (farReturn < nearReturn)
    {
        float distSwap = farReturn;
        farReturn = nearReturn;
        nearReturn = distSwap;
    }

    if (nearReturn > nearMax) {
        nearMax = nearReturn;
        nearNorm = norm;
    }
    if (farReturn < farMin) {
        farMin = farReturn;
        farNorm = norm;
    }
}
