//
//  viewer.cpp
//  Thanda

#include "viewer.hpp"

Viewer::Viewer(int width, int height, Scene s){
    w = width;
    h = height;
    scene = s;
}

void Viewer::initializeGL(){
    if( !glfwInit() )
    {
        fprintf( stderr, "Failed to initialize GLFW\n" );
        getchar();
//        return -1;
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_RESIZABLE,GL_TRUE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window = glfwCreateWindow( w, h, "3D Viewer", NULL, NULL);
    if( window == NULL ){
        fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
        getchar();
        glfwTerminate();
//        return -1;
    }
    glfwMakeContextCurrent(window);

    // Initialize GLEW
    glewExperimental = true; // Needed for core profile
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        getchar();
        glfwTerminate();
//        return -1;
    }

    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    // Hide the mouse and enable unlimited mouvement
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // Set the mouse at the center of the screen
    glfwPollEvents();
    glfwSetCursorPos(window, w/2, h/2);

    // Dark blue background
    glClearColor(0.4f, 0.4f, 0.4f, 0.0f);

    // Enable depth test
    glEnable(GL_DEPTH_TEST);
    // Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS);
    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);
    programID = scene.LoadShaders( "src/resources/Particle.vertexshader", "src/resources/Particle.fragmentshader" );
    /* for particles and billboards */
    // Vertex shader
    CameraRight_worldspace_ID  = glGetUniformLocation(programID, "CameraRight_worldspace");
    CameraUp_worldspace_ID  = glGetUniformLocation(programID, "CameraUp_worldspace");
    ViewProjMatrixID = glGetUniformLocation(programID, "VP");
    TextureID  = glGetUniformLocation(programID, "myTextureSampler");
    /* for cube */
    programIDGeometry = scene.LoadShaders( "src/resources/TransformVertexShader.vertexshader", "src/resources/ColorFragmentShader.fragmentshader" );
    // Get a handle for our "MVP" uniform
    /*end cube*/

    /* MVP matrices */
    geomMatrixID = glGetUniformLocation(programIDGeometry, "MVP");
}

void Viewer::drawParticles(int ParticlesCount){
    // These functions are specific to glDrawArrays*Instanced*.
    // The first parameter is the attribute buffer we're talking about.
    // The second parameter is the "rate at which generic vertex attributes advance when rendering multiple instances"
    // http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribDivisor.xml
    glVertexAttribDivisor(0, 0); // particles vertices : always reuse the same 4 vertices -> 0
    glVertexAttribDivisor(1, 1); // positions : one per quad (its center)                 -> 1
    glVertexAttribDivisor(2, 1); // color : one per quad                                  -> 1

    // This is equivalent to :
    // for(i in ParticlesCount) : glDrawArrays(GL_TRIANGLE_STRIP, 0, 4),
    // but faster.
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, ParticlesCount);
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
}

void Viewer::allocateParticleBuffers(GLuint billboard_vertex_buffer, GLuint particles_position_buffer, GLuint particles_color_buffer){
    // 1rst attribute buffer : vertices
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
    glVertexAttribPointer(
                 0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
                 3,                  // size
                 GL_FLOAT,           // type
                 GL_FALSE,           // normalized?
                 0,                  // stride
                 (void*)0            // array buffer offset
    );

      // 2nd attribute buffer : positions of particles' centers
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
    glVertexAttribPointer(
                 1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
                 4,                                // size : x + y + z + size => 4
                 GL_FLOAT,                         // type
                 GL_FALSE,                         // normalized?
                 0,                                // stride
                 (void*)0                          // array buffer offset
     );

     // 3rd attribute buffer : particles' colors
     glEnableVertexAttribArray(2);
     glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
     glVertexAttribPointer(
                 2,                                // attribute. No particular reason for 1, but must match the layout in the shader.
                 4,                                // size : r + g + b + a => 4
                 GL_UNSIGNED_BYTE,                 // type
                 GL_TRUE,                          // normalized?    *** YES, this means that the unsigned char[4] will be accessible with a vec4 (floats) in the shader ***
                 0,                                // stride
                 (void*)0                          // array buffer offset
     );
}

void Viewer::display(){
//    geom.g_cube_vertex_buffer_data
    initializeGL();
    GLuint elementbuffer;
    glGenBuffers(1, &elementbuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, cube.cub_idx.size() * sizeof(GLuint), &cube.cub_idx[0] , GL_STATIC_DRAW);

    GLuint cubevertexbuffer;
    glGenBuffers(1, &cubevertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, cubevertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, cube.g_cube_vertex_buffer_data.size() * sizeof(GLfloat), &cube.g_cube_vertex_buffer_data[0], GL_STATIC_DRAW);

    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    // Create solver instance.
    PBFSolver* fluid = new PBFSolver(scene.containerBoundsMin, scene.containerBoundsMax);

    static GLfloat* g_particule_position_size_data = new GLfloat[fluid->MaxParticles * 4];
    static GLubyte* g_particule_color_data         = new GLubyte[fluid->MaxParticles * 4];

    GLuint Texture = scene.loadDDS("src/resources/particle.DDS");
    // The VBO containing the 4 vertices of the particles.
    // Thanks to instancing, they will be shared by all particles.
    static const GLfloat g_vertex_buffer_data[] = {
         -0.5f, -0.5f, 0.0f,
          0.5f, -0.5f, 0.0f,
         -0.5f,  0.5f, 0.0f,
          0.5f,  0.5f, 0.0f,
    };
    GLuint billboard_vertex_buffer;
    glGenBuffers(1, &billboard_vertex_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

    GLuint particles_position_buffer;
    glGenBuffers(1, &particles_position_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
    // Initialize with empty (NULL) buffer : it will be updated later, each frame.
    glBufferData(GL_ARRAY_BUFFER, fluid->MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);

    GLuint particles_color_buffer;
    glGenBuffers(1, &particles_color_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
    // Initialize with empty (NULL) buffer : it will be updated later, each frame.
    glBufferData(GL_ARRAY_BUFFER, fluid->MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW);

    // Generate particles.
    fluid->genParticles(scene.particle_separation, scene.particleBounds.x, scene.particleBounds.y, scene.particleBounds.z);

    Camera camera;
    //construct mac grid
    int counter = 0;
    int csome = 0;
    double lastTime = glfwGetTime();
    do{
        // Clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        double currentTime = glfwGetTime();
        double delta = currentTime - lastTime;
        lastTime = currentTime;

        delta = 0.1f;
        positions.clear();
        camera.computeMatricesFromInputs(window);
        glm::mat4 ProjectionMatrixParticles = camera.getProjectionMatrix();
        glm::mat4 ViewMatrixParticles = camera.getViewMatrix();
        // We will need the camera's position in order to sort the particles
        // w.r.t the camera's distance.
        // There should be a getCameraPosition() function in common/controls.cpp,
        // but this works too.
        glm::vec3 CameraPosition(glm::inverse(ViewMatrixParticles)[3]);

        glm::mat4 ViewProjectionMatrixParticles = ProjectionMatrixParticles * ViewMatrixParticles;

        //apply gravity forces and update pos
        for (int i = 0; i < fluid->ParticlesContainer.size(); ++i){
            fluid->ApplyForces(fluid->ParticlesContainer[i], delta * 0.5f);
            //fluid->CalculateViscosityForce(fluid->ParticlesContainer[i], delta * 0.5f);
        }

        // Update simulation.
        fluid->uGrid.update(fluid->ParticlesContainer);
        for(int j = 0; j < fluid->ParticlesContainer.size(); ++j){
            fluid->FindNeighbors(&(fluid->ParticlesContainer[j]));
        }
        //simulation loop
        for(int i = 0; i < 6; i++){
            for(Particle &p : fluid->ParticlesContainer){
                fluid->CalculateLagrangeMultiplier(&p);
            }
            for(int j = 0; j < fluid->ParticlesContainer.size(); ++j){
                fluid->del_p[j] = fluid->CalculateDeltaP(&fluid->ParticlesContainer[j] );
                if (fluid->ParticlesContainer.at(j).pos_star.y < EPSILON){
                    fluid->ParticlesContainer.at(j).speed *= (vec3(1.f, -1.f, 1.f));
                    fluid->ParticlesContainer.at(j).pos_star += fluid->ParticlesContainer.at(j).speed
                    * (float)delta ;
                }
            }
            for(int i = 0; i < fluid->ParticlesContainer.size(); ++i){
                fluid->ParticlesContainer[i].pos_star += fluid->del_p[i];
            }
        }

        for(int i = 0; i < fluid->ParticlesContainer.size(); ++i){
            fluid->ParticlesContainer[i].speed = 1.f / (float)delta * (fluid->ParticlesContainer[i].pos_star - fluid->ParticlesContainer[i].pos);
            fluid->ParticlesContainer[i].pos = fluid->ParticlesContainer[i].pos_star;
        }

        // Bounds clamping
        for(Particle &p : fluid->ParticlesContainer){
            p.pos[0] = p.pos[0] < -fluid->uGrid.dimensions[0] / 2 ? -fluid->uGrid.dimensions[0] / 2 : p.pos[0];
                        p.pos[1] = p.pos[1] < -fluid->uGrid.dimensions[1] / 2 ? -fluid->uGrid.dimensions[1] / 2 : p.pos[1];
                        p.pos[2] = p.pos[2] < -fluid->uGrid.dimensions[2] / 2 ? -fluid->uGrid.dimensions[2] / 2 : p.pos[2];

                        p.pos[0] = p.pos[0] > fluid->uGrid.dimensions[0] / 2 ? fluid->uGrid.dimensions[0] / 2 - 1 : p.pos[0];
                        p.pos[1] = p.pos[1] > fluid->uGrid.dimensions[1] / 2 ? fluid->uGrid.dimensions[1] / 2 - 1 : p.pos[1];
                        p.pos[2] = p.pos[2] > fluid->uGrid.dimensions[2] / 2 ? fluid->uGrid.dimensions[2] / 2 - 1 : p.pos[2];

        }

        int ParticlesCount = 0;
        for(int i=0; i< fluid->ParticlesContainer.size(); i++){

            Particle& p = fluid->ParticlesContainer[i]; // shortcut

            p.cameradistance = glm::length2( p.pos - CameraPosition );

            // Fill the GPU buffer
            g_particule_position_size_data[4*ParticlesCount+0] = p.pos.x;
            g_particule_position_size_data[4*ParticlesCount+1] = p.pos.y;
            g_particule_position_size_data[4*ParticlesCount+2] = p.pos.z;

            g_particule_position_size_data[4*ParticlesCount+3] = p.size;

            g_particule_color_data[4*ParticlesCount+0] = p.r;
            g_particule_color_data[4*ParticlesCount+1] = p.g;
            g_particule_color_data[4*ParticlesCount+2] = p.b;
            g_particule_color_data[4*ParticlesCount+3] = p.a;
            positions.push_back(openvdb::Vec3f(p.pos.x,p.pos.y,p.pos.z));
            ParticlesCount++;

        }

        using namespace openvdb::tools;

        // Initialize the OpenVDB and OpenVDB Points library.  This must be called at least
        // once per program and may safely be called multiple times.
        openvdb::initialize();
        openvdb::points::initialize();


        // Create a linear transform with voxel size of 10.0
        const float voxelSize = 10.0f;
        openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform();

        // Create the PointDataGrid, position attribute is mandatory
        PointDataGrid::Ptr pointDataGrid = createPointDataGrid<PointDataGrid>(
                    positions, TypedAttributeArray<openvdb::Vec3f>::attributeType(), *transform);


#ifdef DEBUG
        // Output leaf nodes
        std::cout << "Leaf Nodes: " << pointDataGrid->tree().leafCount() << std::endl;

        // Output point count
        std::cout << "Point Count: " << pointCount(pointDataGrid->tree()) << std::endl;
#endif
        // Create a VDB file object.
        string file_name = "pbf_";
        string ext = ".vdb";
        openvdb::io::File file(file_name+std::to_string(counter)+ext);
        // Add the grid pointer to a container.
        openvdb::GridPtrVec grids;
        grids.push_back(pointDataGrid);
        // Write out the contents of the container.
//        file.write(grids);
        file.close();
        counter++;

       int MaxParticles = fluid->ParticlesContainer.size();
       glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
       glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
       glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLfloat) * 4, g_particule_position_size_data);

       glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
       glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
       glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLubyte) * 4, g_particule_color_data);

       glEnable(GL_BLEND);
       glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

       // Use our shader
       glUseProgram(programID);

       // Bind our texture in Texture Unit 0
       glActiveTexture(GL_TEXTURE0);
       glBindTexture(GL_TEXTURE_2D, Texture);
       // Set our "myTextureSampler" sampler to user Texture Unit 0
       glUniform1i(TextureID, 0);

       glUniform3f(CameraRight_worldspace_ID, ViewMatrixParticles[0][0], ViewMatrixParticles[1][0], ViewMatrixParticles[2][0]);
       glUniform3f(CameraUp_worldspace_ID   , ViewMatrixParticles[0][1], ViewMatrixParticles[1][1], ViewMatrixParticles[2][1]);

       glUniformMatrix4fv(ViewProjMatrixID, 1, GL_FALSE, &ViewProjectionMatrixParticles[0][0]);

        allocateParticleBuffers(billboard_vertex_buffer, particles_position_buffer, particles_color_buffer);

        drawParticles(ParticlesCount);
        //cube
        glDisable(GL_BLEND);
        glUseProgram(programIDGeometry);
        glm::mat4 cubeModelMatrix(glm::scale(glm::mat4(1.f),glm::vec3(scene.containerBounds.x, scene.containerBounds.y, scene.containerBounds.z )));
        cube.modelMatrix = cubeModelMatrix;
        glm::mat4 cubeMVP = ProjectionMatrixParticles * ViewMatrixParticles * cubeModelMatrix;

        glUniformMatrix4fv(geomMatrixID, 1, GL_FALSE, &cubeMVP[0][0]);
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, cubevertexbuffer);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0  );

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);

        glDrawElements(GL_LINES, 36, GL_UNSIGNED_INT, 0);
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);


        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();
        csome++;

    }while( (glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS ) &&
            glfwWindowShouldClose(window) == 0 && csome < 100);
    glfwTerminate();
}

