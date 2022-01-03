#include <algorithm>
#include <cassert>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <iostream>

#include <GLFW/glfw3.h>
#define GLAD_GL_IMPLEMENTATION
#include <glad/gl.h>
#undef GLAD_GL_IMPLEMENTATION
#include <glm/glm.hpp>
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG
#include <stb_image.h>
#undef STB_IMAGE_IMPLEMENTATION
#include "graphics.h"

// Unnamed namespace for global variables
namespace {
// Cameras
graphics::camera::Camera* currentCamera = nullptr;
// Control variables
bool isReset = false;
bool isShot = false;
bool isWindowSizeChanged = false;
bool isLightChanged = true;
bool mouseBinded = true;
int currentLight = 0;
int currentShader = 1;
int alignSize = 256;
float shootForce = 10000.0f;
// TODO (optional): Configs
// You should change line 32-35 if you add more shader / light / camera / mesh.
constexpr int LIGHT_COUNT = 3;
constexpr int CAMERA_COUNT = 1;
constexpr int CUE_BALL_COUNT = 16;
constexpr int PLANE_COUNT = 21;
constexpr int HOLE_COUNT = 6;
constexpr int MESH_COUNT = PLANE_COUNT + HOLE_COUNT + CUE_BALL_COUNT;
constexpr int SHADER_PROGRAM_COUNT = 3;
}  // namespace

int uboAlign(int i) { return ((i + 1 * (alignSize - 1)) / alignSize) * alignSize; }

void keyCallback(GLFWwindow* window, int key, int, int action, int) {
  // There are three actions: press, release, hold
  if (action != GLFW_PRESS) return;
  // Press ESC to close the window.
  if (key == GLFW_KEY_ESCAPE) {
    glfwSetWindowShouldClose(window, GLFW_TRUE);
    return;
  } else if (key == GLFW_KEY_F9) {
    // Disable / enable mouse cursor.
    if (mouseBinded)
      glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    else
      glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    mouseBinded = !mouseBinded;
  }
  switch (key) {
    // TODO: Detect key-events, to:
    //       1. switch among directional light, point light, and spot light, or
    //       2. switch between phong shader and gouraurd shader
    // Note: 1 key for 1 variable change
    case GLFW_KEY_H:
      currentLight = 0;
      isLightChanged = true;
      break;
    case GLFW_KEY_P:
      currentLight = 1;
      isLightChanged = true;
      break;
    case GLFW_KEY_O:
      currentLight = 2;
      isLightChanged = true;
      break;
    case GLFW_KEY_G: currentShader = 2; break;
    case GLFW_KEY_B: currentShader = 1; break;
    case GLFW_KEY_R: isReset = true; break;
    case GLFW_KEY_E: isShot = true; break;
    default: break;
  }
}

void resizeCallback(GLFWwindow* window, int width, int height) {
  OpenGLContext::framebufferResizeCallback(window, width, height);
  assert(currentCamera != nullptr);
  currentCamera->updateProjection(OpenGLContext::getAspectRatio());
  isWindowSizeChanged = true;
}

void resetCueBallPanel(GLFWwindow* window, simulation::Physics& physics) {
  ImGui::SetNextWindowSize(ImVec2(300.0f, 100.0f), ImGuiCond_Once);
  ImGui::SetNextWindowCollapsed(0, ImGuiCond_Once);
  ImGui::SetNextWindowPos(ImVec2(20.0f, 140.0f), ImGuiCond_Once);
  ImGui::SetNextWindowBgAlpha(0.2f);
  if (ImGui::Begin("Reset Cue Ball Panel")) {
    ImGui::Text("Cue ball x offset:");
    simulation::CueBall* cueBall = &(physics.cueBalls[0]);
    float* XOffset = cueBall->getPositionXOffsetPointer();
    ImGui::SliderFloat(" ", XOffset, -9.5f, 9.5f);
    float cueBallRadius = cueBall->getRadius();
    cueBall->resetCueBall(glm::vec3(*XOffset, cueBallRadius, 10.0f));
    if (ImGui::Button("OK")) {
      cueBall->setExist(true);
      mouseBinded = true;
      physics.isDead = false;
      glfwSetKeyCallback(window, keyCallback);
      glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
      *XOffset = 0.0f;
    }
  }

  ImGui::End();
}

void forceControlPanel() {
  ImGui::SetNextWindowSize(ImVec2(300.0f, 100.0f), ImGuiCond_Once);
  ImGui::SetNextWindowCollapsed(0, ImGuiCond_Once);
  ImGui::SetNextWindowPos(ImVec2(20.0f, 20.0f), ImGuiCond_Once);
  ImGui::SetNextWindowBgAlpha(0.2f);
  if (ImGui::Begin("Force Control Panel")) {
    ImGui::Text("Press F9 to enable/disable mouse cursor.");
    ImGui::Text("Shoot Force:");
    ImGui::SliderFloat(" ", &shootForce, 1000.0f, 20000.0f);
  }

  ImGui::End();
}

void renderUI(GLFWwindow* window, simulation::Physics& physics) {
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();
  if(physics.isDead) resetCueBallPanel(window, physics);
  forceControlPanel();
  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

int main() {
  // Initialize OpenGL context, details are wrapped in class.
  OpenGLContext::createContext(43, GLFW_OPENGL_CORE_PROFILE);
  GLFWwindow* window = OpenGLContext::getWindow();
  glfwSetWindowTitle(window, "Billiard Simulation");
  glfwSetKeyCallback(window, keyCallback);
  glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
  glfwSetFramebufferSizeCallback(window, resizeCallback);
#ifndef NDEBUG
  OpenGLContext::printSystemInfo();
  // This is useful if you want to debug your OpenGL API calls.
  OpenGLContext::enableDebugCallback();
#endif
  // Initialize dear-ImGui
  ImGui::CreateContext();
  ImGui::StyleColorsDark();
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init("#version 410 core");
  // Initialize shader
  std::vector<graphics::shader::ShaderProgram> shaderPrograms(SHADER_PROGRAM_COUNT);
  std::string filenames[SHADER_PROGRAM_COUNT] = {"shadow", "phong", "gouraud"};
  for (int i = 0; i < SHADER_PROGRAM_COUNT; ++i) {
    graphics::shader::VertexShader vs;
    graphics::shader::FragmentShader fs;
    vs.fromFile("../assets/shader/" + filenames[i] + ".vert");
    fs.fromFile("../assets/shader/" + filenames[i] + ".frag");
    shaderPrograms[i].attach(&vs, &fs);
    shaderPrograms[i].link();
    shaderPrograms[i].detach(&vs, &fs);
    shaderPrograms[i].use();
    // TODO: bind the uniform variables
    // Hint:
    //       1. you can set other uniforms you want in this for-loop
    //       2. check ShaderProgram class to know how to bind more easily
    //       3. It's safe to find and bind a non-exist uniform, it will just become NOP
    // Note:
    //       1. glUniform
    //        https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/glUniform.xhtml
    //       2. glGetUniformLocation
    //        https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/glGetUniformLocation.xhtml
    //       3. glUniformBlockBinding
    //        https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/glUniformBlockBinding.xhtml
    //       4. glGetUniformBlockIndex
    //        https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/glGetUniformBlockIndex.xhtml
    //       5. Check uniformBlockBinding and setUniform member function of ShaderProgram class

    shaderPrograms[i].uniformBlockBinding("model", 0);
    shaderPrograms[i].uniformBlockBinding("camera", 1);
    shaderPrograms[i].uniformBlockBinding("light", 2);
    shaderPrograms[i].setUniform("diffuseTexture", 0);
    shaderPrograms[i].setUniform("shadowMap", 1);
    shaderPrograms[i].setUniform("diffuseCubeTexture", 2);
    shaderPrograms[i].setUniform("isCube", 0);
  }
  graphics::buffer::UniformBuffer meshUBO, cameraUBO, lightUBO;
  // Calculate UBO alignment size
  glGetIntegerv(GL_UNIFORM_BUFFER_OFFSET_ALIGNMENT, &alignSize);
  constexpr int perMeshSize = 2 * sizeof(glm::mat4);
  constexpr int perCameraSize = sizeof(glm::mat4) + sizeof(glm::vec4);
  constexpr int perLightSize = sizeof(glm::mat4) + 2 * sizeof(glm::vec4);
  int perMeshOffset = uboAlign(perMeshSize);
  int perCameraOffset = uboAlign(perCameraSize);
  int perLightOffset = uboAlign(perLightSize);
  meshUBO.allocate(MESH_COUNT * perMeshOffset, GL_DYNAMIC_DRAW);
  cameraUBO.allocate(CAMERA_COUNT * perCameraOffset, GL_DYNAMIC_DRAW);
  lightUBO.allocate(LIGHT_COUNT * perLightOffset, GL_DYNAMIC_DRAW);
  // Default to first data
  meshUBO.bindUniformBlockIndex(0, 0, perMeshSize);
  cameraUBO.bindUniformBlockIndex(1, 0, perCameraSize);
  lightUBO.bindUniformBlockIndex(2, 0, perLightSize);
  // Get texture information
  int maxTextureSize = 1024;
  // Uncomment the following line if your GPU is very poor
  glGetIntegerv(GL_MAX_TEXTURE_SIZE, &maxTextureSize);
  maxTextureSize = std::min(maxTextureSize, 4096);
  // Camera
  std::vector<graphics::camera::CameraPTR> cameras;
  cameras.emplace_back(graphics::camera::QuaternionCamera::make_unique(glm::vec3(0, 3, 30)));
  assert(cameras.size() == CAMERA_COUNT);
  // TODO (Just an example for you, no need to modify here): Bind camera object's uniform buffer
  // Hint:
  //       1. what should we bind -> what will be used in shader: camera's view-projection matrix's & camera
  //          position's pointer
  //       2. where to bind -> remind VBO figure: we have to know the offset, size of the obj wanted to bind
  //       3. how to bind -> check spec slide to know binding procedure & trace the obj/class in the template to
  //          call class methods
  for (int i = 0; i < CAMERA_COUNT; ++i) {
    int offset = i * perCameraOffset;
    cameras[i]->initialize(OpenGLContext::getAspectRatio());
    cameraUBO.load(offset, sizeof(glm::mat4), cameras[i]->getViewProjectionMatrixPTR());
    cameraUBO.load(offset + sizeof(glm::mat4), sizeof(glm::vec4), cameras[i]->getPositionPTR());
  }
  currentCamera = cameras[0].get();
  // Lights
  glm::vec2 cutoff = glm::vec2(cos(glm::radians(5.0f)), glm::cos(glm::radians(22.0f)));
  std::vector<graphics::light::LightPTR> lights;
  lights.emplace_back(graphics::light::DirectionalLight::make_unique(glm::vec3(8, 6, 6)));
  lights.emplace_back(graphics::light::PointLight::make_unique(glm::vec3(8, 6, 6)));
  lights.emplace_back(graphics::light::Spotlight::make_unique(currentCamera->getFront(), cutoff));
  assert(lights.size() == LIGHT_COUNT);
  // TODO: Bind light object's buffer
  // Hint: look what we did when binding other UBO
  for (int i = 0; i < LIGHT_COUNT; ++i) {
    int offset = i * perLightOffset;
    lightUBO.load(offset, sizeof(glm::mat4), lights[i]->getLightSpaceMatrixPTR());
    lightUBO.load(offset + sizeof(glm::mat4), sizeof(glm::vec4), lights[i]->getLightVectorPTR());
    lightUBO.load(offset + sizeof(glm::mat4) + sizeof(glm::vec4), sizeof(glm::vec4),
                  lights[i]->getLightCoefficientsPTR());
  }
  // Texture
  graphics::texture::ShadowMap shadow(maxTextureSize);
  graphics::texture::Texture2D cueBallTextures[CUE_BALL_COUNT];
  graphics::texture::Texture2D colorOrange, colorBlack, transparent, wood, cloth;
  graphics::texture::TextureCubeMap dice;
  colorOrange.fromColor(glm::vec4(1, 0.5, 0, 1));
  colorBlack.fromColor(glm::vec4(0, 0, 0, 1));
  // TODO: Read texture(and set color) for objects respectively
  // Hint: check the calss of the variable(wood, colorOrange, dice) we've created for you
  //       fromFile member function
  wood.fromFile("../assets/texture/hardwood_floor.jpg");
  cloth.fromFile("../assets/texture/cloth.jpg");
  dice.fromFile("../assets/texture/posx.jpg", "../assets/texture/negx.jpg", "../assets/texture/posy.jpg",
                "../assets/texture/negy.jpg", "../assets/texture/posz.jpg", "../assets/texture/negz.jpg");
  
  for (int i = 0; i < CUE_BALL_COUNT; i++) {
    cueBallTextures[i].fromFile("../assets/texture/" + std::to_string(i) + ".jpeg");
  }

  /* ===== Generate cue balls ===== */
  simulation::Physics physics = simulation::Physics();
  for (int i = 0; i < CUE_BALL_COUNT; i++) {
    simulation::CueBall cueBall = simulation::CueBall();
    cueBall.setId(i);
    physics.cueBalls.push_back(cueBall);
  }
  physics.reset();
  /* ============================== */

  // Meshes
  std::vector<graphics::shape::ShapePTR> meshes;
  std::vector<graphics::texture::Texture*> diffuseTextures;
  {
    // Planes
    std::vector<GLfloat> vertexData;
    std::vector<GLuint> indexData;
    glm::mat4 model;
    for (int i = 0; i < PLANE_COUNT; i++) {
      float isWood = (i > 12) ? true : false;
      simulation::MPlane tableplane = physics.tablePlanes[i];
      float planeWidth = tableplane.getWidth();
      float planeHeight = tableplane.getHeight();
      glm::vec3 planePosition = tableplane.getPosition();
      glm::quat planeRotation = tableplane.getRotation();
      
      graphics::shape::Plane::generateVertices(vertexData, indexData, 40, planeWidth / 2, planeHeight / 2, isWood);
      auto plane = graphics::shape::Plane::make_unique(vertexData, indexData);
      model = glm::translate(glm::mat4(1), planePosition);
      model = glm::rotate(model, glm::angle(planeRotation), glm::axis(planeRotation));
      plane->setModelMatrix(model);
      meshes.emplace_back(std::move(plane));
      if (isWood) {
        diffuseTextures.emplace_back(&wood);
      } else {
        diffuseTextures.emplace_back(&cloth);
      }

      vertexData.clear();
      indexData.clear();
    }

    // Holes
    std::vector<glm::vec3> holePosition = {glm::vec3(-10.0f, -0.05f, -20.0f), glm::vec3(10.0f, -0.05f, -20.0f),
                                           glm::vec3(-10.45f, -0.5f, 0.0f),   glm::vec3(10.45f, -0.5f, 0.0f),
                                           glm::vec3(-10.0f, -0.05f, 20.0f),   glm::vec3(10.0f, -0.05f, 20.0f)};

    for (int i = 0; i < HOLE_COUNT; i++) {
      auto sphere = graphics::shape::Sphere::make_unique();
      model = glm::translate(glm::mat4(1), holePosition[i]);
      model = glm::scale(model, glm::vec3(1.0f));
      sphere->setModelMatrix(model);
      meshes.emplace_back(std::move(sphere));
      diffuseTextures.emplace_back(&colorBlack);
    }
    //graphics::shape::Plane::generateVertices(vertexData, indexData, 40, 20, 20, false);
    //auto sphere = graphics::shape::Sphere::make_unique();
    //auto cube = graphics::shape::Cube::make_unique();
    //auto ground = graphics::shape::Plane::make_unique(vertexData, indexData);

    /*glm::mat4 model = glm::translate(glm::mat4(1), glm::vec3(3, 0, 4));
    model = glm::scale(model, glm::vec3(2));
    model = glm::rotate(model, glm::half_pi<float>(), glm::vec3(1, 0, 0));
    sphere->setModelMatrix(model);*/

    /*model = glm::translate(glm::mat4(1), glm::vec3(-3, -1, 0));
    model = glm::scale(model, glm::vec3(2));
    cube->setModelMatrix(model);*/

    //glm::mat4 model = glm::translate(glm::mat4(1), glm::vec3(0, -3, 0));
    //ground->setModelMatrix(model);

    //meshes.emplace_back(std::move(ground));
    //diffuseTextures.emplace_back(&cloth);
    /*meshes.emplace_back(std::move(sphere));
    diffuseTextures.emplace_back(&colorOrange);
    meshes.emplace_back(std::move(cube));
    diffuseTextures.emplace_back(&wood);*/
  }

  /* ===== Do it after push cue balls to meshes =====
  assert(meshes.size() == MESH_COUNT);
  assert(diffuseTextures.size() == MESH_COUNT);
  

  // TODO: Bind light object's buffer
  // Hint: look what we did when binding other UBO
  for (int i = 0; i < MESH_COUNT; ++i) {
    int offset = i * perMeshOffset;
    meshUBO.load(offset, sizeof(glm::mat4), meshes[i]->getModelMatrixPTR());
    meshUBO.load(offset + sizeof(glm::mat4), sizeof(glm::mat4), meshes[i]->getNormalMatrixPTR());
  }
  */

  // This will not change in rendering loop
  shadow.bind(1);
  dice.bind(2);

  // Main rendering loop
  while (!glfwWindowShouldClose(window)) {
    // Polling events.
    glfwPollEvents();

    // Reset all cue balls
    if (isReset) {
      physics.reset();
      isReset = false;
    }

    // Update camera's uniforms if camera moves.
    bool isCameraMove = mouseBinded ? currentCamera->move(window) : false;
    if (isCameraMove || isWindowSizeChanged) {
      isWindowSizeChanged = false;
      cameraUBO.load(0, sizeof(glm::mat4), currentCamera->getViewProjectionMatrixPTR());
      cameraUBO.load(sizeof(glm::mat4), sizeof(glm::vec4), currentCamera->getPositionPTR());
      if (lights[currentLight]->getType() == graphics::light::LightType::Spot) {
        lights[currentLight]->update(currentCamera->getViewMatrix());
        int offset = currentLight * perLightOffset;
        glm::vec4 front = currentCamera->getFront();
        lightUBO.load(offset, sizeof(glm::mat4), lights[currentLight]->getLightSpaceMatrixPTR());
        lightUBO.load(offset + sizeof(glm::mat4), sizeof(glm::vec4), glm::value_ptr(front));
      }
    }
    // TODO: Switch light uniforms if light changes
    // Hint:
    //       1. we've load all the lights' unifroms eariler, so here we just tell shader where to start binding
    //       the next light info
    //       2. you should not bind the same light every time, because we are in a while-loop
    if (isLightChanged) {
      int offset = currentLight * perLightOffset;
      lightUBO.bindUniformBlockIndex(2, offset, perLightSize);
      if (lights[currentLight]->getType() == graphics::light::LightType::Spot) {
        lights[currentLight]->update(currentCamera->getViewMatrix());
        glm::vec4 front = currentCamera->getFront();
        lightUBO.load(offset, sizeof(glm::mat4), lights[currentLight]->getLightSpaceMatrixPTR());
        lightUBO.load(offset + sizeof(glm::mat4), sizeof(glm::vec4), glm::value_ptr(front));
      }
      isLightChanged = false;
    }
    // Render shadow first
    glViewport(0, 0, shadow.getSize(), shadow.getSize());
    glCullFace(GL_FRONT);
    shaderPrograms[0].use();
    shadow.bindFramebuffer();
    glClear(GL_DEPTH_BUFFER_BIT);


    /* ===== Do Physics Simulation here ===== */
    static float lastUpdateTime = 0.0f;
    float currentTime = glfwGetTime();
    float deltaTime = currentTime - lastUpdateTime;
    lastUpdateTime = currentTime;
    physics.setDeltaTime(deltaTime);
    physics.resolveCollision();
    physics.computeAllForce();
    // Shot
    if (isShot) {
      simulation::CueBall* cueBall = &physics.cueBalls[0];
      glm::vec3 shotDirection = cueBall->getPosition() - glm::vec3(currentCamera->getPosition());
      shotDirection.y = 0;
      shotDirection = glm::normalize(shotDirection);
      cueBall->addForce(shotDirection * shootForce);
      isShot = false;
    }
    physics.integrate();
    /* ====================================== */


    /* ===== Push cueBalls to meshes ===== */
    for (int i=0; i<CUE_BALL_COUNT; i++) {
      auto cueBall = graphics::shape::Sphere::make_unique();
      glm::vec3 cueBallPosition = physics.cueBalls[i].getPosition();
      glm::quat cueBallRotation = physics.cueBalls[i].getRotation();
      float cueBallRadius = physics.cueBalls[i].getRadius();
      glm::mat4 model = glm::translate(glm::mat4(1), cueBallPosition);

      // Not sure the relationship between radius and scale
      model = glm::scale(model, glm::vec3(cueBallRadius));

      // Maybe set rotation here
      //model = glm::rotate(model, glm::half_pi<float>(), glm::vec3(1, 0, 0));
      model *= glm::mat4_cast(cueBallRotation);

      cueBall->setModelMatrix(model);
      meshes.emplace_back(std::move(cueBall));
      diffuseTextures.emplace_back(&cueBallTextures[i]);
    }
    /* ============================== */

    assert(meshes.size() == MESH_COUNT);
    assert(diffuseTextures.size() == MESH_COUNT);

    // TODO: Bind light object's buffer
    // Hint: look what we did when binding other UBO
    for (int i = 0; i < MESH_COUNT; ++i) {
      int offset = i * perMeshOffset;
      meshUBO.load(offset, sizeof(glm::mat4), meshes[i]->getModelMatrixPTR());
      meshUBO.load(offset + sizeof(glm::mat4), sizeof(glm::mat4), meshes[i]->getNormalMatrixPTR());
    }

    for (int i = 0; i < MESH_COUNT; ++i) {
      // Update model matrix
      meshUBO.bindUniformBlockIndex(0, i * perMeshOffset, perMeshSize);
      meshes[i]->draw();
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glCullFace(GL_BACK);
    glViewport(0, 0, OpenGLContext::getWidth(), OpenGLContext::getHeight());
    // GL_XXX_BIT can simply "OR" together to use.
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // Render all objects
    shaderPrograms[currentShader].use();
    for (int i = 0; i < MESH_COUNT; ++i) {
      if (meshes[i]->getType() == graphics::shape::ShapeType::Cube) {
        shaderPrograms[currentShader].setUniform("isCube", 1);
      } else {
        shaderPrograms[currentShader].setUniform("isCube", 0);
      }
      meshUBO.bindUniformBlockIndex(0, i * perMeshOffset, perMeshSize);
      diffuseTextures[i]->bind(0);
      meshes[i]->draw();
    }

    /* ===== Pop all cueBalls ===== */
    for (int i = 0; i < CUE_BALL_COUNT; i++) {
      meshes.pop_back();
      diffuseTextures.pop_back();
    }
    /* ============================ */
    
    renderUI(window, physics);

#ifdef __APPLE__
    // Some platform need explicit glFlush
    glFlush();
#endif
    glfwSwapBuffers(window);
  }
  return 0;
}

