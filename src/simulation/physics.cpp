#include "simulation/physics.h"

glm::vec3 position_list[16] = {
    glm::vec3(0.0f, 3.0f, 0.0f),
    glm::vec3(1.0f, 3.0f, 0.0f),
    glm::vec3(2.0f, 3.0f, 0.0f), 
    glm::vec3(3.0f, 3.0f, 0.0f), 
    glm::vec3(4.0f, 3.0f, 0.0f), 
    glm::vec3(5.0f, 3.0f, 0.0f), 
    glm::vec3(6.0f, 3.0f, 0.0f), 
    glm::vec3(7.0f, 3.0f, 0.0f),
    glm::vec3(8.0f, 3.0f, 0.0f), 
    glm::vec3(9.0f, 3.0f, 0.0f), 
    glm::vec3(10.0f, 3.0f, 0.0f), 
    glm::vec3(11.0f, 3.0f, 0.0f),
    glm::vec3(12.0f, 3.0f, 0.0f), 
    glm::vec3(13.0f, 3.0f, 0.0f), 
    glm::vec3(14.0f, 3.0f, 0.0f), 
    glm::vec3(15.0f, 3.0f, 0.0f),
};

namespace simulation {

Physics::Physics()
	: cueBallCount(16),
      deltaTime(0.01f),
      gravity(glm::vec3(0.0f, -9.8f, 0.0f)) {}

void Physics::computeAllForce() {
  for (int i = 0; i < cueBallCount; i++) {
    computeCueBallForce(cueBalls[i]);
  }
}

void Physics::computeCueBallForce(CueBall& cueBall) {
    cueBall.setForceField(gravity);
}

void Physics::integrate() {
  for (int i = 0; i < cueBallCount; i++) {
    CueBall* cueBall = &cueBalls[i];
    glm::vec3 force = cueBall->getForce();
    float mass = cueBall->getMass();
    cueBall->setAcceleration(force / mass);
    glm::vec3 acceleration = cueBall->getAcceleration();
    cueBall->addVelocity(acceleration * deltaTime);
    glm::vec3 velocity = cueBall->getVelocity();
    cueBall->addPosition(velocity * deltaTime);
  }
}

void Physics::reset() {
  for (int i = 0; i < cueBallCount; i++) {
    CueBall* cueBall = &cueBalls[i];
    cueBall->resetCueBall(position_list[i]);
  }
}

} // namespace simulation