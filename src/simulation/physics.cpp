#include "simulation/physics.h"

#include "simulation/mPlane.h"

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

    // adjust cueBall's pos and velocity due to collided with ground.
    // TODO: split collision removal logic
    glm::vec3 groundNormal = glm::normalize(glm::vec3(0.0f, 1.0f, 0.0f));
    mPlane ground(glm::vec3(0.0f, -3.0f, 0.0f), groundNormal, 5.0f, 5.0f); // currently hard-coded ground

    glm::vec3 hitPos = ground.projectPointOntoPlane(cueBall.getPosition());
    glm::vec3 centerToHitPosVec = hitPos - cueBall.getPosition();
    float e = 0.01f;
    float minDistance = cueBall.getRadius() + e;
    if (glm::length(centerToHitPosVec) < minDistance) {  // hitted
      cueBall.addPosition(ground.getNormal() * (minDistance - glm::length(centerToHitPosVec)));
      glm::vec3 toGroundVelocity = ground.getNormal() * glm::dot(ground.getNormal(), cueBall.getVelocity());
      float k = -1.9f;
      cueBall.addForce(k * cueBall.getMass() * toGroundVelocity / deltaTime);
    }
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