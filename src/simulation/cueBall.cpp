#include "simulation/cueBall.h"

namespace simulation {

// Constructor
CueBall::CueBall() {}
CueBall::CueBall(glm::vec3 _position) { position = _position; }

//==========================================
//  getter
//==========================================

float CueBall::getMass() const { return mass; }
float CueBall::getRadius() const { return radius; }
glm::vec3 CueBall::getPosition() const { return position; }
glm::vec3 CueBall::getVelocity() const { return velocity; }
glm::vec3 CueBall::getAcceleration() const { return force / mass; }
glm::vec3 CueBall::getForce() const { return force; }

//==========================================
//  setter
//==========================================
void CueBall::setMass(const float _mass) { mass = _mass; }

void CueBall::setRadius(const float _radius) { mass = _radius; }

void CueBall::setPosition(const glm::vec3 &_position) { position = _position; }

void CueBall::setVelocity(const glm::vec3 &_velocity) { velocity = _velocity; }

void CueBall::setAcceleration(const glm::vec3 &_acceleration) { force = _acceleration * mass; }

void CueBall::setForce(const glm::vec3 &_force) { force = _force; }

void CueBall::setForceField(const glm::vec3 &_forceField) { this->setAcceleration(_forceField); }

//==========================================
//  method
//==========================================

void CueBall::addPosition(const glm::vec3 &_position) { position += _position; }

void CueBall::addVelocity(const glm::vec3 &_velocity) { velocity += _velocity; }

void CueBall::addAcceleration(const glm::vec3 &_acceleration) { force += _acceleration * mass; }

void CueBall::addForce(const glm::vec3 &_force) { force += _force; }

void CueBall::resetCueBall(const glm::vec3 &_position) {
  this->setPosition(_position);
  this->setForce(glm::vec3(0.0f, 0.0f, 0.0f));
  this->setVelocity(glm::vec3(0.0f, 0.0f, 0.0f));
}
}  // namespace simulation
