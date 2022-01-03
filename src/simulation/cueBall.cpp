#include "simulation/cueBall.h"

namespace simulation {

// Constructor
CueBall::CueBall() : CueBall(glm::zero<glm::vec3>()) {}
CueBall::CueBall(const glm::vec3 &_position, const glm::quat &_rotation, float _radius, float _mass) {
  position = _position;
  rotation = _rotation;
  radius = _radius;
  mass = _mass;
  inertiaBody = glm::mat3(0.4f * _mass * _radius * _radius);
  inertiaBodyInverse = glm::mat3(1.0f / (0.4f * _mass * _radius * _radius));
}

//==========================================
//  getter
//==========================================

float CueBall::getRadius() const { return radius; }
float CueBall::getMass() const { return mass; }
glm::mat3 CueBall::getInertiaBody() const { return inertiaBody; }
glm::mat3 CueBall::getInertiaBodyInverse() const { return inertiaBodyInverse; }

glm::vec3 CueBall::getPosition() const { return position; }
glm::quat CueBall::getRotation() const { return rotation; }
glm::vec3 CueBall::getLinearMomentum() const { return linearMomentum; }
glm::vec3 CueBall::getAngularMomentum() const { return angularMomentum; }

glm::mat3 CueBall::getInertiaInverse() const {
  auto &&rotateMat = getRotationMatrix();
  return rotateMat * inertiaBodyInverse * glm::transpose(rotateMat);
}
glm::mat3 CueBall::getRotationMatrix() const { return glm::mat3_cast(rotation); }
glm::vec3 CueBall::getVelocity() const { return linearMomentum / mass; }
glm::vec3 CueBall::getAngularVelocity() const { return getInertiaInverse() * angularMomentum; }
glm::vec3 CueBall::getAcceleration() const { return force / mass; }

glm::vec3 CueBall::getForce() const { return force; }
glm::vec3 CueBall::getTorque() const { return torque; }

//==========================================
//  setter
//==========================================

void CueBall::setPosition(const glm::vec3 &_position) { position = _position; }
void CueBall::setRotation(const glm::quat &_rotation) { rotation = _rotation; }
void CueBall::setLinearMomentum(const glm::vec3 &_linearMomentum) { linearMomentum = _linearMomentum; }
void CueBall::setAngularMomentum(const glm::vec3 &_angularMomentum) { angularMomentum = _angularMomentum; }

void CueBall::setVelocity(const glm::vec3 &_velocity) { linearMomentum = mass * _velocity; }
void CueBall::setAcceleration(const glm::vec3 &_acceleration) { force = _acceleration * mass; }

void CueBall::setForce(const glm::vec3 &_force) { force = _force; }
void CueBall::setForceField(const glm::vec3 &_forceField) { this->setAcceleration(_forceField); }
void CueBall::setTorque(const glm::vec3 &_torque) { torque = _torque; }

//==========================================
//  method
//==========================================

void CueBall::addPosition(const glm::vec3 &_position) { position += _position; }
void CueBall::addRotation(const glm::quat &_rotation) { rotation *= _rotation; }
void CueBall::addLinearMomentum(const glm::vec3 &_linearMomentum) { linearMomentum += _linearMomentum; }
void CueBall::addAngularMomentum(const glm::vec3 &_angularMomentum) { angularMomentum += _angularMomentum; }

void CueBall::addVelocity(const glm::vec3 &_velocity) { linearMomentum += mass * _velocity; }
void CueBall::addAcceleration(const glm::vec3 &_acceleration) { force += _acceleration * mass; }

void CueBall::addForce(const glm::vec3 &_force) { force += _force; }
void CueBall::addTorque(const glm::vec3 &_torque) { torque += _torque; }

void CueBall::resetCueBall(const glm::vec3 &_position) {
  this->setPosition(_position);
  this->setRotation(glm::identity<glm::quat>());
  this->setLinearMomentum(glm::zero<glm::vec3>());
  this->setAngularMomentum(glm::zero<glm::vec3>());
  this->setForce(glm::zero<glm::vec3>());
  this->setTorque(glm::zero<glm::vec3>());
}

}  // namespace simulation
