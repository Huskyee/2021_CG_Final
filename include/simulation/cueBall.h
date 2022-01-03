#pragma once
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/constants.hpp>

namespace simulation {

class CueBall {
 private:
  int id;
  float positionXOffset = 0.0f;
  bool isExist = true;

  float radius;
  float mass;
  glm::mat3 inertiaBody;
  glm::mat3 inertiaBodyInverse;

  glm::vec3 position = glm::zero<glm::vec3>();
  glm::quat rotation = glm::identity<glm::quat>();
  glm::vec3 linearMomentum = glm::zero<glm::vec3>();
  glm::vec3 angularMomentum = glm::zero<glm::vec3>();

  glm::vec3 force = glm::zero<glm::vec3>();
  glm::vec3 torque = glm::zero<glm::vec3>();

 public:
  // Constructor
  CueBall();
  CueBall(const glm::vec3 &_position,
          const glm::quat &_rotation = glm::identity<glm::quat>(),
          float _radius = 0.5f,
          float _mass = 5.0f);

  //==========================================
  //  getter
  //==========================================

  int getId() const;
  float *getPositionXOffsetPointer();
  bool getExist() const;

  float getRadius() const;
  float getMass() const;
  glm::mat3 getInertiaBody() const;
  glm::mat3 getInertiaBodyInverse() const;

  glm::vec3 getPosition() const;
  glm::quat getRotation() const;
  glm::vec3 getLinearMomentum() const;
  glm::vec3 getAngularMomentum() const;

  glm::mat3 getInertiaInverse() const;
  glm::mat3 getRotationMatrix() const;
  glm::vec3 getVelocity() const;
  glm::vec3 getAngularVelocity() const;
  glm::vec3 getAcceleration() const;

  glm::vec3 getForce() const;
  glm::vec3 getTorque() const;

  //==========================================
  //  setter
  //==========================================

  void setId(const int _id);
  void setExist(const bool _isExist);

  void setPosition(const glm::vec3 &_position);
  void setRotation(const glm::quat &_rotation);
  void setLinearMomentum(const glm::vec3 &_linearMomentum);
  void setAngularMomentum(const glm::vec3 &_angularMomentum);

  void setVelocity(const glm::vec3 &_velocity);
  void setAcceleration(const glm::vec3 &_acceleration);

  void setForce(const glm::vec3 &_force);
  void setForceField(const glm::vec3 &_forcefield);
  void setTorque(const glm::vec3 &_torque);

  //==========================================
  //  method
  //==========================================

  void addPosition(const glm::vec3 &_position);
  void addRotation(const glm::quat &_rotation);
  void addLinearMomentum(const glm::vec3 &_linearMomentum);
  void addAngularMomentum(const glm::vec3 &_angularMomentum);

  void addVelocity(const glm::vec3 &_velocity);
  void addAcceleration(const glm::vec3 &_acceleration);

  void addForce(const glm::vec3 &_force);
  void addTorque(const glm::vec3 &_torque);

  void resetCueBall(const glm::vec3 &_position);
};

}  // namespace simulation