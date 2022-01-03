#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

namespace simulation {
class CueBall {
 private:
  int id;
  float mass = 5.0f;
  float radius = 0.5f;
  float positionXOffset = 0.0f;
  bool isExist = true;
  glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f);
  glm::vec3 velocity = glm::vec3(0.0f, 0.0f, 0.0f);
  glm::vec3 force = glm::vec3(0.0f, 0.0f, 0.0f);
  glm::quat rotation = glm::angleAxis(0.0f, glm::vec3(1.0f, 0.0f, 0.0f));

 public:
  // Constructor
  CueBall();
  CueBall(glm::vec3 _position);

  //==========================================
  //  getter
  //==========================================

  int getId() const;
  float getMass() const;
  float getRadius() const;
  float *getPositionXOffsetPointer();
  bool getExist() const;
  glm::vec3 getPosition() const;
  glm::vec3 getVelocity() const;
  glm::vec3 getAcceleration() const;
  glm::vec3 getForce() const;
  glm::quat getRotation() const;

  //==========================================
  //  setter
  //==========================================

  void setId(const int _id);
  void setMass(const float _mass);
  void setRadius(const float _radius);
  void setExist(const bool _isExist);
  void setPosition(const glm::vec3 &_position);
  void setVelocity(const glm::vec3 &_velocity);
  void setAcceleration(const glm::vec3 &_acceleration);
  void setForce(const glm::vec3 &_force);
  void setRotation(const glm::quat &_rotation);
  void setForceField(const glm::vec3 &_forcefield);

  //==========================================
  //  method
  //==========================================

  void addPosition(const glm::vec3 &_position);
  void addVelocity(const glm::vec3 &_velocity);
  void addAcceleration(const glm::vec3 &_acceleration);
  void addForce(const glm::vec3 &_force);
  void addRotation(const glm::quat &_rotation);
  void resetCueBall(const glm::vec3 &_position);
};
}  // namespace simulation