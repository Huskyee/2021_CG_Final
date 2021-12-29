#pragma once
#include <glm/glm.hpp>

namespace simulation {
class CueBall {
 private:
  float mass = 5.0f;
  float radius = 0.5f;
  glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f);
  glm::vec3 velocity = glm::vec3(0.0f, 0.0f, 0.0f);
  glm::vec3 force = glm::vec3(0.0f, 0.0f, 0.0f);

 public:
  // Constructor
  CueBall();
  CueBall(glm::vec3 _position);

  //==========================================
  //  getter
  //==========================================

  float getMass() const;
  float getRadius() const;
  glm::vec3 getPosition() const;
  glm::vec3 getVelocity() const;
  glm::vec3 getAcceleration() const;
  glm::vec3 getForce() const;

  //==========================================
  //  setter
  //==========================================

  void setMass(const float _mass);
  void setRadius(const float _radius);
  void setPosition(const glm::vec3 &_position);
  void setVelocity(const glm::vec3 &_velocity);
  void setAcceleration(const glm::vec3 &_acceleration);
  void setForce(const glm::vec3 &_force);
  void setForceField(const glm::vec3 &_forcefield);

  //==========================================
  //  method
  //==========================================

  void addPosition(const glm::vec3 &_position);
  void addVelocity(const glm::vec3 &_velocity);
  void addAcceleration(const glm::vec3 &_acceleration);
  void addForce(const glm::vec3 &_force);
  void resetCueBall(const glm::vec3 &_position);
};
}  // namespace simulation