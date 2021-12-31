#pragma once
#include <glm/glm.hpp>
#include <glm/ext/quaternion_float.hpp>
#include <glm/ext/quaternion_trigonometric.hpp>

namespace simulation {

class MPlane {
 private:
  glm::vec3 position;
  glm::quat rotation;
  float width;
  float height;

 public:
  MPlane();
  MPlane(const glm::vec3 &_position, const glm::quat &_rotation, float _width, float _height);

  glm::vec3 getPosition() const;
  glm::vec3 getNormal() const;
  glm::quat getRotation() const;
  float getWidth() const;
  float getHeight() const;

  void setPosition(const glm::vec3 &_position);
  void setRotation(const glm::quat &_rotation);
  void setWidth(float _width);
  void setHeight(float _height);

  glm::vec3 projectPointOntoPlane(const glm::vec3 &_point) const;
  bool isPointOnPlane(const glm::vec3 &_point) const;
};

} // namespace simulation