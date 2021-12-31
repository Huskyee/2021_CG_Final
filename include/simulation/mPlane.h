#pragma once
#include <glm/glm.hpp>
#include <glm/ext/quaternion_float.hpp>
#include <glm/ext/quaternion_trigonometric.hpp>

namespace simulation {

class MPlane {
 private:
  glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f);
  glm::quat rotation = glm::angleAxis(0.0f, glm::vec3(0.0f, 1.0f, 0.0f));
  float width = 5.0f;
  float height = 5.0f;

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