#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

namespace simulation {
class mPlane {
 private:
  glm::vec3 position = glm::zero<glm::vec3>();
  glm::vec3 normal = glm::vec3(0.0f, 1.0f, 0.0f);
  float width = 5.0f;
  float height = 5.0f;

 public:
  mPlane();
  mPlane(const glm::vec3 &_position, const glm::vec3 &_normal, float _width, float _height);

  glm::vec3 getPosition() const;
  glm::vec3 getNormal() const;
  float getWidth() const;
  float getHeight() const;

  void setPosition(const glm::vec3 &_position);
  void setNormal(const glm::vec3 &_normal);
  void setWidth(float _width);
  void setHeight(float _height);

  glm::vec3 projectPointOntoPlane(const glm::vec3 &point) const;
};
} // namespace simulation