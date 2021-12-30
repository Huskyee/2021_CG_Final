#include "simulation/mPlane.h"

namespace simulation {

mPlane::mPlane() {}

mPlane::mPlane(const glm::vec3& _position, const glm::vec3& _normal, float _width, float _height) {
  position = _position;
  normal = _normal;
  width = _width;
  height = _height;
}


glm::vec3 mPlane::getPosition() const { return position; }

glm::vec3 mPlane::getNormal() const { return normal; }

float mPlane::getWidth() const { return width; }

float mPlane::getHeight() const { return height; }


void mPlane::setPosition(const glm::vec3& _position) { position = _position; }

void mPlane::setNormal(const glm::vec3& _normal) { normal = _normal; }

void mPlane::setWidth(float _width) { width = _width; }

void mPlane::setHeight(float _height) { height = _height; }

// currently ignore plane's width and height.
glm::vec3 mPlane::projectPointOntoPlane(const glm::vec3& point) const {
  glm::vec3 centerToPointVec = point - position;
  glm::vec3 onPlanePart = centerToPointVec - glm::dot(centerToPointVec, normal) * normal;
  return onPlanePart + position;
}
} // namespace simulation