#include "simulation/mPlane.h"
#include <glm/gtc/quaternion.hpp>

namespace simulation {

MPlane::MPlane() {}

MPlane::MPlane(const glm::vec3& _position, const glm::quat& _rotation, float _width, float _height) {
  position = _position;
  rotation = _rotation;
  width = _width;
  height = _height;
}


glm::vec3 MPlane::getPosition() const { return position; }

glm::quat MPlane::getRotation() const { return rotation; }

glm::vec3 MPlane::getNormal() const { return glm::mat3_cast(rotation)[1]; }

float MPlane::getWidth() const { return width; }

float MPlane::getHeight() const { return height; }


void MPlane::setPosition(const glm::vec3& _position) { position = _position; }

void MPlane::setRotation(const glm::quat& _rotation) { rotation = _rotation; }

void MPlane::setWidth(float _width) { width = _width; }

void MPlane::setHeight(float _height) { height = _height; }


// this will ignore plane's width and height
glm::vec3 MPlane::projectPointOntoPlane(const glm::vec3& _point) const {
  glm::vec3 centerToPointVec = _point - position;
  glm::vec3 normal = getNormal();
  glm::vec3 onPlanePart = centerToPointVec - glm::dot(centerToPointVec, normal) * normal;
  return onPlanePart + position;
}

// this will take plane's width and height into account
bool MPlane::isPointOnPlane(const glm::vec3& _point) const {
	float e = 0.00001f;
	if (glm::length(_point - projectPointOntoPlane(_point)) >= e) return false;

	glm::vec3 centerToPointVec = _point - position;
	glm::vec3 pointLocalPos = glm::transpose(glm::mat3_cast(rotation)) * centerToPointVec;
	return glm::abs(pointLocalPos.x) <= width / 2
		&& glm::abs(pointLocalPos.z) <= height / 2;
}

} // namespace simulation