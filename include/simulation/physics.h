#pragma once
#include <vector>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/constants.hpp>
#include <math.h>
#include <iostream>
#include "cueBall.h"
#include "mPlane.h"

namespace simulation {
class Physics {
 public:
  int cueBallCount;
  float deltaTime;
  bool isDead;
  glm::vec3 gravity;
  std::vector<CueBall> cueBalls;
  std::vector<MPlane> tablePlanes;

  Physics();
  void computeAllForce();
  void computeCueBallForce(CueBall& cueBall);
  void computeCueBallPairForce(CueBall &cueBallA, CueBall &cueBallB, float d = 0.002f, float c = 0.9f);
  void computeCueBallTableForce(CueBall &cueBall, const MPlane &plane, float d = 0.002f, float k = 1.9f);
  void integrate();
  void reset();
  void resolveCollision();
  void resolveCollision(CueBall &cueBallA, CueBall &cueBallB, float d = 0.001f);
  void resolveCollision(CueBall &cueBall, const MPlane &plane, float d = 0.001f);
  void holeDetection(CueBall *cueBall);
};
} // namespace simulation