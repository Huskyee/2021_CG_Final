#pragma once
#include <vector>
#include "cueBall.h"

namespace simulation {
class Physics {
 public:
  int cueBallCount;
  float deltaTime;
  glm::vec3 gravity;
  std::vector<CueBall> cueBalls;

  Physics();
  void computeAllForce();
  void computeCueBallForce(CueBall& cueBall);
  void integrate();
  void reset();
};
} // namespace simulation