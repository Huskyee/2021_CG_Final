#include "simulation/physics.h"
#include "simulation/mPlane.h"


//glm::vec3 position_list[16] = {
//    glm::vec3(0.0f, 3.0f, 0.0f),
//    glm::vec3(1.0f, 3.0f, 0.0f),
//    glm::vec3(2.0f, 3.0f, 0.0f), 
//    glm::vec3(3.0f, 3.0f, 0.0f), 
//    glm::vec3(4.0f, 3.0f, 0.0f), 
//    glm::vec3(5.0f, 3.0f, 0.0f), 
//    glm::vec3(6.0f, 3.0f, 0.0f), 
//    glm::vec3(7.0f, 3.0f, 0.0f),
//    glm::vec3(8.0f, 3.0f, 0.0f), 
//    glm::vec3(9.0f, 3.0f, 0.0f), 
//    glm::vec3(10.0f, 3.0f, 0.0f), 
//    glm::vec3(11.0f, 3.0f, 0.0f),
//    glm::vec3(12.0f, 3.0f, 0.0f), 
//    glm::vec3(13.0f, 3.0f, 0.0f), 
//    glm::vec3(14.0f, 3.0f, 0.0f), 
//    glm::vec3(15.0f, 3.0f, 0.0f),
//};

float root_3 = pow(3, 0.5);
float ballRadius = 0.5f;

glm::vec3 cueBallPositionList[16] = {
    glm::vec3(0.0f, ballRadius, 10.0f),                                        // Cue ball
    glm::vec3(0.0f, ballRadius, -10.0f),                                       // 1
    glm::vec3(-ballRadius, ballRadius, -10.0f - root_3 * ballRadius),          // 2
    glm::vec3(ballRadius, ballRadius, -10.0f - root_3 * ballRadius),           // 3
    glm::vec3(-2 * ballRadius, ballRadius, -10.0f - 2 * root_3 * ballRadius),  // 4
    glm::vec3(0.0f, ballRadius, -10.0f - 2 * root_3 * ballRadius),             // 5
    glm::vec3(2 * ballRadius, ballRadius, -10.0f - 2 * root_3 * ballRadius),   // 6
    glm::vec3(-3 * ballRadius, ballRadius, -10.0f - 3 * root_3 * ballRadius),  // 7
    glm::vec3(-ballRadius, ballRadius, -10.0f - 3 * root_3 * ballRadius),      // 8
    glm::vec3(ballRadius, ballRadius, -10.0f - 3 * root_3 * ballRadius),       // 9
    glm::vec3(3 * ballRadius, ballRadius, -10.0f - 3 * root_3 * ballRadius),   // 10
    glm::vec3(-4 * ballRadius, ballRadius, -10.0f - 4 * root_3 * ballRadius),  // 11
    glm::vec3(-2 * ballRadius, ballRadius, -10.0f - 4 * root_3 * ballRadius),  // 12
    glm::vec3(0.0f, ballRadius, -10.0f - 4 * root_3 * ballRadius),             // 13
    glm::vec3(2 * ballRadius, ballRadius, -10.0f - 4 * root_3 * ballRadius),   // 14
    glm::vec3(4 * ballRadius, ballRadius, -10.0f - 4 * root_3 * ballRadius),   // 15
};

namespace simulation {

Physics::Physics()
	: cueBallCount(16),
      deltaTime(0.01f),
      gravity(glm::vec3(0.0f, -9.8f, 0.0f))
{
    //      y
    //      |
    //      |
    //      |_______ x
    //     /
    //    /
    //   z   
    
    // Bottom
    tablePlanes.emplace_back(glm::vec3(0.0f, 0.0f, 0.0f), glm::angleAxis(0.0f, glm::vec3(0.0f, 0.0f, 0.0f)), 20.0f,
                             40.0f);

    // Right (+x)
    tablePlanes.emplace_back(glm::vec3(10.0f, 0.5f, 0.0f),
                             glm::angleAxis(glm::half_pi<float>(), glm::vec3(0.0f, 0.0f, 1.0f)), 1.0f, 39.0f);

    // Left (-x)
    tablePlanes.emplace_back(glm::vec3(-10.0f, 0.5f, 0.0f),
                             glm::angleAxis(-glm::half_pi<float>(), glm::vec3(0.0f, 0.0f, 1.0f)), 1.0f, 39.0f);

    // Front (+z)
    tablePlanes.emplace_back(glm::vec3(0.0f, 0.5f, 20.0f),
                             glm::angleAxis(-glm::half_pi<float>(), glm::vec3(1.0f, 0.0f, 0.0f)), 19.0f, 1.0f);

    // Back (-z)
    tablePlanes.emplace_back(glm::vec3(0.0f, 0.5f, -20.0f),
                             glm::angleAxis(glm::half_pi<float>(), glm::vec3(1.0f, 0.0f, 0.0f)), 19.0f, 1.0f);
}

void Physics::computeAllForce() {
  for (int i = 0; i < cueBallCount; i++) {
    computeCueBallForce(cueBalls[i]);
  }

  float d = 0.002f;
  float k = 1.9f;
  float c = 0.9f;
  // ball and ball
  for (int i = 0; i < cueBallCount; i++) {
      for (int j = i + 1; j < cueBallCount; j++) {
           computeCueBallPairForce(cueBalls[i], cueBalls[j], d, c);
      }
  }
  // ball and table
  for (int i = 0; i < tablePlanes.size(); i++) {
      for (int j = 0; j < cueBallCount; j++) {
          computeCueBallTableForce(cueBalls[j], tablePlanes[i], d, k);
      }
  }
}

void Physics::computeCueBallForce(CueBall& cueBall) {
    cueBall.setForceField(gravity);
    cueBall.setTorque(glm::zero<glm::vec3>());
}

void Physics::computeCueBallPairForce(CueBall& cueBallA, CueBall& cueBallB, float d, float c) {
    glm::vec3 ABVec = cueBallB.getPosition() - cueBallA.getPosition();
    glm::vec3 ABDir = glm::normalize(ABVec);
    float ABLength = glm::length(ABVec);
    float minDistance = cueBallA.getRadius() + cueBallB.getRadius() + d;
    if (ABLength < minDistance) {  // collided

        // ==== apply force due to collision between the balls ====
        float AToBSpeed = glm::dot(ABDir,
            // obtain the newest velocity by considering the force that not applied yet.
            // may not be a good solution
            cueBallA.getVelocity() + cueBallA.getForce() / cueBallA.getMass() * deltaTime);
        float BToASpeed = glm::dot(-ABDir,
            cueBallB.getVelocity() + cueBallB.getForce() / cueBallB.getMass() * deltaTime);
        // at least one is moving toward the other, otherwise no need to apply force
        if (AToBSpeed > 0.0f || BToASpeed > 0.0f) {
            float AMass = cueBallA.getMass();
            float BMass = cueBallB.getMass();
            float ASpeed = AToBSpeed;   // respect to ABDir
            float BSpeed = -BToASpeed;  // respect to ABDir
            float newASpeed = c * (ASpeed * (AMass - BMass) + 2 * BMass * BSpeed) / (AMass + BMass);
            float newBSpeed = c * (BSpeed * (BMass - AMass) + 2 * AMass * ASpeed) / (AMass + BMass);
            cueBallA.addForce(AMass * (newASpeed - ASpeed) * ABDir / deltaTime);
            cueBallB.addForce(BMass * (newBSpeed - BSpeed) * ABDir / deltaTime);
        }
        // =========================================================

    }
}

void Physics::computeCueBallTableForce(CueBall& cueBall, const MPlane& plane, float d, float k) {
  /* glm::vec3 planeNormal = plane.getNormal();
  // use center of ball to obtain contact point and detect collision, not accurate
  // Fixed: use the closest point near the plane
  glm::vec3 contact = plane.projectPointOntoPlane(cueBall.getPosition());
  if (!plane.isPointOnPlane(contact)) return;
  glm::vec3 centerToContactVec = contact - cueBall.getPosition();
  float centerToContactLength = glm::length(centerToContactVec);
  float minDistance = cueBall.getRadius() + d;
  bool isBallStuckInPlane = glm::dot(-centerToContactVec, planeNormal) < 0.0f;
  if (isBallStuckInPlane || centerToContactLength < minDistance) {
      float toPlaneSpeed = glm::dot(-planeNormal, cueBall.getVelocity());
      if (toPlaneSpeed > 0.0f) {
          cueBall.addForce(k * cueBall.getMass() * toPlaneSpeed * planeNormal / deltaTime);
      }
  }*/

  //constexpr float eEPSILON = 0.01f;
  //constexpr float coefResist = 0.8f;
  //constexpr float coefFriction = 0.03f;

  //float cueBallRaduis = cueBall.getRadius();
  //glm::vec3 planeNormal = glm::normalize(plane.getNormal());
  //glm::vec3 planePosition = plane.getPosition();
  //glm::vec3 cueBallForce = cueBall.getForce();
  //glm::vec3 cueBallVelocity = cueBall.getVelocity();
  //glm::vec3 cueBallContactPoint = cueBall.getPosition() - cueBallRaduis * planeNormal;
  //// Decompose cue ball velocity vector into vn, vt
  //// vn: parallel to plane normal
  //// vt: perpendicular to plane normal
  //glm::vec3 vn = float(glm::dot(cueBallVelocity, -planeNormal) / pow(glm::length(-planeNormal), 2)) * (-planeNormal);
  //glm::vec3 vt = cueBallVelocity - vn;

  //bool closeToThePlane = glm::dot(planeNormal, cueBallContactPoint - planePosition) < eEPSILON;
  //bool headingIn = glm::dot(planeNormal, cueBallVelocity) < 0;
  //// Before: v = vn + vt
  //// After: v = -kr * vn + vt
  //if (closeToThePlane && headingIn) {
  //  cueBall.setVelocity(-coefResist * vn + vt);
  //}

  //// Contact force = -(N dot f) * N
  //// Friction = -kf(-N dot f) * vt
  //bool onThePlane = abs(glm::dot(planeNormal, cueBallContactPoint - planePosition)) < eEPSILON;
  //bool hasContactForce = glm::dot(planeNormal, cueBallForce) < 0;
  //bool movingAlongThePlane = abs(glm::dot(planeNormal, cueBallVelocity)) < eEPSILON;
  //if (onThePlane) {
  //  if (hasContactForce) {
  //    glm::vec3 contactForce = -glm::dot(planeNormal, cueBallForce) * planeNormal;
  //    cueBall.addForce(contactForce);
  //  }
  //  if (movingAlongThePlane) {
  //    glm::vec3 friction = coefFriction * glm::dot(planeNormal, cueBallForce) * vt;
  //    cueBall.addForce(friction);
  //  }
  //}
  
  constexpr float eEPSILON = 0.01f;
  constexpr float coefResist = 0.8f;
  constexpr float coefKineticFriction = 20.0f;
  constexpr float coefStaticFriction = coefKineticFriction * 2.0f;
  constexpr float thresholdSlidingSpeed = 0.01f;
  constexpr float coefSpinningSpeedDown = 1.0f;

  auto&& planeNormal = plane.getNormal();
  auto&& cueBallPosition = cueBall.getPosition();
  auto&& cueBallRadius = cueBall.getRadius();
  auto&& cueBallVelocity = cueBall.getVelocity();
  auto&& cueBallForce = cueBall.getForce();
  auto&& contactPoint = cueBallPosition - cueBallRadius * planeNormal;
  auto&& centerToContactPoint = contactPoint - cueBallPosition;
  auto&& contactForce = glm::dot(cueBallForce, planeNormal) * planeNormal;
  contactForce = glm::dot(contactForce, planeNormal) < 0.0f ? contactForce : glm::zero<glm::vec3>();

  auto&& vn = glm::dot(cueBallVelocity, planeNormal) * planeNormal;
  auto&& vt = cueBallVelocity - vn;

  auto&& cueBallAngularVelocity = cueBall.getAngularVelocity();
  auto&& contactPointVelocityAngularTerm = glm::cross(cueBallAngularVelocity, centerToContactPoint);
  auto&& contactPointOnPlaneVelocityAngularTerm =
      contactPointVelocityAngularTerm - glm::dot(contactPointVelocityAngularTerm, planeNormal) * planeNormal;
  auto&& contactPointOnPlaneVelocity = vt + contactPointOnPlaneVelocityAngularTerm;

  bool onPlane = glm::dot(contactPoint - plane.getPosition(), planeNormal) < eEPSILON;
  bool headingIn = glm::dot(cueBallVelocity, planeNormal) < 0.0f;
  bool hasContactForce = glm::length(contactForce) > 0.0f;
  bool sliding = glm::length(contactPointOnPlaneVelocity) > thresholdSlidingSpeed;
  bool rotatingAroundPlaneNormal = glm::dot(cueBallAngularVelocity, planeNormal) != 0.0f;

  if (onPlane) {
    if (headingIn) {
      cueBall.setVelocity(vt + -coefResist * vn);
    }
    if (hasContactForce) {
      cueBall.addForce(-contactForce);
    }
    if (sliding) { // sliding, use kinetic friction
      auto&& friction = -coefKineticFriction * contactPointOnPlaneVelocity;
      auto&& torque = glm::cross(centerToContactPoint, friction);
      cueBall.addForce(friction);
      cueBall.addTorque(torque);
    } else { // rolling or stationary, use static friction
      auto&& friction = -coefStaticFriction * vt; // vt is 0 in case of stationary, so no friction will be applied
      auto&& torque = glm::cross(centerToContactPoint, -friction);
      cueBall.addForce(friction);
      cueBall.addTorque(torque);
    }
    if (rotatingAroundPlaneNormal) {
      auto&& angularVelocityAlongPlaneNormal = glm::dot(cueBallAngularVelocity, planeNormal) * planeNormal;
      auto&& torque = -coefSpinningSpeedDown * angularVelocityAlongPlaneNormal;
      cueBall.addTorque(torque);
    }
  }
}

void Physics::integrate() {
  for (auto&& cueBall : cueBalls) {
    cueBall.addLinearMomentum(cueBall.getForce() * deltaTime);
    cueBall.addAngularMomentum(cueBall.getTorque() * deltaTime);
    cueBall.addPosition(cueBall.getVelocity() * deltaTime);
    auto&& angularVelocity = cueBall.getAngularVelocity();
    if (glm::length(angularVelocity) > 0.0f)
        cueBall.addRotation(glm::angleAxis(
            glm::length(angularVelocity) * deltaTime,
            glm::normalize(angularVelocity)
        ));
  }
}

void Physics::reset() {
  for (int i = 0; i < cueBallCount; i++) {
    CueBall* cueBall = &cueBalls[i];
    cueBall->resetCueBall(cueBallPositionList[i]);
  }
}

void Physics::resolveCollision() {
    int iteration = 10;
    // this d should be smaller than the one used when applying force due to collision
    float d = 0.001f;

    for (int k = 0; k < iteration; k++) {
        // ball and ball
        for (int i = 0; i < cueBallCount; i++) {
            for (int j = i + 1; j < cueBallCount; j++) {
                resolveCollision(cueBalls[i], cueBalls[j], d);
            }
        }
        // ball and table
        for (int i = 0; i < tablePlanes.size(); i++) {
            for (int j = 0; j < cueBallCount; j++) {
                resolveCollision(cueBalls[j], tablePlanes[i], d);
            }
        }
    }
}

void Physics::resolveCollision(CueBall& cueBallA, CueBall& cueBallB, float d) {
    float massA = cueBallA.getMass();
    float massB = cueBallB.getMass();
    float massTotal = massA + massB;
    glm::vec3 ABVec = cueBallB.getPosition() - cueBallA.getPosition();
    glm::vec3 ABDir = glm::normalize(ABVec);
    float ABLength = glm::length(ABVec);
    float minDistance = cueBallA.getRadius() + cueBallB.getRadius() + d;
    if (ABLength < minDistance) {
        float moveDistance = minDistance - ABLength;
        // may use a better distribution
        cueBallA.addPosition(-ABDir * moveDistance * (massB / massTotal));
        cueBallB.addPosition(ABDir * moveDistance * (massA / massTotal));
    }
}

void Physics::resolveCollision(CueBall& cueBall, const MPlane& plane, float d) {
    glm::vec3 planeNormal = plane.getNormal();
    // use center of ball to obtain contact point and detect collision, not accurate
    glm::vec3 contact = plane.projectPointOntoPlane(cueBall.getPosition());
    if (!plane.isPointOnPlane(contact)) return;
    glm::vec3 centerToContactVec = contact - cueBall.getPosition();
    float centerToContactLength = glm::length(centerToContactVec);
    float minDistance = cueBall.getRadius() + d;
    bool isBallStuckInPlane = glm::dot(-centerToContactVec, planeNormal) < 0.0f;
    if (isBallStuckInPlane || centerToContactLength < minDistance) {
        // may cause weird behavior in case the ball is stuck in plane
        // TODO: fix it
        float moveDistance = isBallStuckInPlane
            ? centerToContactLength + minDistance : minDistance - centerToContactLength;
        cueBall.addPosition(planeNormal * moveDistance);
    }
}

} // namespace simulation