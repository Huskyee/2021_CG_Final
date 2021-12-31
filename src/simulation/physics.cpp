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
  //float c = 0.7f;
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
    glm::vec3 planeNormal = plane.getNormal();
    // use center of ball to obtain contact point and detect collision, not accurate
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
    }
}

void Physics::integrate() {
  for (int i = 0; i < cueBallCount; i++) {
    CueBall* cueBall = &cueBalls[i];
    glm::vec3 force = cueBall->getForce();
    float mass = cueBall->getMass();
    cueBall->setAcceleration(force / mass);
    glm::vec3 acceleration = cueBall->getAcceleration();
    cueBall->addVelocity(acceleration * deltaTime);

    float d = 0.001f;
    if (glm::length(cueBall->getVelocity()) < d)
        cueBall->setVelocity(glm::vec3(0.0f, 0.0f, 0.0f));

    glm::vec3 velocity = cueBall->getVelocity();
    cueBall->addPosition(velocity * deltaTime);
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