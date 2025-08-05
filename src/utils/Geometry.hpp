#ifndef CDWPTTCHI2RECTOOL_UTILS_GEOMETRY_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_GEOMETRY_HPP_

#include <utility>

#include "utils/vec3.hpp"

std::pair<double, double> computeWpHits(const vec3& ipos, const vec3& dir);

bool trajectoryHitCylinder(double radius, double height, const vec3& pos, const vec3& dir, double& d1, double& d2);
bool trajectoryHitCylinderEdge(double radius, const vec3& pos, const vec3& dir, double& d1, double& d2);
bool trajectoryHitDisk(double radius, double height, const vec3& pos, const vec3& dir, double& d);
bool trajectoryHitPlane(double height, const vec3& pos, const vec3& dir, double& d);
bool trajectoryHitSphere(double radius, const vec3& pos, const vec3& dir, double& d1, double& d2);

#endif // CDWPTTCHI2RECTOOL_UTILS_GEOMETRY_HPP_