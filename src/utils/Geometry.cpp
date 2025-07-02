#include "utils/Geometry.hpp"

#include <cmath>
#include <vector>

#include "SniperKernel/SniperLog.h"

bool trajectoryHitCylinder(double radius, double height, const vec3& pos, const vec3& dir, double& d1, double& d2) {
    std::vector<double> hitpoints;
    std::vector<int> label;

    vec3 hit1, hit2;
    if (trajectoryHitCylinderEdge(radius, pos, dir, d1, d2)) {
        hit1 = pos + d1 * dir;
        if (std::abs(hit1.z) < height / 2.0) { // Cylinder along z-axis
            hitpoints.push_back(d1);
            label.push_back(1);
        }
        hit2 = pos + d2 * dir;
        if (std::abs(hit2.z) < height / 2.0) {
            hitpoints.push_back(d2);
            label.push_back(2);
        }
    }

    if (trajectoryHitDisk(radius, height / 2.0, pos, dir, d1)) {
        hitpoints.push_back(d1);
        label.push_back(3);
    }

    if (trajectoryHitDisk(radius, -height / 2.0, pos, dir, d2)) {
        hitpoints.push_back(d2);
        label.push_back(4);
    }

    if (hitpoints.size() > 2) {
        LogWarn << "Number of intersections in the cylinder is more than two\n";
        for (std::size_t i = 0; i < hitpoints.size(); ++i) {
            LogWarn << "Hit length: " << hitpoints[i] << ", label: " << label[i] << '\n';
        }
        LogWarn << "Position: " << pos.x << ' ' << pos.y << ' ' << pos.z << " (" << theta(pos) << ", " << phi(pos) << ")\n";
        LogWarn << "Direction: " << dir.x << ' ' << dir.y << ' ' << dir.z << " (" << theta(dir) << ", " << phi(dir) << ")\n";
        return false;
    }

    else if (hitpoints.size() == 1) {
        d1 = hitpoints[0];
        d2 = hitpoints[0];
        return true;
    }

    else if (hitpoints.size() == 2) {
        d1 = std::min(hitpoints[0], hitpoints[1]);
        d2 = std::max(hitpoints[0], hitpoints[1]);
        return true;
    }

    else {
        return false;
    }

}

bool trajectoryHitCylinderEdge(double radius, const vec3& pos, const vec3& dir, double& d1, double& d2) {
    vec3 cross_pos = cross(vec3(0.0, 0.0, 1.0), pos);
    vec3 cross_dir = cross(vec3(0.0, 0.0, 1.0), dir);
    
    double a = mag2(cross_dir);
    if (a < 1.0e-5) return false;
    double b = 2.0 * dot(cross_pos, cross_dir);
    double c = mag2(cross_pos) - radius * radius;
    double delta = b * b - 4 * a * c;

    if (delta < 0.0) return false;
    d1 = (-b - std::sqrt(delta)) / (2.0 * a);
    d2 = (-b + std::sqrt(delta)) / (2.0 * a); 
    return true;
}

bool trajectoryHitDisk(double radius, double height, const vec3& pos, const vec3& dir, double& d) {
    if (!trajectoryHitPlane(height, pos, dir, d)) return false;
    vec3 hit = pos + d * dir;
    return (mag2(hit - vec3(0.0, 0.0, height)) <= radius * radius);
}

bool trajectoryHitPlane(double height, const vec3& pos, const vec3& dir, double& d) {
    if (std::abs(dir.z) < 1.0e-5) return false; // Trajectory is parallel to the plane (or is contained in it) and plane normal is along the z-axis
    d = (height - pos.z) / dir.z;
    return true;
}

bool trajectoryHitSphere(double radius, const vec3& pos, const vec3& dir, double& d1, double& d2) {
    double a = mag2(dir);
    double b = 2.0 * dot(pos, dir);
    double c = mag2(pos) - radius * radius;
    double delta = b * b - 4 * a * c;

    if (delta < 0.0) return false;
    d1 = (-b - std::sqrt(delta)) / (2.0 * a);
    d2 = (-b + std::sqrt(delta)) / (2.0 * a);
    return true;
}