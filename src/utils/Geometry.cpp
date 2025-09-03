#include "utils/Geometry.hpp"

#include <cmath>
#include <vector>

#include "SniperKernel/SniperLog.h"

std::pair<double, double> computeWpHits(const vec3& ipos, const vec3& dir) {
    double d_entry = 0.0;
    double d_exit = 0.0;

    double a = dir.x * dir.x + dir.y * dir.y; // = mag2(cross(vec3(0, 0, 1), dir)) = | z x dir |^2
    if (a < 1.0e-5) { // parallel to cylinder z-axis ==> hits only the end-caps
        d_entry = (dir.z > 0.0) ? (-22000.0 - ipos.z) / dir.z : (22000.0 - ipos.z) / dir.z;
        d_exit = (dir.z > 0.0) ? (22000.0 - ipos.z) / dir.z : (-22000.0 - ipos.z) / dir.z;
        return {d_entry, d_exit};
    }
    double b = 2.0 * (ipos.x * dir.x + ipos.y * dir.y); // = 2 * dot(cross(vec3(0, 0, 1), ipos), cross(vec3(0, 0, 1), dir)) = 2 * | z x ipos | * | z x dir |
    double c = ipos.x * ipos.x + ipos.y * ipos.y - 21750.0 * 21750.0; // = mag2(cross(vec3(0, 0, 1), ipos)) - R^2 = | z x ipos |^2 - R^2
    double delta = b * b - 4 * a * c;
    double sqrt_delta = std::sqrt(delta);

    d_entry = (-b - std::sqrt(delta)) / (2.0 * a);
    d_exit = (-b + std::sqrt(delta)) / (2.0 * a); 

    double z_entry = ipos.z + d_entry * dir.z;
    double z_exit = ipos.z + d_exit * dir.z;

    if (z_entry < -22000.0 || 22000.0 < z_entry) {
        d_entry = (dir.z > 0.0) ? (-22000.0 - ipos.z) / dir.z : (22000.0 - ipos.z) / dir.z;
    }
    if (z_exit < -22000.0 || 22000.0 < z_exit) {
        d_exit = (dir.z > 0.0) ? (22000.0 - ipos.z) / dir.z : (-22000.0 - ipos.z) / dir.z;
    }

    return {d_entry, d_exit};
}