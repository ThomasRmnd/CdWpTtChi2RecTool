#ifndef CDWPTTCHI2RECTOOL_UTILS_GEOMETRY_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_GEOMETRY_HPP_

#include <utility>

#include "utils/vec3.hpp"

/**
 * @brief Calculate the entry and exit point of a track on the WP edges (cylinder + end caps).
 * 
 * This function assumes that the initial position is already inside the WP.
 * 
 * @param ipos The initial positions of the track.
 * @param dir The direction of the track.
 * 
 * @return A pair of distances respectively corresponding to the entry and exit point.
 * 
 * - Entry = ipos + pair.first * dir.
 * 
 * - Exit = ipos + pair.second * dir.
 */
std::pair<double, double> computeWpHits(const vec3& ipos, const vec3& dir);

#endif // CDWPTTCHI2RECTOOL_UTILS_GEOMETRY_HPP_