/**
 * Copyright (c) 2015 Carnegie Mellon University, Daniel Maturana <dimatura@cmu.edu>
 *
 * For License information please see the LICENSE file in the root directory.
 *
 */

#ifndef RAY_HPP_TVFS3AKQ
#define RAY_HPP_TVFS3AKQ

#include <tuple>

#include <Eigen/Core>
#include <Eigen/Dense>

namespace ca { namespace scrollgrid {

template<typename Scalar>
class Ray2 {
 public:
  typedef Eigen::Matrix<Scalar, 2, 1> Vec2;

 public:
  /**
   * @param origin: origin of ray
   * @param direction: *unit* direction of ray
   *
   */
  Ray2(const Vec2& origin,
       const Vec2& direction) :
         origin(origin),
         direction(direction.normalized()),
         tmin(Scalar(0)),
         tmax(std::numeric_limits<Scalar>::max()),
         invdir(Scalar(1.) / direction.array()),
         sign(invdir.x() < 0,
              invdir.y() < 0)
  {
  }

 public:
  Vec2 point_at(Scalar t) const {
    return origin + t*direction;
  }

 public:
  Vec2 origin, direction;
  Scalar tmin, tmax; /// ray min and max distances
  Vec2 invdir; // for convenience in AABB intersection
  std::tuple<int, int> sign;

};

template<typename Scalar>
class Ray3 {
 public:
  typedef Eigen::Matrix<Scalar, 3, 1> Vec3;

 public:
  /**
   * @param origin: origin of ray
   * @param direction: *unit* direction of ray
   *
   */
  Ray3(const Vec3& origin,
       const Vec3& direction) :
         origin(origin),
         direction(direction.normalized()),
         tmin(Scalar(0)),
         tmax(std::numeric_limits<Scalar>::max()),
         invdir(Scalar(1.) / direction.array()),
         sign(invdir.x() < 0,
              invdir.y() < 0,
              invdir.z() < 0)
  {
  }

 public:
  Vec3 point_at(Scalar t) const {
    return origin + t*direction;
  }

 public:
  Vec3 origin, direction;
  Scalar tmin, tmax; /// ray min and max distances
  Vec3 invdir; // for convenience in AABB intersection
  std::tuple<int, int, int> sign;

};

typedef Ray2<float> Ray2f;
typedef Ray3<float> Ray3f;

} } /* ca */
#endif /* end of include guard: RAY_HPP_TVFS3AKQ */
