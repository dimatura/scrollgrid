/**
 * Copyright (c) 2015 Carnegie Mellon University, Daniel Maturana <dimatura@cmu.edu>
 *
 * For License information please see the LICENSE file in the root directory.
 *
 */

#ifndef BOX_HPP_4WFLKG9Q
#define BOX_HPP_4WFLKG9Q

#include <ros/ros.h>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace ca { namespace scrollgrid {

/**
 * Describes an axis-aligned volume in space.
 */
template<typename Scalar, int Dim>
class Box {
 public:
  typedef Eigen::Matrix<Scalar, Dim, 1> Vec;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  Box() :
      center_(Vec::Zero()),
      radius_(Vec::Zero()),
      bounds_(Vec::Zero(), Vec::Zero())
  {
  }

  Box(const Vec& min_pt,
      const Vec& max_pt) :
      center_( (min_pt+max_pt)/2 ),
      radius_( (max_pt-min_pt)/2 ),
      bounds_( min_pt, max_pt )
  {
    ROS_ASSERT(check_bounds());
  }

  Box(const Box& other) :
      center_(other.center_),
      radius_(other.radius_),
      bounds_(other.bounds_)
  {
    ROS_ASSERT(check_bounds());
  }

  Box& operator=(const Box& other) {
    if (this==&other) { return *this; }
    center_ = other.center_;
    radius_ = other.radius_;
    bounds_ = other.bounds_;
    return *this;
  }

  const Vec& min_pt() const { return bounds_.first; }

  const Vec& max_pt() const { return bounds_.second; }

  const Vec& bound(int ix) const {
    ROS_ASSERT (ix >= 0 && ix < 2);
    if (ix==0) { return bounds_.first; }
    return bounds_.second;
  }

  void set_bound(int ix, const Vec& v) {
    ROS_ASSERT (ix >= 0 && ix < 2);
    ROS_ASSERT( check_bounds());
    if (ix==0) {
      bounds_.first = v;
    } else {
      bounds_.second = v;
    }
  }

  void set_min_pt(const Vec& min_pt) {
    bounds_.first = min_pt;
    //ROS_ASSERT( check_bounds());
    this->reset_center_radius();
  }

  void set_max_pt(const Vec& max_pt) {
    bounds_.second = max_pt;
    //ROS_ASSERT( check_bounds());
    this->reset_center_radius();
  }

  const Vec& center() const {
    return center_;
  }

  const Vec& radius() const {
    return radius_;
  }

  void set_center(const Vec& v) {
    Vec delta = v - center_;
    this->translate(delta);
  }

  void set_radius(const Vec& v) {
    radius_ = v;
    bounds_.first = center_ - radius_;
    bounds_.second = center_ + radius_;
  }

  void translate(const Vec& v) {
    center_ += v;
    bounds_.first += v;
    bounds_.second += v;
  }

  bool contains(const Vec& v) const {
    return (bounds_.first.array() <= v.array()).all() &&
        (v.array() <= bounds_.second.array()).all();
  }

 private:
  // call after bounds change
  void reset_center_radius() {
    center_ = (bounds_.first+bounds_.second)/2;
    radius_ = (bounds_.second-bounds_.first)/2;
  }

  bool check_bounds() {
    return bounds_.first.x() <= bounds_.second.x() &&
           bounds_.first.y() <= bounds_.second.y() &&
           bounds_.first.z() <= bounds_.second.z();
  }

 private:
  Vec center_;
  Vec radius_;
  std::pair<Vec, Vec> bounds_;

};

}
} /* ca */

#endif /* end of include guard: BOX_HPP_4WFLKG9Q */
