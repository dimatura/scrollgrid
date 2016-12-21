/**
 * Copyright (c) 2015 Carnegie Mellon University, Daniel Maturana <dimatura@cmu.edu>
 *
 * For License information please see the LICENSE file in the root directory.
 *
 */

#ifndef BOX_HPP_4WFLKG9Q
#define BOX_HPP_4WFLKG9Q

#include <type_traits>

#include <ros/ros.h>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace ca { namespace scrollgrid {

template<class Scalar> class Ray3;

/**
 * outcode for cohen-sutherland clipping.
 */
enum class OutcodeSide : int {
  left = (1 << 0),
  right = (1 << 1),
  bottom = (1 << 2),
  top = (1 << 3),
  near = (1 << 4),
  far = (1 << 5)
};


/**
 * Describes an axis-aligned volume in space.
 */
template<typename ScalarT, int dim>
class Box {
 public:
  static constexpr int Dim = dim;
  typedef std::integral_constant<int, Dim> dim_t;
  typedef std::integral_constant<int, 2> two_t;
  typedef std::integral_constant<int, 3> three_t;
  typedef ScalarT Scalar;
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

  bool clip_line(const Eigen::Matrix<Scalar, 2, 1>& pt1,
                 const Eigen::Matrix<Scalar, 2, 1>& pt2,
                 Eigen::Matrix<Scalar, 2, 1>& out_pt1,
                 Eigen::Matrix<Scalar, 2, 1>& out_pt2) {

    static_assert(std::is_same<dim_t, two_t>::value, "wrong dim on clip_line");

    bool accept = false;
    bool done = false;

    int codeout(0);

    Scalar x1(pt1.x());
    Scalar x2(pt2.x());
    Scalar y1(pt1.y());
    Scalar y2(pt2.y());
    Scalar min_x(this->min_pt().x());
    Scalar max_x(this->max_pt().x());
    Scalar min_y(this->min_pt().y());
    Scalar max_y(this->max_pt().y());

    int code1 = this->compute_outcode({x1, y1}, dim_t());
    int code2 = this->compute_outcode({x2, y2}, dim_t());

    //constexpr Scalar eps = ClipLineEps<Scalar>::value;

    int ctr = 0;

    do {
      if (!(code1 | code2)) {
        // accept because both endpoints are in screen or on the border, trivial accept
        accept = done = true;
      } else if (code1 & code2) {
        //the line isn't visible on screen, trivial reject
        done = true;
      } else {
        //if no trivial reject or accept, continue the loop
        Scalar x, y;
        codeout = code1 ? code1 : code2;
        if (codeout & static_cast<int>(OutcodeSide::top)) {
          x = x1 + (x2 - x1) * (max_y - y1) / (y2 - y1);
          y = max_y;
        } else if (codeout & static_cast<int>(OutcodeSide::bottom)) {
          x = x1 + (x2 - x1) * (min_y -y1) / (y2 - y1);
          y = min_y;
        } else if (codeout & static_cast<int>(OutcodeSide::right)) {
          y = y1 + (y2 - y1) * (max_x - x1) / (x2 - x1);
          x = max_x;
        } else { // intersect left
          y = y1 + (y2 - y1) * (min_x - x1) / (x2 - x1);
          x = min_x;
        }

        if (codeout == code1) { //first endpoint was clipped
          x1 = x;
          y1 = y;
          code1 = this->compute_outcode({x1, y1}, dim_t());
        } else { //second endpoint was clipped
          x2 = x;
          y2 = y;
          code2 = this->compute_outcode({x2, y2}, dim_t());
        }
      }

      ++ctr;
      if (ctr > 4) {
        ROS_FATAL("clipping failed");
        break;
      }
    } while (!done);

    if (accept) {
      out_pt1.x() = x1;
      out_pt2.x() = x2;
      out_pt1.y() = y1;
      out_pt2.y() = y2;
      return true;
    }

    out_pt1.setZero();
    out_pt2.setZero();

    return false;
  }

  bool clip_line(const Eigen::Matrix<Scalar, 3, 1>& pt1,
                 const Eigen::Matrix<Scalar, 3, 1>& pt2,
                 Eigen::Matrix<Scalar, 3, 1>& out_pt1,
                 Eigen::Matrix<Scalar, 3, 1>& out_pt2) {

    static_assert(std::is_same<dim_t, three_t>::value, "wrong dim on clip_line");

    bool accept = false;
    bool done = false;

    int codeout;

    Scalar x1(pt1.x());
    Scalar x2(pt2.x());
    Scalar y1(pt1.y());
    Scalar y2(pt2.y());
    Scalar z1(pt1.z());
    Scalar z2(pt2.z());
    Scalar min_x(this->min_pt().x());
    Scalar max_x(this->max_pt().x());
    Scalar min_y(this->min_pt().y());
    Scalar max_y(this->max_pt().y());
    Scalar min_z(this->min_pt().z());
    Scalar max_z(this->max_pt().z());

    int code1 = this->compute_outcode({x1, y1, z1}, dim_t());
    int code2 = this->compute_outcode({x2, y2, z2}, dim_t());

    //constexpr Scalar eps = ClipLineEps<Scalar>::value;

    int ctr = 0;

    do {
      if (!(code1 | code2)) {
        // accept because both endpoints are in screen or on the border, trivial accept
        accept = done = true;
      } else if (code1 & code2) {
        //the line isn't visible on screen, trivial reject
        done = true;
      } else {
        //if no trivial reject or accept, continue the loop
        Scalar x, y, z;
        codeout = code1 ? code1 : code2;
        if (codeout & static_cast<int>(OutcodeSide::top)) {
          Scalar t = (max_y - y1) / (y2 - y1);
          x = x1 + t*(x2 - x1);
          y = max_y;
          z = z1 + t*(z2 - z1);
        } else if (codeout & static_cast<int>(OutcodeSide::bottom)) {
          Scalar t = (min_y - y1) / (y2 - y1);
          x = x1 + t*(x2 - x1);
          y = min_y;
          z = z1 + t*(z2 - z1);
        } else if (codeout & static_cast<int>(OutcodeSide::right)) {
          Scalar t = (max_x - x1) / (x2 - x1);
          x = max_x;
          y = y1 + t*(y2 - y1);
          z = z1 + t*(z2 - z1);
        } else if (codeout & static_cast<int>(OutcodeSide::left)) {
          // TODO is this better than just doing it vectorially
          Scalar t = (min_x - x1) / (x2 - x1);
          x = min_x;
          y = y1 + t*(y2 - y1);
          z = z1 + t*(z2 - z1);
        } else if (codeout & static_cast<int>(OutcodeSide::far)) {
          Scalar t = (max_z - z1) / (z2 - z1);
          x = x1 + t*(x2 - x1);
          y = y1 + t*(y2 - y1);
          z = max_z;
        } else if (codeout & static_cast<int>(OutcodeSide::near)) {
          Scalar t = (min_z - z1) / (z2 - z1);
          x = x1 + t*(x2 - x1);
          y = y1 + t*(y2 - y1);
          z = min_z;
        } else {
          ROS_FATAL("bad outcode");
          ROS_ASSERT(false);
        }

        if (codeout == code1) {
          // first endpoint was clipped
          x1 = x;
          y1 = y;
          z1 = z;
          code1 = this->compute_outcode({x1, y1, z1}, dim_t());
        } else {
          // second endpoint was clipped
          x2 = x;
          y2 = y;
          z2 = z;
          code2 = this->compute_outcode({x2, y2, z2}, dim_t());
        }
      }

      ++ctr;
      if (ctr > 6) {
        ROS_FATAL("clipping failed");
        break;
      }
    } while (!done);

    if (accept) {
      out_pt1.x() = x1;
      out_pt2.x() = x2;
      out_pt1.y() = y1;
      out_pt2.y() = y2;
      out_pt1.z() = z1;
      out_pt2.z() = z2;
      return true;
    }

    out_pt1.setZero();
    out_pt2.setZero();

    return false;
  }



 /**
   * Axis-aligned bounding box intersection test.
   * Reference:
   * An Efficient and Robust Ray–Box Intersection Algorithm, Williams et al. 2004
   * tmin and tmax are updated in place
   */
  bool aabb_ray_intersect(ca::scrollgrid::Ray3<Scalar> &r) {
    Scalar tmin = (this->bound(   std::get<0>(r.sign) ).x() - r.origin.x()) * r.invdir.x();
    Scalar tmax = (this->bound( 1-std::get<0>(r.sign) ).x() - r.origin.x()) * r.invdir.x();

    Scalar tymin = (this->bound(  std::get<1>(r.sign) ).y() - r.origin.y()) * r.invdir.y();
    Scalar tymax = (this->bound(1-std::get<1>(r.sign) ).y() - r.origin.y()) * r.invdir.y();

    if ((tmin > tymax) || (tymin > tmax)) { return false; }
    if (tymin > tmin) { tmin = tymin; }
    if (tymax < tmax) { tmax = tymax; }

    Scalar tzmin = (this->bound(  std::get<2>(r.sign)).z() - r.origin.z()) * r.invdir.z();
    Scalar tzmax = (this->bound(1-std::get<2>(r.sign)).z() - r.origin.z()) * r.invdir.z();

    if ((tmin > tzmax) || (tzmin > tmax)) { return false; }
    if (tzmin > tmin) { tmin = tzmin; }
    if (tzmax < tmax) { tmax = tzmax; }
    if (tmin > r.tmin) { r.tmin = tmin; }
    if (tmax < r.tmax) { r.tmax = tmax; }
    return true;
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

  int compute_outcode(const Eigen::Matrix<Scalar, 2, 1>& pt, two_t tag) {

    Scalar min_x(this->min_pt().x());
    Scalar max_x(this->max_pt().x());
    Scalar min_y(this->min_pt().y());
    Scalar max_y(this->max_pt().y());

    int code = 0;
    code |= (pt.x() > max_x)*static_cast<int>(OutcodeSide::right);
    code |= (pt.x() < min_x)*static_cast<int>(OutcodeSide::left);
    code |= (pt.y() > max_y)*static_cast<int>(OutcodeSide::top);
    code |= (pt.y() < min_y)*static_cast<int>(OutcodeSide::bottom);

    return code;
  }

  /**
   * outcode for cohen-sutherland 3d
   */
  int compute_outcode(const Eigen::Matrix<Scalar, 3, 1>& pt, three_t tag) {
    Scalar min_x(this->min_pt().x());
    Scalar max_x(this->max_pt().x());
    Scalar min_y(this->min_pt().y());
    Scalar max_y(this->max_pt().y());
    Scalar min_z(this->min_pt().z());
    Scalar max_z(this->max_pt().z());

    int code = 0;
    code |= (pt.x() > max_x)*static_cast<int>(OutcodeSide::right);
    code |= (pt.x() < min_x)*static_cast<int>(OutcodeSide::left);
    code |= (pt.y() > max_y)*static_cast<int>(OutcodeSide::top);
    code |= (pt.y() < min_y)*static_cast<int>(OutcodeSide::bottom);
    code |= (pt.z() > max_z)*static_cast<int>(OutcodeSide::far);
    code |= (pt.z() < min_z)*static_cast<int>(OutcodeSide::near);

    return code;
  }


 private:
  Vec center_;
  Vec radius_;
  std::pair<Vec, Vec> bounds_;

};

typedef Box<float, 2> Box2f;
typedef Box<float, 3> Box3f;

} } /* ca */

#endif /* end of include guard: BOX_HPP_4WFLKG9Q */
