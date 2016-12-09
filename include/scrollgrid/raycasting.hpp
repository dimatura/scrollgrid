/**
 * Copyright (c) 2015 Carnegie Mellon University, Daniel Maturana <dimatura@cmu.edu>
 *
 * For License information please see the LICENSE file in the root directory.
 *
 */

#ifndef RAYCASTING_HPP_OLVFBMND
#define RAYCASTING_HPP_OLVFBMND

#include <Eigen/Core>
#include <Eigen/Dense>

#include <pcl_util/point_types.hpp>

#include "scrollgrid/grid_types.hpp"
#include "scrollgrid/box.hpp"
#include "scrollgrid/ray.hpp"
#include "scrollgrid/scrollgrid3.hpp"

#include <ros/ros.h>

namespace ca
{

/**
 * Axis-aligned bounding box intersection test.
 * Reference:
 * An Efficient and Robust Ray–Box Intersection Algorithm, Williams et al. 2004
 * tmin and tmax are updated in place
 */
template <typename Scalar>
bool aabb_ray_intersect(const ca::scrollgrid::Box<Scalar, 3>& box,
                        ca::scrollgrid::Ray3<Scalar> &r) {
  Scalar tmin = (box.bound(   std::get<0>(r.sign) ).x() - r.origin.x()) * r.invdir.x();
  Scalar tmax = (box.bound( 1-std::get<0>(r.sign) ).x() - r.origin.x()) * r.invdir.x();

  Scalar tymin = (box.bound(  std::get<1>(r.sign) ).y() - r.origin.y()) * r.invdir.y();
  Scalar tymax = (box.bound(1-std::get<1>(r.sign) ).y() - r.origin.y()) * r.invdir.y();

  if ((tmin > tymax) || (tymin > tmax)) { return false; }
  if (tymin > tmin) { tmin = tymin; }
  if (tymax < tmax) { tmax = tymax; }

  Scalar tzmin = (box.bound(  std::get<2>(r.sign)).z() - r.origin.z()) * r.invdir.z();
  Scalar tzmax = (box.bound(1-std::get<2>(r.sign)).z() - r.origin.z()) * r.invdir.z();

  if ((tmin > tzmax) || (tzmin > tmax)) { return false; }
  if (tzmin > tmin) { tmin = tzmin; }
  if (tzmax < tmax) { tmax = tzmax; }
  if (tmin > r.tmin) { r.tmin = tmin; }
  if (tmax < r.tmax) { r.tmax = tmax; }
  return true;
}

/**
 * Axis-aligned bounding box intersection test.
 * Reference:
 * An Efficient and Robust Ray–Box Intersection Algorithm, Williams et al. 2004
 * tmin and tmax are updated in place
 */
template <typename Scalar>
bool aabb_ray_intersect(const ca::scrollgrid::Box<Scalar, 2>& box,
                        ca::scrollgrid::Ray2<Scalar> &r) {
  Scalar tmin = (box.bound(   std::get<0>(r.sign) ).x() - r.origin.x()) * r.invdir.x();
  Scalar tmax = (box.bound( 1-std::get<0>(r.sign) ).x() - r.origin.x()) * r.invdir.x();

  Scalar tymin = (box.bound(  std::get<1>(r.sign) ).y() - r.origin.y()) * r.invdir.y();
  Scalar tymax = (box.bound(1-std::get<1>(r.sign) ).y() - r.origin.y()) * r.invdir.y();

  if ((tmin > tymax) || (tymin > tmax)) { return false; }
  if (tymin > tmin) { tmin = tymin; }
  if (tymax < tmax) { tmax = tymax; }

  if (tmin > r.tmin) { r.tmin = tmin; }
  if (tmax < r.tmax) { r.tmax = tmax; }
  return true;
}

enum class OutcodeSide : int {
  left = (1 << 0),
  right = (1 << 1),
  bottom = (1 << 2),
  top = (1 << 3),
  near = (1 << 4),
  far = (1 << 5)
};


/**
 * outcode for cohen-sutherland
 */
template <typename Scalar>
int compute_outcode2(const Eigen::Matrix<Scalar, 2, 1>& pt,
                     const ca::scrollgrid::Box<Scalar, 2>& box) {

  Scalar min_x(box.min_pt().x());
  Scalar max_x(box.max_pt().x());
  Scalar min_y(box.min_pt().y());
  Scalar max_y(box.max_pt().y());

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
template <typename Scalar>
int compute_outcode3(const Eigen::Matrix<Scalar, 3, 1>& pt,
                     const ca::scrollgrid::Box<Scalar, 3>& box) {

  Scalar min_x(box.min_pt().x());
  Scalar max_x(box.max_pt().x());
  Scalar min_y(box.min_pt().y());
  Scalar max_y(box.max_pt().y());
  Scalar min_z(box.min_pt().z());
  Scalar max_z(box.max_pt().z());

  int code = 0;
  code |= (pt.x() > max_x)*static_cast<int>(OutcodeSide::right);
  code |= (pt.x() < min_x)*static_cast<int>(OutcodeSide::left);
  code |= (pt.y() > max_y)*static_cast<int>(OutcodeSide::top);
  code |= (pt.y() < min_y)*static_cast<int>(OutcodeSide::bottom);
  code |= (pt.z() > max_z)*static_cast<int>(OutcodeSide::far);
  code |= (pt.z() < min_z)*static_cast<int>(OutcodeSide::near);

  return code;
}


#if 0
template <typename Scalar>
struct ClipLineEps { };
template <>
struct ClipLineEps<float> { static constexpr float value = 1e-6; };
template <>
struct ClipLineEps<double> { static constexpr double value = 1e-6; };
template <>
struct ClipLineEps<grid_ix_t> { static constexpr grid_ix_t value = 1; };
#endif

template <typename Scalar>
bool clip_line2(const Eigen::Matrix<Scalar, 2, 1>& pt1,
                const Eigen::Matrix<Scalar, 2, 1>& pt2,
                const ca::scrollgrid::Box<Scalar, 2>& box,
                Eigen::Matrix<Scalar, 2, 1>& out_pt1,
                Eigen::Matrix<Scalar, 2, 1>& out_pt2) {

  bool accept = false;
  bool done = false;

  int codeout;

  Scalar x1(pt1.x());
  Scalar x2(pt2.x());
  Scalar y1(pt1.y());
  Scalar y2(pt2.y());
  Scalar min_x(box.min_pt().x());
  Scalar max_x(box.max_pt().x());
  Scalar min_y(box.min_pt().y());
  Scalar max_y(box.max_pt().y());

  int code1 = compute_outcode2({x1, y1}, box);
  int code2 = compute_outcode2({x2, y2}, box);

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
        code1 = compute_outcode2({x1, y1}, box);
      } else { //second endpoint was clipped
        x2 = x;
        y2 = y;
        code2 = compute_outcode2({x2, y2}, box);
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

template <typename Scalar>
bool clip_line3(const Eigen::Matrix<Scalar, 3, 1>& pt1,
                const Eigen::Matrix<Scalar, 3, 1>& pt2,
                const ca::scrollgrid::Box<Scalar, 3>& box,
                Eigen::Matrix<Scalar, 3, 1>& out_pt1,
                Eigen::Matrix<Scalar, 3, 1>& out_pt2) {

  bool accept = false;
  bool done = false;

  int codeout;

  Scalar x1(pt1.x());
  Scalar x2(pt2.x());
  Scalar y1(pt1.y());
  Scalar y2(pt2.y());
  Scalar z1(pt1.z());
  Scalar z2(pt2.z());
  Scalar min_x(box.min_pt().x());
  Scalar max_x(box.max_pt().x());
  Scalar min_y(box.min_pt().y());
  Scalar max_y(box.max_pt().y());
  Scalar min_z(box.min_pt().z());
  Scalar max_z(box.max_pt().z());

  int code1 = compute_outcode3({x1, y1, z1}, box);
  int code2 = compute_outcode3({x2, y2, z2}, box);

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
      } else if (codeout & static_cast<int>(OutcodeSide::near)) {
        Scalar t = (max_z - z1) / (z2 - z1);
        x = x1 + t*(x2 - x1);
        y = y1 + t*(y2 - y1);
        z = max_z;
      } else if (codeout & static_cast<int>(OutcodeSide::far)) {
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
        code1 = compute_outcode3({x1, y1, z1}, box);
      } else {
        // second endpoint was clipped
        x2 = x;
        y2 = y;
        z2 = z;
        code2 = compute_outcode3({x2, y2, z2}, box);
      }
    }

    ++ctr;
    if (ctr > 6) {
      ROS_FATAL("clipping failed");
      ROS_ASSERT(false);
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

class Bresenham2Iterator {
public:
  Bresenham2Iterator(const Vec2Ix& start_pos,
                     const Vec2Ix& end_pos) :
      x0_(start_pos.x()),
      y0_(start_pos.y()),
      x1_(end_pos.x()),
      y1_(end_pos.y()),
      x_(start_pos.x()),
      y_(start_pos.y()),
      sx_(0),
      sy_(0),
      ax_(0),
      ay_(0)
  {
    this->init();
  }

  void init() {
    x_ = x0_;
    y_ = y0_;

    dx_ = x1_ - x0_;
    dy_ = y1_ - y0_;

    //X
    if (dx_>0) {
      sx_ = 1;
    } else if (dx_<0) {
      sx_ = -1;
      dx_ = -dx_;
    } else {
      sx_ = 0;
    }

    //Y
    if (dy_>0) {
      sy_ = 1;
    } else if (dy_<0) {
      sy_ = -1;
      dy_ = -dy_;
    } else {
      sy_ = 0;
    }

    ax_ = 2*dx_;
    ay_ = 2*dy_;

    if (dy_ <= dx_) {
      decy_ = ay_-dx_;
    } else {
      decx_ = ax_-dy_;
    }

    done_ = false;
  }

  void step() {
    if (dy_ <= dx_){
      if (done_) {
        return;
      }
      if (decy_ >= 0) {
        decy_ -= ax_;
        y_ += sy_;
      }
      x_ += sx_;
      decy_ += ay_;
    } else {
      if (done_) {
        return;
      }
      if (decx_ >= 0) {
        decx_ -= ay_;
        x_ += sx_;
      }
      y_ += sy_;
      decx_ += ax_;
    }

    if (x_ == x1_ || y_ == y1_) {
      done_ = true;
    }

  }

  int i() const {
    return x_;
  }

  int j() const {
    return y_;
  }

  Vec2Ix pos() const {
    return Vec2Ix(x_, y_);
  }

  bool done() const {
    return done_;
  }

  virtual ~Bresenham2Iterator() { }
  Bresenham2Iterator(const Bresenham2Iterator& other) = delete;
  Bresenham2Iterator& operator=(const Bresenham2Iterator& other) = delete;

private:
  int x0_, y0_;
  int x1_, y1_;
  int x_, y_;
  int sx_, sy_;
  int ax_, ay_;
  int dx_, dy_;
  int decx_, decy_;
  bool done_;
};

class Bresenham3Iterator {
public:
  Bresenham3Iterator(const Vec3Ix& start_pos,
                     const Vec3Ix& end_pos) :
      x0_(start_pos.x()),
      y0_(start_pos.y()),
      z0_(start_pos.z()),
      x1_(end_pos.x()),
      y1_(end_pos.y()),
      z1_(end_pos.z())
  {
    this->init();
  }

  void init() {
    x_ = x0_;
    y_ = y0_;
    z_ = z0_;
    dx_ = x1_ - x0_,
    dy_ = y1_ - y0_,
    dz_ = z1_ - z0_;

    //x_
    if (dx_ > 0) {
      sx_ = 1;
    } else if (dx_ < 0) {
      sx_ = -1;
      dx_ = -dx_;
    } else {
      sx_ = 0;
    }

    //y_
    if (dy_ > 0) {
      sy_ = 1;
    } else if (dy_ < 0) {
      sy_ = -1;
      dy_ = -dy_;
    } else {
      sy_ = 0;
    }

    //z_
    if (dz_ > 0) {
      sz_ = 1;
    } else if (dz_ < 0) {
      sz_ = -1;
      dz_ = -dz_;
    } else {
      sz_ = 0;
    }

    ax_ = 2*dx_;
    ay_ = 2*dy_;
    az_ = 2*dz_;

    if ((dy_ <= dx_) && (dz_ <= dx_)) {
      decy_ = ay_-dx_;
      decz_ = az_-dx_;
    } else if ((dx_ <= dy_) && (dz_ <= dy_)) {
      decx_ = ax_-dy_;
      decz_ = az_-dy_;
    } else if ((dx_ <= dz_) && (dy_ <= dz_)) {
      decx_ = ax_-dz_;
      decy_ = ay_-dz_;
    }
  }

  void step() {
    if ((dy_ <= dx_) && (dz_ <= dx_)) {
      if (done_) {
        return;
      }
      if (decy_ >= 0) {
        decy_ -= ax_;
        y_ += sy_;
      }
      if (decz_ >= 0) {
        decz_ -= ax_;
        z_ += sz_;
      }

      x_ += sx_; decy_ += ay_; decz_ += az_;
    } else if ((dx_ <= dy_) && (dz_ <= dy_)) {
      if (done_) {
        return;
      }
      if (decx_ >= 0) {
        decx_ -= ay_;
        x_ += sx_;
      }
      if (decz_ >= 0) {
        decz_ -= ay_;
        z_ += sz_;
      }
      y_ += sy_; decx_ += ax_; decz_ += az_;

    } else if ((dx_ <= dz_) && (dy_ <= dz_)) {
      if (done_) {
        return;
      }
      if (decx_ >= 0) {
        decx_ -= az_;
        x_ += sx_;
      }
      if (decy_ >= 0) {
        decy_ -= az_;
        y_ += sy_;
      }
      z_ += sz_; decx_ += ax_; decy_ += ay_;
    }

    if (x_ == x1_ || y_ == y1_ || z_ == z1_) {
      done_ = true;
    }

  }

  bool done() {
    return done_;
  }

  Vec3Ix pos() { return Vec3Ix(x_, y_, z_); }

  virtual ~Bresenham3Iterator() { }
  Bresenham3Iterator(const Bresenham3Iterator& other) = delete;
  Bresenham3Iterator& operator=(const Bresenham3Iterator& other) = delete;

private:
  int x0_, y0_, z0_;
  int x1_, y1_, z1_;
  int x_, y_, z_;
  int sx_, sy_, sz_;
  int ax_, ay_, az_;
  int dx_, dy_, dz_;
  int decx_, decy_, decz_;
  bool done_;
};

/**
 * Trace a straight line from start_pos to end_pos.
 * At each step fun(i, j, k) is called.
 *
 * NOTE start_pos and end_pos should be inside the grid.
 *
 * Reference: graphics gems article
 * TODO consider DDA-type raytracing.
 */
template<class TraceFunctor>
void bresenham_trace3(const Vec3Ix& start_pos,
                      const Vec3Ix& end_pos,
                      const TraceFunctor& fun) {
#define UPDATE_CELL_HIT
#define UPDATE_CELL_PASS \
  if (!fun(x, y, z)) { break; }
#include "bresenham3_macro.def"
#undef UPDATE_CELL_HIT
#undef UPDATE_CELL_PASS
}


template<class GridT, class ArrayT>
void bresenham_trace3_increment(const Vec3Ix& start_pos,
                                const Vec3Ix& end_pos,
                                const GridT& grid3,
                                ArrayT& array3) {
#define UPDATE_CELL_HIT
#define UPDATE_CELL_PASS \
  mem_ix_t mem_ix = grid3.grid_to_mem(x, y, z); \
  array3[mem_ix] += 1;
#include "bresenham3_macro.def"
#undef UPDATE_CELL_HIT
#undef UPDATE_CELL_PASS
}


template<class TraceFunctor>
void bresenham_trace2(const Vec2Ix& start_pos,
                      const Vec2Ix& end_pos,
                      const TraceFunctor& fun) {

  int x = start_pos[0], y = start_pos[1];

  int dx = end_pos[0] - start_pos[0],
      dy = end_pos[1] - start_pos[1];

  int sx, sy;

  //X
  if (dx>0) {
    sx = 1;
  } else if (dx<0) {
    sx = -1;
    dx = -dx;
  } else {
    sx = 0;
  }

  //Y
  if (dy>0) {
    sy = 1;
  } else if (dy<0) {
    sy = -1;
    dy = -dy;
  } else {
    sy = 0;
  }

  int ax = 2*dx,
      ay = 2*dy;

  if (dy <= dx){
    for (int decy=ay-dx;;
         x+=sx, decy+=ay) {
      bool end_cell = false;
      //Bresenham step
      if(!fun(x, y, end_cell)) {
        break;
      }
      if (x==end_pos[0]) {
        break;
      }
      if (decy>=0) {
        decy-=ax;
        y+=sy;
      }
    }
  } else if (dx <= dy){
    for (int decx=ax-dy;;
         y+=sy,decx+=ax) {
      bool end_cell = false;
      //Bresenham step
      if (!fun(x,y,end_cell)){
        break;
      }

      if (y==end_pos[1]) break;
      if (decx>=0) {
        decx-=ay;
        x+=sx;
      }
    }
  }
}

}

#endif /* end of include guard: RAYCASTING_HPP_OLVFBMND */
