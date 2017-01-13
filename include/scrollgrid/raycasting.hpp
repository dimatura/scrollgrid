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

#include "scrollgrid/grid_types.hpp"

#include <ros/ros.h>

namespace ca { namespace scrollgrid {


template <int Dim>
class BresenhamIterator {
};

template <>
class BresenhamIterator<2> {
public:
  BresenhamIterator(const Vec2Ix& start_pos,
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
    if (done_) {
      return;
    }
    if (dy_ <= dx_){
      if (decy_ >= 0) {
        decy_ -= ax_;
        y_ += sy_;
      }
      x_ += sx_;
      decy_ += ay_;
    } else {
      if (decx_ >= 0) {
        decx_ -= ay_;
        x_ += sx_;
      }
      y_ += sy_;
      decx_ += ax_;
    }

    if (x_ == x1_ && y_ == y1_) {
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

  virtual ~BresenhamIterator() { }
  BresenhamIterator(const BresenhamIterator& other) = delete;
  BresenhamIterator& operator=(const BresenhamIterator& other) = delete;

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

template <>
class BresenhamIterator<3> {
public:
  BresenhamIterator(const Vec3Ix& start_pos,
                    const Vec3Ix& end_pos) :
      x0_(start_pos.x()),
      y0_(start_pos.y()),
      z0_(start_pos.z()),
      x_(start_pos.x()),
      y_(start_pos.y()),
      z_(start_pos.z()),
      x1_(end_pos.x()),
      y1_(end_pos.y()),
      z1_(end_pos.z()),
      sx_(0),
      sy_(0),
      sz_(0),
      ax_(0),
      ay_(0),
      az_(0)
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

    done_ = false;
  }

  void step() {
    if (done_) {
      return;
    }
    if ((dy_ <= dx_) && (dz_ <= dx_)) {
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

    if (x_ == x1_ && y_ == y1_ && z_ == z1_) {
      done_ = true;
    }

  }

  bool done() {
    return done_;
  }

  Vec3Ix pos() { return Vec3Ix(x_, y_, z_); }

  virtual ~BresenhamIterator() { }
  BresenhamIterator(const BresenhamIterator& other) = delete;
  BresenhamIterator& operator=(const BresenhamIterator& other) = delete;

private:
  int x0_ = 0, y0_ = 0, z0_ = 0;
  int x1_ = 0, y1_ = 0, z1_ = 0;
  int x_ = 0, y_ = 0, z_ = 0;
  int sx_ = 0, sy_ = 0, sz_ = 0;
  int ax_ = 0, ay_ = 0, az_ = 0;
  int dx_ = 0, dy_ = 0, dz_ = 0;
  int decx_ = 0, decy_ = 0, decz_ = 0;
  bool done_ = false;
};


} }

#endif /* end of include guard: RAYCASTING_HPP_OLVFBMND */
