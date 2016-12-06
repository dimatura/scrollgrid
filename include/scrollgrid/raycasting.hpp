/**
 * Copyright (c) 2015 Carnegie Mellon University, Daniel Maturana <dimatura@cmu.edu>
 *
 * For License information please see the LICENSE file in the root directory.
 *
 */

#ifndef RAYCASTING_HPP_OLVFBMND
#define RAYCASTING_HPP_OLVFBMND

#include <Eigen/Core>

#include <pcl_util/point_types.hpp>

#include "scrollgrid/grid_types.hpp"
#include "scrollgrid/box.hpp"
#include "scrollgrid/ray.hpp"
#include "scrollgrid/scrollgrid3.hpp"
#include "scrollgrid/dense_array3.hpp"

namespace ca
{

/**
 * Axis-aligned bounding box intersection test.
 * Reference:
 * An Efficient and Robust Ray–Box Intersection Algorithm, Williams et al. 2004
 * tmin and tmax are updated in place
 */
template<typename Scalar>
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

class Bresenham2Iterator {
public:
  Bresenham2Iterator(const Vec2Ix& start_pos,
                     const Vec2Ix& end_pos) :
      _x0(start_pos.x()),
      _y0(start_pos.y()),
      _x1(end_pos.x()),
      _y1(end_pos.y()),
      _x(start_pos.x()),
      _y(start_pos.y()),
      _sx(0),
      _sy(0),
      _ax(0),
      _ay(0)
  {
    this->init();
  }

  void init() {
    _x = _x0;
    _y = _y0;

    _dx = _x1 - _x0;
    _dy = _y1 - _y0;

    //X
    if (_dx>0) {
      _sx = 1;
    } else if (_dx<0) {
      _sx = -1;
      _dx = -_dx;
    } else {
      _sx = 0;
    }

    //Y
    if (_dy>0) {
      _sy = 1;
    } else if (_dy<0) {
      _sy = -1;
      _dy = -_dy;
    } else {
      _sy = 0;
    }

    _ax = 2*_dx;
    _ay = 2*_dy;

    if (_dy <= _dx) {
      _decy = _ay-_dx;
    } else {
      _decx = _ax-_dy;
    }

    _done = false;
  }

  void step() {
    if (_dy <= _dx){
      if (_x == _x1) {
        // hit
        _done = true;
        return;
      } else {
        // pass
      }
      if (_decy >= 0) {
        _decy -= _ax;
        _y += _sy;
      }
      _x += _sx;
      _decy += _ay;
    } else {
      if (_y == _y1) {
        // hit
        _done = true;
        return;
      } else {
        // pass
      }

      if (_decx >= 0) {
        _decx -= _ay;
        _x += _sx;
      }
      _y += _sy;
      _decx += _ax;
    }
  }

  int x() const {
    return _x;
  }

  int y() const {
    return _y;
  }

  Vec2Ix xy() const {
    return Vec2Ix(_x, _y);
  }

  bool done() const {
    return _done;
  }

  virtual ~Bresenham2Iterator() { }
  Bresenham2Iterator(const Bresenham2Iterator& other) = delete;
  Bresenham2Iterator& operator=(const Bresenham2Iterator& other) = delete;

private:
  int _x0, _y0, _x1, _y1;
  int _x, _y;
  int _sx, _sy;
  int _ax, _ay;
  int _dx, _dy;
  int _decx, _decy;
  bool _done;
};

/**
 * Trace a straight line from start_pos to end_pos.
 * At each step fun(i, j, k) is called.
 *
 * NOTE start_pos and end_pos should be inside the grid
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


template<class GridT, class ArrayT>
void bresenham_trace3_hit_pass(const Vec3Ix& start_pos,
                           const Vec3Ix& end_pos,
                           const GridT& grid3,
                           ArrayT& array3) {
#define UPDATE_CELL_HIT
#define UPDATE_CELL_PASS \
  mem_ix_t mem_ix = grid3.grid_to_mem(x, y, z); \
  array3[mem_ix].pass += 1;
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
