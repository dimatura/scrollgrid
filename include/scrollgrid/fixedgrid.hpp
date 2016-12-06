/**
 * Copyright (c) 2015 Carnegie Mellon University, Daniel Maturana <dimatura@cmu.edu>
 *
 * For License information please see the LICENSE file in the root directory.
 *
 */

#ifndef FIXEDGRID_HPP_TL9WFUEA
#define FIXEDGRID_HPP_TL9WFUEA

#include <stdint.h>
#include <math.h>

#include <vector>
#include <algorithm>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ros/console.h>

#include <geom_cast/geom_cast.hpp>
#include <pcl_util/point_types.hpp>

#include "scrollgrid/grid_types.hpp"
#include "scrollgrid/box.hpp"

namespace ca { namespace scrollgrid {

template<class Scalar, int Dim>
class FixedGrid {
public:
  enum { dim = Dim };
  typedef Scalar ScalarType;
  typedef Eigen::Matrix<Scalar, Dim, 1> Vec;
  typedef Eigen::Matrix<Scalar, Dim, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<grid_ix_t, Dim, 1> VecIx;
  typedef Eigen::Matrix<grid_ix_t, Dim, 1> MatIx;

  typedef std::shared_ptr<FixedGrid> Ptr;
  typedef std::shared_ptr<const FixedGrid> ConstPtr;

public:
  FixedGrid() :
      box_(),
      origin_(Vec::Zero()),
      min_world_corner_ix_(Vec3Ix::Zero()),
      dimension_(Vec::Zero()),
      num_cells_(0),
      strides_(VecIx::Zero()),
      resolution_(0.)
  {
  }

  /**
   * @param center: center of the grid in global frame
   * @param dimension: number of grid cells along each coordinate
   * @param resolution: size of each grid cell side. they are cubic.
   * Default assumes so-called "ZYX" layout, i.e. z changes fastest.
   */
  FixedGrid(const Vec& center,
            const VecIx& dimension,
            Scalar resolution) :
      box_(center-(dimension.template cast<Scalar>()*resolution)/2,
           center+(dimension.template cast<Scalar>()*resolution)/2),
      origin_(center-box_.radius()),
      dimension_(dimension),
      num_cells_(dimension.prod()),
      resolution_(resolution) {

    this->calc_strides();
    this->calc_min_corner_ix();
  }

  FixedGrid(const FixedGrid& other) :
      box_(other.box_),
      origin_(other.origin_),
      min_world_corner_ix_(other.min_world_corner_ix_),
      dimension_(other.dimension_),
      num_cells_(other.num_cells_),
      strides_(other.strides_),
      resolution_(other.resolution_) {
  }

  FixedGrid& operator=(const FixedGrid& other) {
    if (*this==other) { return *this; }
    box_ = other.box_;
    origin_ = other.origin_;
    min_world_corner_ix_ = other.min_world_corner_ix_;
    dimension_ = other.dimension_;
    num_cells_ = other.num_cells_;
    strides_ = other.strides_;
    resolution_ = other.resolution_;
    return *this;
  }

  virtual ~FixedGrid() { }

public:
  /**
   * see ctor for params
   */
  void reset(const Vec& center,
             const VecIx& dimension,
             Scalar resolution) {
    box_.set_center(center);
    box_.set_radius((dimension.template cast<Scalar>()*resolution)/2);
    origin_ = center - box_.radius();
    dimension_ = dimension;
    num_cells_ = dimension.prod();
    resolution_ = resolution;
    this->calc_strides();
  }

  void copy_from(FixedGrid& other) {
    this->reset(other.center(),
                other.dimension(),
                other.resolution());
  }

  /**
   * Is pt inside 3D box containing grid?
   * @param pt point in same frame as center (probably world_view)
   */
  bool is_inside_box(const Vec& pt) const {
    return box_.contains(pt);
  }


  Eigen::VectorXi multiple_is_inside_box(const Mat& pts) {
    Eigen::VectorXi out(pts.cols());
    for (int i=0; i < pts.cols(); ++i) {
      const Eigen::Matrix<Scalar, Dim, 1>& pt(pts.col(i));
      int inside = box_.contains(pt);
      out[i] = inside;
    }
    return out;
  }


  template<class PointT>
  bool is_inside_box(const PointT& pt) const {
    return box_.contains(ca::point_cast<Vec>(pt));
  }

  /**
   * is i, j, (k) inside the grid limits?
   */
  bool is_inside_grid(const VecIx& grid_ix) const {
    return ((grid_ix.array() >= 0).all() &&
            (grid_ix.array() < dimension_.array()).all());
  }

  /**
   * Given position in world coordinates, return grid coordinates.
   * Note: does not check if point is inside grid.
   */
  VecIx world_to_grid(const Vec& x) const {
    Vec tmp = ((x - origin_).array() - 0.5*resolution_)/resolution_;
    //return tmp.cast<grid_ix_t>();
    return VecIx(tmp.round());
  }

  Vec grid_to_world(const VecIx& grid_ix) const {
    Vec w((grid_ix.template cast<Scalar>()*resolution_ + origin_).array() + 0.5*resolution_);
    return w;
  }

  Mat multiple_grid_to_world(const MatIx& grid_indices) {
    Mat out(3, grid_indices.cols());
    for (int ix=0; ix < grid_indices.cols(); ++ix) {
      Vec3Ix gix(grid_indices.col(ix));
      out.col(ix) = this->grid_to_world(gix);
    }
    return out;
  }

  /**
   * Note that this wraps the z dimension.
   * TODO c-order/f-order config
   */
  mem_ix_t grid_to_mem(const VecIx& grid_ix) const {
    return strides_.dot(grid_ix);
  }


  MemIxVector multiple_grid_to_mem(const MatIx& grid_indices) {
    MemIxVector out(grid_indices.cols());
    for (int ix=0; ix < grid_indices.cols(); ++ix) {
      Vec3Ix gix(grid_indices.col(ix));
      out[ix] = this->grid_to_mem(gix);
    }
    return out;
  }

  // TODO is it worth unrolling for 2/3 special case
  VecIx mem_to_grid(mem_ix_t mem_ix) const {
    VecIx out(VecIx::Zero());
    for (int i=0; i < Dim; ++i) {
      grid_ix_t ix = mem_ix/strides_[i];
      mem_ix -= ix*strides_[i];
      out(i) = ix;
    }
    return out;
  }

  /**
   * pack into grid_ix_t as [int16, int16, int16, 0]
   * note grid_ix_it is a signed 64-bit type
   * this gives range of [-32768, 32768] for each coordinate
   * so if voxel is 5cm, [-1638.4 m, 1638.4 m] relative to initial center.
   * this is useful as a unique hash.
   * this is better than mem_ix because mem_ix is ambiguous for absolute ijk.
   * if we use linear mem_ix as a hash,
   * then because of scrolling there may be collisions, and mem_ix become
   * invalidated once the corresponding voxel scrolls out.
   *
   * Note: we assume 3D here, even for 2D grids. Set grid_ix[2] to 0 for 2D.
   */
  hash_ix_t grid_to_hash(const Vec3Ix& grid_ix) const {
    VecIx grid_positive(grid_ix - min_world_corner_ix_);
    hash_ix_t hi = static_cast<hash_ix_t>(grid_positive[0]);
    hash_ix_t hj = static_cast<hash_ix_t>(grid_positive[1]);
    hash_ix_t hk = static_cast<hash_ix_t>(grid_positive[2]);
    hash_ix_t h = (hi << 48) | (hj << 32) | (hk << 16);
    return h;
  }

  Vec3Ix hash_to_grid(hash_ix_t hix) const {
    hash_ix_t hi = (hix & 0xffff000000000000) >> 48;
    hash_ix_t hj = (hix & 0x0000ffff00000000) >> 32;
    hash_ix_t hk = (hix & 0x00000000ffff0000) >> 16;
    Vec3Ix grid_ix(hi, hj, hk);
    grid_ix += min_world_corner_ix_;
    return grid_ix;
  }

  Mat3Ix multiple_hash_to_grid(const HashIxVector& hindices) const {
    Mat3Ix out(3, hindices.rows());
    for (size_t i=0; i < hindices.rows(); ++i) {
      Vec3Ix gix = this->hash_to_grid(hindices[i]);
      out.col(i) = gix;
    }
    return out;
  }

  /**
   * Like the above but does not offset by origin.
   * In the fixed case we assume the min ijk for origin is 0,0,0.
   */
  hash_ix_t local_grid_to_hash(const Vec3Ix& grid_ix) const {
    // TODO assumes grid_ix is inside box.
    hash_ix_t hi = static_cast<hash_ix_t>(grid_ix[0]);
    hash_ix_t hj = static_cast<hash_ix_t>(grid_ix[1]);
    hash_ix_t hk = static_cast<hash_ix_t>(grid_ix[2]);
    hash_ix_t h = (hi << 48) | (hj << 32) | (hk << 16);
    return h;
  }

  Vec3Ix hash_to_local_grid(hash_ix_t hix) const {
    hash_ix_t hi = (hix & 0xffff000000000000) >> 48;
    hash_ix_t hj = (hix & 0x0000ffff00000000) >> 32;
    hash_ix_t hk = (hix & 0x00000000ffff0000) >> 16;
    Vec3Ix grid_ix(hi, hj, hk);
    return grid_ix;
  }

  Mat3Ix multiple_local_hash_to_grid(const HashIxVector& hindices) const {
    Mat3Ix out(3, hindices.rows());
    for (size_t i=0; i < hindices.rows(); ++i) {
      Vec3Ix gix = this->local_hash_to_grid(hindices[i]);
      out.col(i) = gix;
    }
    return out;
  }

 public:
  grid_ix_t dim(int i) const { return dimension_[i]; }
  grid_ix_t first(int i) const { return 0; }
  grid_ix_t last(int i) const { return dimension_[i]; }

  const VecIx& dimension() const { return dimension_; }
  const Vec& radius() const { return box_.radius(); }
  const Vec& origin() const { return origin_; }
  Vec min_pt() const { return box_.min_pt(); }
  Vec max_pt() const { return box_.max_pt(); }
  const Vec& center() const { return box_.center(); }
  Scalar resolution() const { return resolution_; }
  ca::scrollgrid::Box<Scalar, Dim> box() const { return box_; }

  grid_ix_t num_cells() const { return num_cells_; }

 private:
  void calc_strides() {
    for (int i=0; i < Dim; ++i) {
      strides_(Dim-i-1) = dimension_.tail(i).prod();
    }
  }

  void calc_min_corner_ix() {
    // set the "minimum possible" ijk, assuming the grid
    // won't stray "too far" from the initial position.
    // too far == more than 2^15 grid cells.
    // so if you voxel resolution is 1 cm, 327.68 m.
    Scalar m = -static_cast<Scalar>(std::numeric_limits<uint16_t>::max()/2)*resolution_;
    Vec mn; mn.fill(m);
    mn += box_.center();
    min_world_corner_ix_ = this->world_to_grid(mn);
  }

 private:
  // 2d box enclosing grid. In whatever coordinates were given (probably
  // world_view)
  ca::scrollgrid::Box<Scalar, Dim> box_;

  // static origin of the grid coordinate system.
  Vec origin_;

  // minimum world corner in ijk. used for hash
  // Note that it is always 3D, even for 2D grids.
  // For 2D grids we ignore the third dimension.
  Vec3Ix min_world_corner_ix_;

  // number of grid cells along each axis
  VecIx dimension_;

  // number of cells
  grid_ix_t num_cells_;

  // grid strides to translate from linear to 3D layout.
  // C-ordering, ie x slowest, z fastest.
  VecIx strides_;

  // size of grid cells
  Scalar resolution_;
};

typedef FixedGrid<float, 2> FixedGrid2f;
typedef FixedGrid<float, 3> FixedGrid3f;

} }

#endif
