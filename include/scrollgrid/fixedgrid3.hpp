/**
 * @author  Daniel Maturana
 * @year    2015
 *
 * @attention Copyright (c) 2015
 * @attention Carnegie Mellon University
 * @attention All rights reserved.
 *
 **@=*/


#ifndef FIXEDGRID3_HPP_6KPVVRZF
#define FIXEDGRID3_HPP_6KPVVRZF

#include <math.h>
#include <stdint.h>

#include <vector>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ros/console.h>

#include <pcl_util/point_types.hpp>
#include <geom_cast/geom_cast.hpp>

#include "scrollgrid/grid_types.hpp"
#include "scrollgrid/box.hpp"

namespace ca
{

template<class Scalar>
class FixedGrid3 {
public:
  typedef Scalar ScalarType; // TODO what is my convention for this?
  typedef Eigen::Matrix<Scalar, 3, 1> Vec3;
  typedef Eigen::Matrix<Scalar, 3, Eigen::Dynamic> Mat3;

  typedef boost::shared_ptr<FixedGrid3> Ptr;
  typedef boost::shared_ptr<const FixedGrid3> ConstPtr;

public:
  FixedGrid3() :
      box_(),
      origin_(0, 0, 0),
      min_world_corner_ijk_(0, 0, 0),
      dimension_(0, 0, 0),
      num_cells_(0),
      strides_(0, 0, 0),
      resolution_(0)
  { }

  /**
   * @param center: center of the grid in global frame
   * @param dimension: number of grid cells along each coordinate
   * @param resolution: size of each grid cell side. they are cubic.
   * Default assumes "ZYX" layout, i.e. z changes fastest.
   */
  FixedGrid3(const Vec3& center,
             const Vec3Ix& dimension,
             Scalar resolution,
             bool x_fastest=false) :
      box_(center-(dimension.cast<Scalar>()*resolution)/2,
           center+(dimension.cast<Scalar>()*resolution)/2),
      origin_(center-box_.radius()),
      dimension_(dimension),
      num_cells_(dimension.prod()),
      resolution_(resolution)
  {
    if (x_fastest) {
      strides_ = Vec3Ix(1, dimension[0], dimension.head<2>().prod());
    } else {
      strides_ = Vec3Ix(dimension.tail<2>().prod(), dimension[2], 1);
    }
    this->calc_min_corner_ijk();

  }

  virtual ~FixedGrid3() { }

  FixedGrid3(const FixedGrid3& other) :
      box_(other.box_),
      origin_(other.origin_),
      min_world_corner_ijk_(other.min_world_corner_ijk_),
      dimension_(other.dimension_),
      num_cells_(other.num_cells_),
      strides_(other.strides_),
      resolution_(other.resolution_)

  {
  }

  FixedGrid3& operator=(const FixedGrid3& other) {
    if (this==&other) { return *this; }
    box_ = other.box_;
    origin_ = other.origin_;
    min_world_corner_ijk_ = other.min_world_corner_ijk_;
    dimension_ = other.dimension_;
    num_cells_ = other.num_cells_;
    strides_ = other.strides_;
    resolution_ = other.resolution_;
    return *this;
  }

public:

  /**
   * see ctor for params
   */
  void reset(const Vec3& center,
             const Vec3Ix& dimension,
             Scalar resolution,
             bool x_fastest=false) {
    box_.set_center(center);
    box_.set_radius((dimension.template cast<Scalar>()*resolution)/2);
    origin_ = center - box_.radius();
    dimension_ = dimension;
    num_cells_ = dimension.prod();
    if (x_fastest) {
      strides_ = Vec3Ix(1, dimension[0], dimension.head<2>().prod());
    } else {
      strides_ = Vec3Ix(dimension.tail<2>().prod(), dimension[2], 1);
    }
    resolution_ = resolution;
  }


  /**
   * Is pt inside 3D box containing grid?
   * @param pt point in same frame as center (probably world_view)
   */
  bool is_inside_box(const Vec3& pt) const {
    return box_.contains(pt);
  }

  /**
   */
  Eigen::VectorXi multiple_is_inside_box(const Eigen::Matrix<Scalar, 3, Eigen::Dynamic>& pts) {
    Eigen::VectorXi out(pts.cols());
    for (int i=0; i < pts.cols(); ++i) {
      const Eigen::Matrix<Scalar, 3, 1>& pt(pts.col(i));
      int inside = box_.contains(pt);
      out[i] = inside;
    }
    return out;
  }

  bool is_inside_box(Scalar x, Scalar y, Scalar z) const {
    return box_.contains(Vec3(x, y, z));
  }

  template<class PointT>
  bool is_inside_box(const PointT& pt) const {
    return box_.contains(ca::point_cast<Eigen::Vector3d>(pt));
  }

  /**
   * is i, j, k inside the grid limits?
   */
  bool is_inside_grid(const Vec3Ix& grid_ix) const {
    return ((grid_ix.array() >= 0).all() &&
            (grid_ix.array() < dimension_.array()).all());
  }

  bool is_inside_grid(grid_ix_t i, grid_ix_t j, grid_ix_t k) const {
    return this->is_inside_grid(Vec3Ix(i, j, k));
  }

  /**
   * Given position in world coordinates, return ijk grid coordinates.
   * Note: does not check if point is inside grid.
   */
  Vec3Ix world_to_grid(const Vec3& xyz) const {
    Vec3 tmp = ((xyz - origin_).array() - 0.5*resolution_)/resolution_;
    //return tmp.cast<grid_ix_t>();
    return Vec3Ix(round(tmp.x()), round(tmp.y()), round(tmp.z()));
  }

  Vec3Ix world_to_grid(Scalar x, Scalar y, Scalar z) const {
    return this->world_to_grid(Vec3(x, y, z));
  }

  /**
   * Given ijk grid coordinates resurn xyz world coordinates.
   * xyz is center of corresponding voxel.
   */
  Vec3 grid_to_world(const Vec3Ix& grid_ix) const {
    Vec3 w((grid_ix.cast<Scalar>()*resolution_ + origin_).array() + 0.5*resolution_);
    return w;
  }

  Vec3 grid_to_world(grid_ix_t i, grid_ix_t j, grid_ix_t k) const {
    return this->grid_to_world(Vec3Ix(i, j, k));
  }

  Mat3 multiple_grid_to_world(const Mat3Ix& grid_indices) {

    Mat3 out( 3, grid_indices.cols() );
    for (int ix=0; ix < grid_indices.cols(); ++ix) {
      Vec3Ix gix(grid_indices.col(ix));
      out.col(ix) = this->grid_to_world( gix );
    }
    return out;

  }

  /**
   * given grid ijk coordinate, return linear memory index.
   */
  mem_ix_t grid_to_mem(grid_ix_t i, grid_ix_t j, grid_ix_t k) const {
    return this->grid_to_mem(Vec3Ix(i, j, k));
  }

  mem_ix_t grid_to_mem(const Vec3Ix& grid_ix) const {
    return strides_.dot(grid_ix);
  }

  MemIxVector multiple_grid_to_mem(const Mat3Ix& grid_indices) {
    MemIxVector out( grid_indices.cols() );
    for (int ix=0; ix < grid_indices.cols(); ++ix) {
      Vec3Ix gix(grid_indices.col(ix));
      out[ix] = this->grid_to_mem( gix );
    }
    return out;
  }

  /**
   * given linear memory index, return ijk coordinate.
   */
  Vec3Ix mem_to_grid(mem_ix_t mem_ix) const {
    grid_ix_t i = mem_ix/strides_[0];
    mem_ix -= i*strides_[0];
    grid_ix_t j = mem_ix/strides_[1];
    mem_ix -= j*strides_[1];
    grid_ix_t k = mem_ix;

    return Vec3Ix(i, j, k);
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
   */
  hash_ix_t grid_to_hash(const Vec3Ix& grid_ix) const {
    // grid2 should be all positive
    Vec3Ix grid2(grid_ix - min_world_corner_ijk_);

    hash_ix_t hi = static_cast<hash_ix_t>(grid2[0]);
    hash_ix_t hj = static_cast<hash_ix_t>(grid2[1]);
    hash_ix_t hk = static_cast<hash_ix_t>(grid2[2]);
    hash_ix_t h = (hi << 48) | (hj << 32) | (hk << 16);
    return h;
  }

  Vec3Ix hash_to_grid(hash_ix_t hix) const {
    hash_ix_t hi = (hix & 0xffff000000000000) >> 48;
    hash_ix_t hj = (hix & 0x0000ffff00000000) >> 32;
    hash_ix_t hk = (hix & 0x00000000ffff0000) >> 16;
    Vec3Ix grid_ix(hi, hj, hk);
    grid_ix += min_world_corner_ijk_;
    return grid_ix;
  }

  Mat3Ix multiple_hash_to_grid(const HashIxVector& hindices) const {
    Mat3Ix out( 3, hindices.rows() );
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
    Mat3Ix out( 3, hindices.rows() );
    for (size_t i=0; i < hindices.rows(); ++i) {
      Vec3Ix gix = this->local_hash_to_grid(hindices[i]);
      out.col(i) = gix;
    }
    return out;
  }

 public:
  const ca::scrollgrid::Box<Scalar, 3>& box() const { return box_; }
  grid_ix_t dim_i() const { return dimension_[0]; }
  grid_ix_t dim_j() const { return dimension_[1]; }
  grid_ix_t dim_k() const { return dimension_[2]; }
  // this is the case because of the way origin is set at construction.
  // TODO should we allow arbitrary origins?
  grid_ix_t first_i() const { return 0; }
  grid_ix_t first_j() const { return 0; }
  grid_ix_t first_k() const { return 0; }
  grid_ix_t last_i() const { return dimension_[0]; }
  grid_ix_t last_j() const { return dimension_[1]; }
  grid_ix_t last_k() const { return dimension_[2]; }
  const Vec3Ix& dimension() const { return dimension_; }
  const Vec3& radius() const { return box_.radius(); }
  const Vec3& origin() const { return origin_; }
  Vec3 min_pt() const { return box_.min_pt(); }
  Vec3 max_pt() const { return box_.max_pt(); }
  const Vec3& center() const { return box_.center(); }
  Scalar resolution() const { return resolution_; }

  // basically equivalent to scroll_offset = (0, 0, 0)
  grid_ix_t num_cells() const { return num_cells_; }

 private:
  void calc_min_corner_ijk() {
    // set the "minimum possible" ijk, assuming the grid
    // won't stray "too far" from the initial position.
    // too far == more than 2^15 grid cells.
    // so if you voxel resolution is 1 cm, 327.68 m.
    Scalar m = -static_cast<Scalar>(std::numeric_limits<uint16_t>::max()/2)*resolution_;
    Vec3 m3(m, m, m);
    //std::cerr << "m3 = " << m3.transpose() << std::endl;
    m3 += box_.center();
    //std::cerr << "m3pcenter = " << m3.transpose() << std::endl;
    min_world_corner_ijk_ = this->world_to_grid(m3);
    //std::cerr << "min_world_corner_ijk_ = " << min_world_corner_ijk_.transpose() << std::endl;
  }


 private:
  // 3d box enclosing grid. In whatever coordinates were given (probably
  // world_view)
  ca::scrollgrid::Box<Scalar, 3> box_;

  // static origin of the grid coordinate system.
  // it's center - box.radius
  Vec3 origin_;

  // minimum world corner in ijk. used for hash
  Vec3Ix min_world_corner_ijk_;

  // number of grid cells along each axis
  Vec3Ix dimension_;

  // number of cells
  grid_ix_t num_cells_;

  // grid strides to translate from linear to 3D layout.
  // C-ordering, ie x slowest, z fastest.
  Vec3Ix strides_;

  // size of grid cells
  Scalar resolution_;

};

typedef FixedGrid3<double> FixedGrid3d;
typedef FixedGrid3<float> FixedGrid3f;

} /* ca */

#endif /* end of include guard: FIXEDGRID3_HPP_6KPVVRZF */
