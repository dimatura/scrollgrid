#ifndef DENSE_ARRAY3_HPP_JEO7CAXQ
#define DENSE_ARRAY3_HPP_JEO7CAXQ

#include <math.h>
#include <stdint.h>

#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ros/console.h>

#include <geom_cast/geom_cast.hpp>
#include <pcl_util/point_types.hpp>

#include "pc_semantic_micro/geom_util.hpp"

namespace ca
{

/**
 * A dense 3D array.
 * Maps ijk to an CellT.
 * No notion of origin, scrolling etc.
 * The *Grid3 classes handle that.
 *
 */
template<class CellT>
class DenseArray3 {
public:
  typedef CellT CellType;
  typedef CellT * ArrayType;
  typedef CellT * iterator;
  typedef const CellT * const_iterator;

  typedef boost::shared_ptr<DenseArray3> Ptr;
  typedef boost::shared_ptr<const DenseArray3> ConstPtr;

public:
  // TODO perhaps stride/memory layout should
  // be further configurable.
  // TODO maybe *grid3 should just map to ijk and then
  // this class will map to linear memory address.
  // problem: then it needs scrolling info
  // I guess we can separate
  // - world_to_grid (infinite grid)
  // - grid_to_storage (linear mem or hash key)
  // - storage (dense array or sparse table)
  //

  DenseArray3() :
      dimension_(0, 0, 0),
      num_cells_(0),
      strides_(0, 0, 0),
      grid_(NULL),
      begin_(NULL),
      end_(NULL)
  { }

  DenseArray3(const Vec3Ix& dimension) :
      dimension_(dimension),
      num_cells_(dimension.prod()),
      strides_(dimension.tail<2>().prod(), dimension[2], 1),
      grid_(new CellT[num_cells_]),
      begin_(&grid_[0]),
      end_(&grid_[0]+num_cells_)
  { }

  virtual ~DenseArray3() {
    if (grid_) { delete[] grid_; }
  }

  void reset(const Vec3Ix& dimension) {
    dimension_ = dimension;
    num_cells_ = dimension.prod();
    // TODO configurable
    strides_ = Vec3Ix(dimension.tail<2>().prod(), dimension[2], 1);
    if (grid_) { delete[] grid_; }
    // note new T[n]() != new T[n];
    // the former initializes to zero
    grid_ = new CellT[num_cells_]();
    begin_ = &(grid_[0]);
    end_ = &(grid_[0])+num_cells_;
  }

public:

  /**
   * these are for convenience, prbly shouldnt be used
   */
  grid_ix_t grid_to_mem(grid_ix_t i, grid_ix_t j, grid_ix_t k) const {
    return this->grid_to_mem(Vec3Ix(i, j, k));
  }

  grid_ix_t grid_to_mem(const Vec3Ix& grid_ix) const {
    return strides_.dot(grid_ix);
  }

  Vec3Ix mem_to_grid_local(grid_ix_t mem_ix) const {
    grid_ix_t i = mem_ix/strides_[0];
    mem_ix -= i*strides_[0];
    grid_ix_t j = mem_ix/strides_[1];
    mem_ix -= j*strides_[1];
    grid_ix_t k = mem_ix;

    return Vec3Ix(i, j, k);
  }

public:

  CellType& get(grid_ix_t i, grid_ix_t j, grid_ix_t k) {
    grid_ix_t mem_ix = this->grid_to_mem(i, j, k);
    return grid_[mem_ix];
  }

  const CellType& get(grid_ix_t i, grid_ix_t j, grid_ix_t k) const {
    grid_ix_t mem_ix = this->grid_to_mem(i, j, k);
    return grid_[mem_ix];
  }

  CellType& get(const Vec3Ix& grid_ix) {
    grid_ix_t mem_ix = this->grid_to_mem(grid_ix);
    return grid_[mem_ix];
  }

  const CellType& get(const Vec3Ix& grid_ix) const {
    grid_ix_t mem_ix = this->grid_to_mem(grid_ix);
    return grid_[mem_ix];
  }

  CellType& get(grid_ix_t mem_ix) {
    return grid_[mem_ix];
  }

  const CellType& get(grid_ix_t mem_ix) const {
    return grid_[mem_ix];
  }

public:

  CellType& operator[](grid_ix_t mem_ix) {
    return grid_[mem_ix];
  }

  const CellType& operator[](grid_ix_t mem_ix) const {
    return grid_[mem_ix];
  }

public:
  // properties
  grid_ix_t dim_i() const { return dimension_[0]; }
  grid_ix_t dim_j() const { return dimension_[1]; }
  grid_ix_t dim_k() const { return dimension_[2]; }
  Vec3Ix dimension() const { return dimension_; }
  grid_ix_t num_cells() const { return num_cells_; }

private:
  DenseArray3(const DenseArray3& other);
  DenseArray3& operator=(const DenseArray3& other);

private:
  // number of grid cells along each axis
  Vec3Ix dimension_;

  // number of cells
  grid_ix_t num_cells_;

  // grid strides to translate from linear to 3D layout.
  // C-ordering, ie x slowest, z fastest.
  Vec3Ix strides_;

  ArrayType grid_;
  iterator begin_, end_;
};

} /* ca */

#endif /* end of include guard: DENSE_ARRAY3_HPP_JEO7CAXQ */
