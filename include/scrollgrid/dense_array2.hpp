/**
 * Copyright (c) 2015 Carnegie Mellon University, Daniel Maturana <dimatura@cmu.edu>
 *
 * For License information please see the LICENSE file in the root directory.
 *
 */


#ifndef DENSE_ARRAY2_HPP_ZJGDW1JR
#define DENSE_ARRAY2_HPP_ZJGDW1JR

#include <cmath>
#include <cstdint>

#include <vector>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ros/console.h>

#include <geom_cast/geom_cast.hpp>
#include <pcl_util/point_types.hpp>

#include "scrollgrid/grid_types.hpp"

namespace ca
{

/**
 * A dense 2D array.
 * Maps ij to an CellT.
 * No notion of origin, scrolling etc.
 * The *Grid2 classes handle that.
 */
template<class CellT>
class DenseArray2 {
public:
  typedef CellT CellType;
  typedef CellT * ArrayType;
  typedef CellT * iterator;
  typedef const CellT * const_iterator;

  typedef std::shared_ptr<DenseArray2> Ptr;
  typedef std::shared_ptr<const DenseArray2> ConstPtr;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  DenseArray2() :
      dimension_(0, 0),
      num_cells_(0),
      strides_(0, 0),
      grid_(nullptr),
      begin_(nullptr),
      end_(nullptr),
      own_memory_(false)
  { }

  DenseArray2(const Vec2Ix& dimension) :
      dimension_(dimension),
      num_cells_(dimension.prod()),
      strides_(dimension[1], 1),
      grid_(new CellT[num_cells_]),
      begin_(&grid_[0]),
      end_(&grid_[0]+num_cells_),
      own_memory_(true)
  { }

  virtual ~DenseArray2() {
    if (grid_ && own_memory_) { delete[] grid_; }
  }

  void reset(const ca::Vec2Ix& dimension) {
    if (grid_) { delete[] grid_; }
    dimension_ = dimension;
    num_cells_ = dimension.prod();
    strides_ = ca::Vec2Ix(dimension[1], 1);
    grid_ = new CellT[num_cells_]();
    begin_ = &grid_[0];
    end_ = &grid_[0]+num_cells_;
  }

  void reset(const Vec2Ix& dimension, ArrayType grid_data) {
    if (grid_ != nullptr) {
      delete[] grid_;
    }
    dimension_ = dimension;
    num_cells_ = dimension.prod();
    // TODO configurable strides
    strides_ = ca::Vec2Ix(dimension[1], 1);
    grid_ = grid_data;
    begin_ = &(grid_[0]);
    end_ = &(grid_[0]) + num_cells_;
    own_memory_ = false;
  }

public:

  /**
   * NOTE be careful when using these!
   * They do no take into account any sort of wrapping
   * coming from scrollgrid.
   */
  grid_ix_t grid_to_mem(grid_ix_t i, grid_ix_t j) const {
    return this->grid_to_mem(Vec2Ix(i, j));
  }

  grid_ix_t grid_to_mem(const Vec2Ix& grid_ix) const {
    return strides_.dot(grid_ix);
  }

  /**
   * Note: no bound checking.
   */
  CellType& get(grid_ix_t i, grid_ix_t j) {
    grid_ix_t mem_ix = this->grid_to_mem(i, j);
    return grid_[mem_ix];
  }

  const CellType& get(grid_ix_t i, grid_ix_t j) const {
    grid_ix_t mem_ix = this->grid_to_mem(i, j);
    return grid_[mem_ix];
  }

  CellType& get(const Vec2Ix& grid_ix) {
    grid_ix_t mem_ix = this->grid_to_mem(grid_ix);
    return grid_[mem_ix];
  }

  const CellType& get(const Vec2Ix& grid_ix) const {
    grid_ix_t mem_ix = this->grid_to_mem(grid_ix);
    return grid_[mem_ix];
  }

  CellType& get(grid_ix_t mem_ix) {
    return grid_[mem_ix];
  }

  const CellType& get(grid_ix_t mem_ix) const {
    return grid_[mem_ix];
  }

  CellType& operator[](mem_ix_t mem_ix) {
    return grid_[mem_ix];
  }

  const CellType& operator[](mem_ix_t mem_ix) const {
    return grid_[mem_ix];
  }


public:
  // properties
  grid_ix_t dim_i() const { return dimension_[0]; }
  grid_ix_t dim_j() const { return dimension_[1]; }
  Vec2Ix dimension() const { return dimension_; }
  grid_ix_t num_cells() const { return num_cells_; }
  ArrayType data() const { return &grid_[0]; }
  grid_ix_t stride(int i) { return strides_[i]; }
  Vec2Ix strides() const { return strides_; }

private:
  DenseArray2(const DenseArray2& other);
  DenseArray2& operator=(const DenseArray2& other);

private:
  bool own_memory_;
  // number of grid cells along each axis
  Vec2Ix dimension_;

  // number of cells
  grid_ix_t num_cells_;

  // grid strides to translate from linear to 2D layout.
  // C-ordering, ie x slowest, z fastest.
  Vec2Ix strides_;

  ArrayType grid_;
  iterator begin_, end_;
};

typedef DenseArray2<float> DenseArray2f;
typedef DenseArray2<double> DenseArray2d;
typedef DenseArray2<uint8_t> DenseArray2u;

} /* ca */

#endif /* end of include guard: DENSE_ARRAY2_HPP_ZJGDW1JR */
