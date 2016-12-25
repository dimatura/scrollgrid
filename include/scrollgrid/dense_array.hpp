/**
 * Copyright (c) 2015 Carnegie Mellon University, Daniel Maturana <dimatura@cmu.edu>
 *
 * For License information please see the LICENSE file in the root directory.
 *
 */

#ifndef DENSE_ARRAY_HPP_KYZBZH9N
#define DENSE_ARRAY_HPP_KYZBZH9N

#include <cmath>
#include <cstdint>

#include <stdexcept>
#include <memory>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ros/console.h>

#include <pcl_util/point_types.hpp>

#include "scrollgrid/grid_types.hpp"

namespace ca {

/**
 * A dense ND array.
 * Maps ijk to an CellT.
 * No notion of origin, scrolling etc.
 * The *GridN classes handle that.
 * TODO: maybe use std::shared_array
 */
template<class CellT, int Dim_>
class DenseArray {
public:
  enum { Dim = Dim_ };
  typedef CellT CellType;
  typedef CellT * ArrayType;
  typedef CellT * iterator;
  typedef const CellT * const_iterator;
  typedef Eigen::Matrix<grid_ix_t, Dim, 1> VecIx;

  typedef std::shared_ptr<DenseArray> Ptr;
  typedef std::shared_ptr<const DenseArray> ConstPtr;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

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

  DenseArray() :
      dimension_(VecIx::Zero()),
      num_cells_(0),
      strides_(VecIx::Zero()),
      grid_(nullptr),
      begin_(nullptr),
      end_(nullptr),
      owns_memory_(false)
  { }

  DenseArray(const VecIx& dimension) :
      dimension_(dimension),
      num_cells_(dimension.prod()),
      grid_(new CellT[num_cells_]()),
      begin_(&grid_[0]),
      end_(&grid_[0]+num_cells_),
      owns_memory_(true)
  {
    for (int i=0; i < Dim; ++i) { strides_(Dim-i-1) = dimension_.tail(i).prod(); }
  }

  /**
   * Used when wrapping an external chunk of data.
   */
  DenseArray(const VecIx& dimension, ArrayType grid_data) :
      dimension_(dimension),
      num_cells_(dimension.prod()),
      grid_(grid_data),
      begin_(&grid_[0]),
      end_(&grid_[0]+num_cells_),
      owns_memory_(false)
  {
    for (int i=0; i < Dim; ++i) { strides_(Dim-i-1) = dimension_.tail(i).prod(); }
  }

  void CopyFrom(const DenseArray& other) {
    this->reset(other.dimension());
    std::copy(other.begin(), other.end(), this->begin());
  }

  virtual ~DenseArray() {
    if (owns_memory_ && grid_) {
      delete[] grid_;
      grid_ = nullptr;
    }
  }

  void reset(const VecIx& dimension) {
    dimension_ = dimension;
    num_cells_ = dimension.prod();
    // TODO configurable
    for (int i=0; i < Dim; ++i) { strides_(Dim-i-1) = dimension_.tail(i).prod(); }
    if (owns_memory_ && grid_) { delete[] grid_; }
    // note new T[n]() != new T[n];
    // the former initializes to zero
    grid_ = new CellT[num_cells_]();
    begin_ = &(grid_[0]);
    end_ = &(grid_[0])+num_cells_;
  }

  void reset(const VecIx& dimension, ArrayType grid_data) {
    if (owns_memory_ && grid_) { delete[] grid_; }
    dimension_ = dimension;
    num_cells_ = dimension.prod();
    // TODO configurable
    for (int i=0; i < Dim; ++i) { strides_(Dim-i-1) = dimension_.tail(i).prod(); }
    owns_memory_ = false;
    grid_ = grid_data;
    begin_ = &(grid_[0]);
    end_ = &(grid_[0])+num_cells_;
  }

  size_t allocated_bytes() {
    return sizeof(CellT)*num_cells_;
  }

public:

  void fill(const CellT& val) {
    std::fill(this->begin(), this->end(), val);
  }

public:
  iterator begin() const { return begin_; }
  iterator end() const { return end_; }

public:

  /**
   * NOTE be careful when using these!
   * They do no take into account any sort of wrapping
   * coming from scrollgrid.
   * CAREFUL.
   */
  grid_ix_t local_grid_to_mem(const VecIx& grid_ix) const {
    return strides_.dot(grid_ix);
  }

  VecIx local_mem_to_grid(grid_ix_t mem_ix) const {
    VecIx out(VecIx::Zero());
    for (int i=0; i < Dim; ++i) {
      grid_ix_t ix = mem_ix/strides_[i];
      mem_ix -= ix*strides_[i];
      out(i) = ix;
    }
    return out;
  }

public:

  CellType& get(const VecIx& grid_ix) {
    grid_ix_t mem_ix = this->local_grid_to_mem(grid_ix);
    return grid_[mem_ix];
  }

  const CellType& get(const VecIx& grid_ix) const {
    grid_ix_t mem_ix = this->local_grid_to_mem(grid_ix);
    return grid_[mem_ix];
  }

public:

  /**
   * Bound check
   */
  CellType& get_safe(grid_ix_t mem_ix) {
    if (mem_ix < 0 || mem_ix >= num_cells_) {
      throw std::out_of_range("bad mem_ix");
    }
    return grid_[mem_ix];
  }

  /**
   * Bound check
   */
  const CellType& get_safe(grid_ix_t mem_ix) const {
    if (mem_ix < 0 || mem_ix >= num_cells_) {
      throw std::out_of_range("bad mem_ix");
    }
    return grid_[mem_ix];
  }

public:

  /**
   * No bound check
   */
  CellType& operator[](grid_ix_t mem_ix) {
    return grid_[mem_ix];
  }

  /**
   * No bound check
   */
  const CellType& operator[](grid_ix_t mem_ix) const {
    return grid_[mem_ix];
  }

public:
  // properties
  grid_ix_t dim(int i) const { return dimension_[i]; }
  VecIx dimension() const { return dimension_; }
  grid_ix_t num_cells() const { return num_cells_; }
  ArrayType data() const { return &grid_[0]; }
  grid_ix_t stride(int i) { return strides_[i]; }
  VecIx strides() const { return strides_; }

private:
  // number of grid cells along each axis
  VecIx dimension_;

  // number of cells
  grid_ix_t num_cells_;

  // grid strides to translate from linear to 3D layout.
  // C-ordering, ie x slowest, z fastest.
  VecIx strides_;

  ArrayType grid_;
  iterator begin_, end_;

  // does this object own the grid_ mem
  bool owns_memory_;

//private:
 // DenseArray3(const DenseArray3& other);
 // DenseArray3& operator=(const DenseArray3& other);
};

typedef DenseArray<float, 2> DenseArray2f;
typedef DenseArray<float, 3> DenseArray3f;

}/* ca */

#endif /* end of include guard: DENSE_ARRAY_HPP_KYZBZH9N */
