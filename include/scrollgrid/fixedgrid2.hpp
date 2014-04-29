#ifndef FIXEDGRID2_HPP_UYCWT1KR
#define FIXEDGRID2_HPP_UYCWT1KR

#include <stdint.h>
#include <math.h>

#include <vector>
#include <algorithm>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ros/console.h>

#include <geom_cast/geom_cast.hpp>
#include <pcl_util/point_types.hpp>

#include "scrollgrid/grid_types.hpp"
#include "scrollgrid/box.hpp"

namespace ca
{

template<class CellT>
class FixedGrid2 {
public:
  typedef CellT CellType;
  typedef Eigen::Vector2d Vec2d;

  typedef CellT * ArrayType;
  typedef CellT * iterator;
  typedef const CellT * const_iterator;

  typedef boost::shared_ptr<FixedGrid2> Ptr;
  typedef boost::shared_ptr<const FixedGrid2> ConstPtr;

public:
  FixedGrid2() :
      box_(),
      origin_(0, 0),
      dimension_(0, 0),
      num_cells_(0),
      stride_(0),
      resolution_(0),
      grid_(NULL),
      begin_(NULL),
      end_(NULL)
  { }

  FixedGrid2(const Vec2d& center,
             const Vec2Ix& dimension,
             double resolution) :
      box_(center-(dimension.cast<double>()*resolution)/2,
           center+(dimension.cast<double>()*resolution)/2),
      origin_(center-box_.radius()),
      dimension_(dimension),
      num_cells_(dimension.prod()),
      stride_(dimension[1]),
      resolution_(resolution),
      grid_(new CellT[num_cells_]),
      begin_(&grid_[0]),
      end_(&grid_[0]+num_cells_)
  { }

  virtual ~FixedGrid2() {
    if (grid_ != NULL) { delete[] grid_; }
    grid_ = begin_ = end_ = NULL;
  }

private:
  FixedGrid2(const FixedGrid2& other);
  FixedGrid2& operator=(const FixedGrid2& other);

public:

  void reset(const Vec2d& center,
             const Vec2Ix& dimension,
             double resolution) {
    box_.set_center(center);
    box_.set_radius((dimension.cast<double>()*resolution)/2);
    origin_ = center - box_.radius();
    dimension_ = dimension;
    num_cells_ = dimension.prod();
    stride_ = dimension[1];
    resolution_ = resolution;
    if (grid_ != NULL) { delete[] grid_; }
    grid_ = new CellT[num_cells_];
    begin_ = &grid_[0];
    end_ = &grid_[0]+num_cells_;
  }

  void set_zero() {
    CellT empty_cell;
    std::fill(&grid_[0], &grid_[0]+num_cells_, empty_cell);
  }

  void copy_from(ca::FixedGrid2<CellT>& other) {
    this->reset(other.center(),
                other.dimension(),
                other.resolution());
    std::copy(other.begin(), other.end(), &grid_[0]);
  }

  /**
   *  purely memory-based view of the data.
   */
  iterator begin() { return begin_; }
  const_iterator cbegin() const { return begin_; }

  iterator end() { return end_; }
  const_iterator cend() { return end_; }

  /**
   * Initialize with same physical dimensions as other.
   * Does not copy contents of other.
   */
  void init_like(const FixedGrid2& other) {
    this->reset(other.center(),
                other.dimension(),
                other.resolution());
  }

  /**
   * Is inside 3D box containing grid?
   * pt is in same frame as center (probably world_view)
   */
  bool is_inside_box(const Vec2d& pt) const {
    return box_.contains(pt);
  }

  template<class PointT>
  bool is_inside_box(const PointT& pt) const {
    return box_.contains(ca::point_cast<Eigen::Vector2d>(pt));
  }

  /**
   * is i, j inside the grid limits?
   */
  bool is_inside_grid(const Vec2Ix& grid_ix) const {
    return ((grid_ix.array() >= 0).all() &&
            (grid_ix.array() < dimension_.array()).all());
  }

  bool is_inside_grid(grid_ix_t i, grid_ix_t j) const {
    return this->is_inside_grid(Vec2Ix(i, j));
  }

  /**
   * Given position in world coordinates, return grid coordinates.
   * Note: does not check if point is inside grid.
   */
  Vec2Ix world_to_grid(const Vec2d& xyz) const {
    Vec2d tmp = ((xyz - origin_).array() - 0.5*resolution_)/resolution_;
    //return tmp.cast<grid_ix_t>();
    return Vec2Ix(round(tmp.x()), round(tmp.y()), round(tmp.z()));
  }

  Vec2Ix world_to_grid(double x, double y, double z) const {
    return this->world_to_grid(Vec2d(x, y, z));
  }

  Vec2d grid_to_world(const Vec2Ix& grid_ix) const {
    Vec2d w((grid_ix.cast<double>()*resolution_ + origin_).array() + 0.5*resolution_);
    return w;
  }

  Vec2d grid_to_world(grid_ix_t i, grid_ix_t j) const {
    return this->grid_to_world(Vec2Ix(i, j));
  }

  grid_ix_t grid_to_mem(grid_ix_t i, grid_ix_t j) const {
    return this->grid_to_mem(Vec2Ix(i, j));
  }

  /**
   * Note that this wraps the z dimension.
   */
  grid_ix_t grid_to_mem(const Vec2Ix& grid_ix) const {
    return stride_*grid_ix[0] + grid_ix[1];
  }

  Vec2Ix mem_to_grid(grid_ix_t mem_ix) const {
    grid_ix_t i = mem_ix/stride_;
    mem_ix -= i*stride_;
    grid_ix_t j = mem_ix;
    mem_ix -= j;

    return Vec2Ix(i, j);
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

 public:
  grid_ix_t dim_i() const { return dimension_[0]; }
  grid_ix_t dim_j() const { return dimension_[1]; }
  const Vec2Ix& dimension() const { return dimension_; }
  const Vec2d& radius() const { return box_.radius(); }
  const Vec2d& origin() const { return origin_; }
  Vec2d min_pt() const { return box_.min_pt(); }
  Vec2d max_pt() const { return box_.max_pt(); }
  const Vec2d& center() const { return box_.center(); }
  double resolution() const { return resolution_; }

  grid_ix_t num_cells() const { return num_cells_; }

 private:
  // 3d box enclosing grid. In whatever coordinates were given (probably
  // world_view)
  ca::Box<double, 2> box_;

  // static origin of the grid coordinate system.
  Vec2d origin_;

  // number of grid cells along each axis
  Vec2Ix dimension_;

  // number of cells
  grid_ix_t num_cells_;

  // grid strides to translate from linear to 3D layout.
  // C-ordering, ie x slowest, z fastest.
  grid_ix_t stride_;

  // size of grid cells
  double resolution_;

  ArrayType grid_;
  iterator begin_, end_;

};

}

#endif /* end of include guard: FIXEDGRID2_HPP_UYCWT1KR */
