#ifndef FIXEDOCCMAP2_H_FB1K65AT
#define FIXEDOCCMAP2_H_FB1K65AT

#include <memory>
#include <limits>

#include <ros/ros.h>

#include <pcl_util/point_types.hpp>

#include <scrollgrid/dense_array.hpp>
#include <scrollgrid/fixedgrid.hpp>
#include <scrollgrid/raycasting.hpp>
#include <scrollgrid/occmap_constants.hpp>

/**
 * FixedGrid + DenseArray = static occmap.
 */
namespace ca { namespace scrollgrid {

static inline
void compute_start_end_grid_ix(const Eigen::Vector2f& origin,
                               const Eigen::Vector2f& x,
                               const ca::FixedGrid2f& grid,
                               ca::Vec2Ix& start_grid_ix,
                               ca::Vec2Ix& end_grid_ix,
                               bool& hit,
                               bool& intersects) {

  // inter1 and inter2 are clipped versions of origin, x
  Eigen::Vector2f inter1(Eigen::Vector2f::Zero());
  Eigen::Vector2f inter2(Eigen::Vector2f::Zero());
  intersects = ca::clip_line2(origin, x, grid.box(), inter1, inter2);

  if (!intersects) {
    start_grid_ix.setZero();
    end_grid_ix.setZero();
    hit = false;
  } else {
    start_grid_ix = grid.world_to_grid(inter1);
    end_grid_ix = grid.world_to_grid(inter2);
    hit = grid.is_inside_box(x);
  }
}

template<class OccMapT>
void multi_update(OccMapT* occmap,
                  const Eigen::Matrix2Xf& vp,
                  const Eigen::Matrix2Xf& p) {
  ROS_ASSERT(p.rows() == vp.rows());
  for (int i=0; i < p.cols(); ++i) {
    Eigen::Vector2f xyf(p.col(i));
    Eigen::Vector2f originf(vp.col(i));
    occmap->update(xyf, originf);
  }
}


template<class T>
class FixedOccMap2 {

public:
  typedef std::shared_ptr<FixedOccMap2> Ptr;

public:
  FixedOccMap2(const FixedGrid2f& grid) {
    grid_ = grid;
    this->init();
  }

  FixedOccMap2(const FixedOccMap2& other) = delete;
  FixedOccMap2& operator=(const FixedOccMap2& other) = delete;

  void init() {
    occstats_.reset(grid_.dimension());
    occstats_.fill(OccMapConstants<T>::UNKNOWN);
  }

  void fill_unknown() {
    occstats_.fill(OccMapConstants<T>::UNKNOWN);
  }

  void update(const Eigen::Vector2f& originf, const Eigen::Vector2f& xyf) {
    bool hit = false, intersects = false;
    ca::Vec2Ix start_grid_ix(ca::Vec2Ix::Zero());
    ca::Vec2Ix end_grid_ix(ca::Vec2Ix::Zero());
    compute_start_end_grid_ix(originf,
                              xyf,
                              grid_,
                              start_grid_ix,
                              end_grid_ix,
                              hit,
                              intersects);
    if (!intersects) {
      return;
    }

    // TODO just hit at end of the ray
    Bresenham2Iterator b2itr(start_grid_ix, end_grid_ix);
    while (!b2itr.done()) {
      b2itr.step();
      ca::Vec2Ix ij(b2itr.pos());
      if (!grid_.is_inside_grid(ij)) {
        break;
      }
      mem_ix_t mem_ix = grid_.grid_to_mem(ij);
      T val = occstats_[mem_ix];
      if (b2itr.done() && hit) {
        occstats_[mem_ix] = OccMapConstants<T>::update_pos(val);
      } else {
        occstats_[mem_ix] = OccMapConstants<T>::update_neg(val);
      }
    }
  }

  void multi_update(const Eigen::Matrix2Xf& vp, const Eigen::Matrix2Xf& p) {
    multi_update<FixedOccMap2>(this, vp, p);
  }

  //TODO memory ownership?
  ca::DenseArray<T, 2> get_array() const {
    return occstats_;
  }

  virtual ~FixedOccMap2() { }

public:
  FixedGrid2f grid() const {
    return grid_;
  }

private:
  FixedGrid2f grid_;
  ca::DenseArray<T, 2> occstats_;

};

typedef FixedOccMap2<float> FixedOccMap2f;
typedef FixedOccMap2<uint8_t> FixedOccMap2u;



} }

#endif /* end of include guard: FIXEDOCCMAP2_H_FB1K65AT */
