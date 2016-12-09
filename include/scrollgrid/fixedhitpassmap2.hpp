#ifndef FIXEDHITPASSMAP2_HPP_FZQPVRHM
#define FIXEDHITPASSMAP2_HPP_FZQPVRHM

#include <memory>
#include <limits>

#include <ros/ros.h>

#include <pcl_util/point_types.hpp>

#include <scrollgrid/dense_array.hpp>
#include <scrollgrid/fixedgrid.hpp>
#include <scrollgrid/raycasting.hpp>

/**
 * FixedGrid + DenseArray = static occmap.
 */
namespace ca { namespace scrollgrid {

template<class T>
class FixedHitPassMap2 {

public:
  typedef std::shared_ptr<FixedHitPassMap2> Ptr;

public:
  FixedHitPassMap2(const FixedGrid2f& grid) {
    grid_ = grid;
    this->init();
  }

  FixedHitPassMap2(const FixedHitPassMap2& other) = delete;
  FixedHitPassMap2& operator=(const FixedHitPassMap2& other) = delete;

  void init() {
    hits_.reset(grid_.dimension());
    passes_.reset(grid_.dimension());
    hits_.fill(0);
    passes_.fill(0);
  }

  void fill_unknown() {
    hits_.fill(0);
    passes_.fill(0);
  }

  void update(const Eigen::Vector2f& xyf, const Eigen::Vector2f& originf) {

    bool hit = false, intersects = false;
    ca::Vec2Ix start_grid_ix(ca::Vec2Ix::Zero());
    ca::Vec2Ix end_grid_ix(ca::Vec2Ix::Zero());
    this->compute_start_end_grid_ix(xyf, originf, start_grid_ix, end_grid_ix, hit, intersects);

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
      if (b2itr.done() && hit) {
        //std::cerr << ij.transpose() << " hit " << std::endl;
        hits_[mem_ix] += 1;
      } else {
        //std::cerr << ij.transpose() << " pass " << std::endl;
        passes_[mem_ix] += 1;
      }
    }
  }

  void multi_update(const Eigen::Matrix2Xf& p, const Eigen::Matrix2Xf& vp) {
    ROS_ASSERT(p.rows() == vp.rows());
    const auto& box(grid_.box());
    for (int i=0; i < p.cols(); ++i) {
      Eigen::Vector2f xyf(p.col(i));
      Eigen::Vector2f originf(vp.col(i));
      this->update(xyf, originf);
    }
  }

  //TODO memory ownership?
  ca::DenseArray<T, 2> get_hits_array() const {
    return hits_;
  }

  ca::DenseArray<T, 2> get_passes_array() const {
    return passes_;
  }

  virtual ~FixedHitPassMap2() { }

public:
  FixedGrid2f grid() const {
    return grid_;
  }

private:

  void compute_start_end_grid_ix(const Eigen::Vector2f& x,
                                 const Eigen::Vector2f& origin,
                                 ca::Vec2Ix& start_grid_ix,
                                 ca::Vec2Ix& end_grid_ix,
                                 bool& hit,
                                 bool& intersects) {

    // inter1 and inter2 are clipped versions of origin, x
    Eigen::Vector2f inter1(Eigen::Vector2f::Zero());
    Eigen::Vector2f inter2(Eigen::Vector2f::Zero());
    intersects = ca::clip_line2(origin, x, grid_.box(), inter1, inter2);

    if (!intersects) {
      start_grid_ix.setZero();
      end_grid_ix.setZero();
      hit = false;
    } else {
      start_grid_ix = grid_.world_to_grid(inter1);
      end_grid_ix = grid_.world_to_grid(inter2);
      hit = grid_.is_inside_box(x);
    }
  }


private:
  FixedGrid2f grid_;
  ca::DenseArray<T, 2> hits_;
  ca::DenseArray<T, 2> passes_;

};

typedef FixedHitPassMap2<int> FixedHitPassMap2i;
typedef FixedHitPassMap2<uint8_t> FixedHitPassMap2u;

} }

#endif /* end of include guard: FIXEDHITPASSMAP2_HPP_FZQPVRHM */
