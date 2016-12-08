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

  bool compute_start_end_grid_ix(const Eigen::Vector2f& xyf,
                                 const Eigen::Vector2f& originf,
                                 ca::Vec2Ix& start_grid_ix,
                                 ca::Vec2Ix& end_grid_ix,
                                 bool& hit) {
    // a stupid little dance
    const auto& box(grid_.box());
    Ray2<float> ray(originf, (xyf-originf));
    if (!ca::aabb_ray_intersect(box, ray)) {
      // no intersection
      return false;
    }
    Eigen::Vector2f start_ray(ray.point_at(ray.tmin+1e-6));
    Eigen::Vector2f end_ray(ray.point_at(ray.tmax-1e-6));
    if ((xyf-originf).norm() < (start_ray-originf).norm()) {
      // ray ends before actually hitting box
      return false;
    }
    // adjust segment endpoints
    if (grid_.is_inside_box(originf)) {
      start_ray = originf;
    }

    hit = false;
    if (grid_.is_inside_box(xyf)) {
      hit = true;
      end_ray = xyf;
    }

    start_grid_ix = grid_.world_to_grid(start_ray);
    end_grid_ix = grid_.world_to_grid(end_ray);

    if (!grid_.is_inside_grid(end_grid_ix) ||
        !grid_.is_inside_grid(start_grid_ix)) {
      return false;
    }
    return true;
  }

  void update(const Eigen::Vector2f& xyf, const Eigen::Vector2f& originf) {

    bool hit = false;
    ca::Vec2Ix start_grid_ix;
    ca::Vec2Ix end_grid_ix;
    bool intersects = compute_start_end_grid_ix(xyf, originf, start_grid_ix, end_grid_ix, hit);
    if (!intersects) {
      return;
    }

    // TODO just hit at end of the ray
    for (Bresenham2Iterator b2itr(start_grid_ix, end_grid_ix);
         !b2itr.done();
         b2itr.step()) {
      Vec2Ix ij(b2itr.pos());
      std::cout << "ij = " << ij << std::endl;
      mem_ix_t mem_ix = grid_.grid_to_mem(ij);
      T cur_val(occstats_[mem_ix]);
      if (hit && b2itr.done()) {
        occstats_[mem_ix] = OccMapConstants<T>::update_pos(cur_val);
      } else {
        occstats_[mem_ix] = OccMapConstants<T>::update_neg(cur_val);
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
