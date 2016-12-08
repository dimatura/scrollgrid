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

  bool compute_start_end_grid_ix(const Eigen::Vector2f& x,
                                 const Eigen::Vector2f& origin,
                                 ca::Vec2Ix& start_grid_ix,
                                 ca::Vec2Ix& end_grid_ix,
                                 bool& hit) {
    // a stupid little dance
    const auto& box(grid_.box());
    Ray2<float> ray(origin, (x-origin));
    if (!ca::aabb_ray_intersect(box, ray)) {
      // no intersection at all - skip
      return false;
    }
    // position of ray at start of intersection with box
    Eigen::Vector2f inter_start(ray.point_at(ray.tmin+1e-6));
    // position of ray at end of intersection with box
    Eigen::Vector2f inter_end(ray.point_at(ray.tmax-1e-6));

    Eigen::Vector2f ray_start = Eigen::Vector2f::Zero();
    Eigen::Vector2f ray_end = Eigen::Vector2f::Zero();

    // assuming an intersection, four cases
    // 1. x inside, origin inside
    // 2. x inside, origin outside
    // 3. x outside, origin inside
    // 4. x outside, origin outside
    //   4a. ray points away from box, segment doesn't touch it
    //   4b. ray points towards box, segment doesn't touch it
    //   4c. ray points towards box, segment crosses it

    float range = (x-origin).norm();

    bool x_inside = grid_.is_inside_box(x);
    bool origin_inside = grid_.is_inside_box(origin);

    if (x_inside && origin_inside) {
      std::cout << "case 1\n";
      // case 1; use whole beam
      ray_start = origin;
      ray_end = x;
    } else if (x_inside && !origin_inside) {
      std::cout << "case 2\n";
      // case 2; start at intersection, same origin
      ray_start = inter_start;
      ray_end = x;
    } else if (!x_inside && origin_inside) {
      std::cout << "case 3\n";
      // case 3
      ray_start = origin;
      ray_end = inter_end;
      std::cerr << "ray.tmin = " << ray.tmin << std::endl;
      std::cerr << "ray.tmax = " << ray.tmax << std::endl;
      std::cerr << "inter_start = " << inter_start.transpose() << std::endl;
      std::cerr << "inter_end = " << inter_end.transpose() << std::endl;
      std::cerr << "ray_start = " << ray_start.transpose() << std::endl;
      std::cerr << "ray_start = " << ray_end.transpose() << std::endl;
    } else if (!x_inside && !origin_inside) {
      if (ray.tmin < 0) {
        std::cout << "case 4a\n";
        // case 4a
        return false;
      }
      if (range < ray.tmin) {
        std::cout << "case 4b\n";
        // case 4b
        return false;
      }
      std::cout << "case 4c\n";
      // case 4c
      ray_start = inter_start;
      ray_end = inter_end;
    }

    hit = x_inside;
    start_grid_ix = grid_.world_to_grid(ray_start);
    end_grid_ix = grid_.world_to_grid(ray_end);

    return true;
  }

  void update(const Eigen::Vector2f& xyf, const Eigen::Vector2f& originf) {

    bool hit = false;
    ca::Vec2Ix start_grid_ix;
    ca::Vec2Ix end_grid_ix;
    bool intersects = compute_start_end_grid_ix(xyf, originf, start_grid_ix, end_grid_ix, hit);

    std::cerr << "intersects = " << intersects << std::endl;
    std::cerr << "start_grid_ix = " << start_grid_ix.transpose() << std::endl;
    std::cerr << "end_grid_ix = " << end_grid_ix.transpose() << std::endl;

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
        std::cerr << ij.transpose() << " hit " << std::endl;
        hits_[mem_ix] += 1;
      } else {
        std::cerr << ij.transpose() << " pass " << std::endl;
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
  FixedGrid2f grid_;
  ca::DenseArray<T, 2> hits_;
  ca::DenseArray<T, 2> passes_;

};

typedef FixedHitPassMap2<int> FixedHitPassMap2i;
typedef FixedHitPassMap2<uint8_t> FixedHitPassMap2u;

} }

#endif /* end of include guard: FIXEDHITPASSMAP2_HPP_FZQPVRHM */
