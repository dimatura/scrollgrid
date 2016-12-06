#ifndef FIXEDOCCMAP3_H_J5TFK2HQ
#define FIXEDOCCMAP3_H_J5TFK2HQ

#include <memory>
#include <limits>

#include <ros/ros.h>

#include <pcl_util/point_types.hpp>
#include <scrollgrid/dense_array3.hpp>
#include <scrollgrid/fixedgrid3.hpp>
#include <scrollgrid/raycasting.hpp>

/**
 * FixedGrid3 + DenseArray3 = static occmap.
 */
namespace ca { namespace scrollgrid {

template<class T>
struct OccMapConstants;

template<>
struct OccMapConstants<float> {
  static constexpr float UNKNOWN = 0.5f;
  static constexpr float OCCUPIED = 0.9f;
  static constexpr float FREE = 0.0f;
  static constexpr float UPDATE_POS = 0.08f;
  static constexpr float UPDATE_NEG = 0.01f;

  float update_pos(float x) {
    return (std::min(x + UPDATE_POS, OCCUPIED));
  }

  float update_neg(float x) {
    return (std::max(x - UPDATE_NEG, FREE));
  }
};

template<>
struct OccMapConstants<uint8_t> {
  static constexpr uint8_t UNKNOWN = 128;
  static constexpr uint8_t OCCUPIED = 250;
  static constexpr uint8_t FREE = 0;
  static constexpr int32_t UPDATE_POS = 20;
  static constexpr int32_t UPDATE_NEG = 2;

  uint8_t update_pos(uint8_t x) {
    // TODO avoid casts
    int32_t new_val = static_cast<int32_t>(x) + UPDATE_POS;
    new_val = std::min(static_cast<int32_t>(OCCUPIED), new_val);
    return static_cast<uint8_t>(new_val);
  }

  uint8_t update_neg(uint8_t x) {
    // TODO avoid casts
    int32_t new_val = static_cast<int32_t>(x) - UPDATE_NEG;
    new_val = std::max(static_cast<int32_t>(FREE), new_val);
    return static_cast<uint8_t>(new_val);
  }

};


template<class T>
class FixedOccMap3 {

public:
  typedef std::shared_ptr<FixedOccMap3> Ptr;

private:
  ca::FixedGrid3f grid_;
  ca::DenseArray3<T> occstats_;


public:
  FixedOccMap3(const ca::FixedGrid3f& grid) {
    grid_ = grid;
    this->Init();
  }

  FixedOccMap3(const FixedOccMap3& other) = delete;
  FixedOccMap3& operator=(const FixedOccMap3& other) = delete;

  void Init() {
    occstats_.reset(grid_.dimension());
    occstats_.fill(static_cast<uint8_t>(OccMapConstants<T>::UNKNOWN));
  }

  void SetUnknown() {
    occstats_.fill(static_cast<uint8_t>(OccMapConstants<T>::UNKNOWN));
  }

  void Update(RowMatrixX3f& p, RowMatrixX3f& vp) {
    //UpdateOccupancy(xyzwvp, grid_, occstats_);

    ROS_ASSERT(p.rows() == vp.cols());

    const auto& box(grid_.box());
    for (int i=0; i < p.rows(); ++i) {
      Eigen::Vector3f xyzf(p.row(i));
      Eigen::Vector3f originf(vp.row(i));
      Ray3<float> ray(originf, (xyzf-originf));
      if (!ca::aabb_ray_intersect(box, ray)) {
        //std::cerr << "no box intersect" << std::endl;
        continue;
      }
      Eigen::Vector3f start_ray = ray.point_at(ray.tmin+1e-6);
      Eigen::Vector3f end_ray = ray.point_at(ray.tmax-1e-6);
      if ((xyzf-originf).norm() < (start_ray-originf).norm()) {
        // ray ends before actually hitting box
        continue;
      }
      // adjust segment endpoints
      if (grid_.is_inside_box(originf)) {
        start_ray = originf;
      }

      bool hit = false;
      if (grid_.is_inside_box(xyzf)) {
        hit = true;
        end_ray = xyzf;
      }

      ca::Vec3Ix start_grid_ix(grid_.world_to_grid(start_ray));
      ca::Vec3Ix end_grid_ix(grid_.world_to_grid(end_ray));

      if (!grid_.is_inside_grid(end_grid_ix) ||
          !grid_.is_inside_grid(start_grid_ix)) {
        continue;
      }

      // TODO just hit at end of the ray
      for (Bresenham3Iterator b3itr(start_grid_ix, end_grid_ix);
           !b3itr.done();
           b3itr.step()) {
        Vec3Ix ijk(b3itr.pos());
        mem_ix_t mem_ix = grid_.grid_to_mem(ijk);
        T cur_val(occstats_[mem_ix]);
        if (hit && b3itr.done()) {
          occstats_[mem_ix] = OccMapConstants<T>::update_pos(cur_val);
        } else {
          occstats_[mem_ix] = OccMapConstants<T>::update_neg(cur_val);
        }
      }
    }
  }

  void CopyToBuffer(T * buffer) const {
    // TODO dangerous!
    // perhaps faster but even more dangerous, is to directly
    // wrap a numpy array buffer
    std::copy(occstats_.begin(), occstats_.end(), buffer);
  }

  virtual ~FixedOccMap3() { }

public:
  ca::FixedGrid3f grid() const {
    return grid_;
  }

};

} }

#endif /* end of include guard: FIXEDOCCMAP3_H_J5TFK2HQ */
