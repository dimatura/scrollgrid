#ifndef FIXEDOCCMAP3_H_J5TFK2HQ
#define FIXEDOCCMAP3_H_J5TFK2HQ

#include <memory>

#include <ros/ros.h>

#include <pcl_util/point_types.hpp>
#include <scrollgrid/dense_array3.hpp>
#include <scrollgrid/fixedgrid3.hpp>

/**
 * FixedGrid3 + DenseArray3 = static occmap.
 * TODO: uint8 for now. Do we need float?
 */
namespace ca { namespace scrollgrid {


class FixedOccMap3 {

public:
  typedef std::shared_ptr<FixedOccMap3> Ptr;

  static constexpr uint8_t UNKNOWN = 128;

private:
  ca::FixedGrid3f grid_;
  ca::DenseArray3<uint8_t> occstats_;


public:
  FixedOccMap3(const ca::FixedGrid3f& grid) {
    grid_ = grid;
    this->Init();
  }

  FixedOccMap3(const FixedOccMap3& other) = delete;
  FixedOccMap3& operator=(const FixedOccMap3& other) = delete;

  void Init() {
    occstats_.reset(grid_.dimension());
    occstats_.fill(static_cast<uint8_t>(UNKNOWN));
  }

  void SetUnknown() {
    occstats_.fill(static_cast<uint8_t>(UNKNOWN));
  }

  void Update(RowMatrixX3f& p, RowMatrixX3f& vp) {
#if 0
    pcl::PointCloud<pcl::PointWithViewpoint> xyzwvp;
    ca::NeaPc2ToPclXyzvp(cloud, xyzwvp);
    //ROS_INFO_STREAM("xyzwvp.size() = " << xyzwvp.size());
    UpdateOccupancy(xyzwvp, grid_, occstats_);
#endif
  }

  void CopyToBuffer(uint8_t * buffer) const {
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
