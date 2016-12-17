#ifndef FIXEDOCCMAP2_H_FB1K65AT
#define FIXEDOCCMAP2_H_FB1K65AT

#include <memory>
#include <limits>

//#include <ros/ros.h>

//#include <pcl_util/point_types.hpp>

#include <scrollgrid/dense_array.hpp>
#include <scrollgrid/fixedgrid.hpp>
#include <scrollgrid/raycasting.hpp>
#include <scrollgrid/occmap_constants.hpp>

/**
 * FixedGrid + DenseArray = static occmap.
 */
namespace ca { namespace scrollgrid {

template <typename CellT>
class CellUpdater {
};

template <>
class CellUpdater<ca::HitPass> {
public:

  static
  void init(ca::HitPass& cell) {
    cell.hits = 0;
    cell.passes = 0;
  }

  static
  void update_pos(ca::HitPass& cell) {
    cell.hits += 1;
  }

  static
  void update_neg(ca::HitPass& cell) {
    cell.passes += 1;
  }
};


template <>
class CellUpdater<ca::BinaryOccupancy<float>> {
public:
  static constexpr float UNKNOWN = 0.5f;
  static constexpr float OCCUPIED = 0.9f;
  static constexpr float FREE = 0.0f;
  static constexpr float UPDATE_POS = 0.08f;
  static constexpr float UPDATE_NEG = 0.01f;

public:
  typedef std::shared_ptr<ca::HitPass> Ptr;

  static
  void init(ca::BinaryOccupancy<float>& cell) {
    cell.val = UNKNOWN;
  }

  static
  void update_pos(ca::BinaryOccupancy<float>& cell) {
    float val = cell.val;
    float new_val = std::min(val + UPDATE_POS, OCCUPIED);
    cell.val = new_val;
  }

  static
  void update_neg(ca::BinaryOccupancy<float>& cell) {
    float val = cell.val;
    float new_val = std::max(val - UPDATE_NEG, FREE);
    cell.val = new_val;
  }
};

template <typename CellT, int Dim>
class OccMap {
};

template <typename CellT>
class OccMap<CellT, 2> {
public:
  void compute_start_end_grid_ix(const Eigen::Vector2f& origin,
                                 const Eigen::Vector2f& x,
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

  void update(const Eigen::Vector2f& originf,
              const Eigen::Vector2f& xyf) {
    bool hit = false, intersects = false;
    ca::Vec2Ix start_grid_ix(ca::Vec2Ix::Zero());
    ca::Vec2Ix end_grid_ix(ca::Vec2Ix::Zero());
    this->compute_start_end_grid_ix(originf,
                                    xyf,
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
      if (b2itr.done() && hit) {
        CellT& cell(storage_[mem_ix]);
        CellT::update_pos(cell);
      } else {
        CellT& cell(storage_[mem_ix]);
        CellT::update_neg(cell);
      }
    }
  }

  void multi_update(const Eigen::Matrix2Xf& vp,
                    const Eigen::Matrix2Xf& p) {
    ROS_ASSERT(p.rows() == vp.rows());
    for (int i=0; i < p.cols(); ++i) {
      Eigen::Vector2f xyf(p.col(i));
      Eigen::Vector2f originf(vp.col(i));
      this->update(originf, xyf);
    }
  }

  FixedGrid2f grid_;
  DenseArray<CellT, 2> storage_;

};


#if 0
template <>
class HitPassCore<float> {
public:
  typedef std::shared_ptr<HitPassCore> Ptr;

public:
  static constexpr float UNKNOWN = 0.5f;
  static constexpr float OCCUPIED = 0.9f;
  static constexpr float FREE = 0.0f;
  static constexpr float UPDATE_POS = 0.08f;
  static constexpr float UPDATE_NEG = 0.01f;

public:
  HitPassCore() { }
  virtual ~HitPassCore() { }
  HitPassCore(const HitPassCore& other) = delete;
  HitPassCore& operator=(const HitPassCore& other) = delete;

public:

  void init(const Eigen::Matrix<grid_ix_t, Dim, 1>& dims) {
    occ_.reset(dims);
    occ_.fill(UNKNOWN);
  }

  void reset() {
    occ_.fill(UNKNOWN);
  }

  void update_pos(mem_ix_t mem_ix) {
    Scalar val = occ_[mem_ix];
    Scalar new_val = std::min(val + UPDATE_POS, OCCUPIED);
    occ_[mem_ix] = new_val;
  }

  void update_neg(mem_ix_t mem_ix) {
    Scalar val = occ_[mem_ix];
    Scalar new_val = std::max(val - UPDATE_NEG, FREE);
    occ_[mem_ix] = new_val;
  }

  ca::DenseArray<Scalar, Dim> occ_;
};
#endif

#if 0
template<typename Scalar, int Dim>
class OccMapCore {
public:
  typedef std::shared_ptr<OccMapCore> Ptr;

public:
  static constexpr float UNKNOWN = 0.5f;
  static constexpr float OCCUPIED = 0.9f;
  static constexpr float FREE = 0.0f;
  static constexpr float UPDATE_POS = 0.08f;
  static constexpr float UPDATE_NEG = 0.01f;

public:

  void init(const Eigen::Matrix<grid_ix_t, Dim, 1>& dims) {
    occ_.reset(dims);
    occ_.fill(UNKNOWN);
  }

  void reset() {
    occ_.fill(UNKNOWN);
  }

  void update_pos(mem_ix_t mem_ix) {
    Scalar val = occ_[mem_ix];
    Scalar new_val = std::min(val + UPDATE_POS, OCCUPIED);
    occ_[mem_ix] = new_val;
  }

  void update_neg(mem_ix_t mem_ix) {
    Scalar val = occ_[mem_ix];
    Scalar new_val = std::max(val - UPDATE_NEG, FREE);
    occ_[mem_ix] = new_val;
  }

private:
  ca::DenseArray<Scalar, Dim> occ_;
};
#endif


#if 0
template<typename Scalar, int Dim>
class HitPassMapCore {
public:
  typedef std::shared_ptr<HitPassMapCore> Ptr;

public:
  void init(const Eigen::Matrix<grid_ix_t, Dim, 1>& dims) {
    hits_.reset(dims);
    passes_.reset(dims);
    hits_.fill(0);
    passes_.fill(0);
  }

  void reset() {
    hits_.fill(0);
    passes_.fill(0);
  }

  void update_pos(mem_ix_t mem_ix) {
    hits_[mem_ix] += 1;
  }

  void update_neg(mem_ix_t mem_ix) {
    passes_[mem_ix] += 1;
  }

  ca::DenseArray<Scalar, Dim> get_hits_array() const {
    return hits_;
  }

  ca::DenseArray<Scalar, Dim> get_passes_array() const {
    return passes_;
  }

private:
  ca::DenseArray<Scalar, Dim> hits_;
  ca::DenseArray<Scalar, Dim> passes_;

};
#endif

#if 0
template <class CoreT>
class OccMap2 {
public:
  OccMap2() {
  }
  virtual ~OccMap2() { }
  OccMap2(const OccMap2& other) = delete;
  OccMap2& operator=(const OccMap2& other) = delete;

  void compute_start_end_grid_ix(const Eigen::Vector2f& origin,
                                 const Eigen::Vector2f& x,
                                 ca::Vec2Ix& start_grid_ix,
                                 ca::Vec2Ix& end_grid_ix,
                                 bool& hit,
                                 bool& intersects) {

    // inter1 and inter2 are clipped versions of origin, x
    Eigen::Vector2f inter1(Eigen::Vector2f::Zero());
    Eigen::Vector2f inter2(Eigen::Vector2f::Zero());

    const ca::scrollgrid::FixedGrid2f& grid(occmap_.grid());
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

  void update(const Eigen::Vector2f& originf,
              const Eigen::Vector2f& xyf) {
    bool hit = false, intersects = false;
    ca::Vec2Ix start_grid_ix(ca::Vec2Ix::Zero());
    ca::Vec2Ix end_grid_ix(ca::Vec2Ix::Zero());
    this->compute_start_end_grid_ix(originf,
                                    xyf,
                                    start_grid_ix,
                                    end_grid_ix,
                                    hit,
                                    intersects);
    if (!intersects) {
      return;
    }

    const ca::scrollgrid::FixedGrid2f& grid(occmap_.grid());
    // TODO just hit at end of the ray
    Bresenham2Iterator b2itr(start_grid_ix, end_grid_ix);
    while (!b2itr.done()) {
      b2itr.step();
      ca::Vec2Ix ij(b2itr.pos());
      if (!grid.is_inside_grid(ij)) {
        break;
      }

      mem_ix_t mem_ix = grid.grid_to_mem(ij);
      if (b2itr.done() && hit) {
        occmap_.update_pos(mem_ix);
      } else {
        occmap_.update_neg(mem_ix);
      }
    }
  }

  void multi_update(const Eigen::Matrix2Xf& vp,
                    const Eigen::Matrix2Xf& p) {
    ROS_ASSERT(p.rows() == vp.rows());
    for (int i=0; i < p.cols(); ++i) {
      Eigen::Vector2f xyf(p.col(i));
      Eigen::Vector2f originf(vp.col(i));
      occmap_.update(xyf, originf);
    }
  }

  OccMapT occmap_;
};
#endif


} }

#endif /* end of include guard: FIXEDOCCMAP2_H_FB1K65AT */
