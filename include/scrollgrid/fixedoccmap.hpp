#ifndef FIXEDOCCMAP_HPP_I9FPFYWV
#define FIXEDOCCMAP_HPP_I9FPFYWV

#include <memory>
#include <limits>

#include <ros/ros.h>

//#include <pcl_util/point_types.hpp>

#include <scrollgrid/dense_array.hpp>
#include <scrollgrid/fixedgrid.hpp>
#include <scrollgrid/raycasting.hpp>

/**
 * FixedGrid + DenseArray = static occmap.
 */
namespace ca { namespace scrollgrid {

class HitPassUpdater {
public:

  static
  void init(ca::HitPassI4& cell) {
    cell.hits = 0;
    cell.passes = 0;
  }

  static
  void update_pos(ca::HitPassI4& cell) {
    cell.hits += 1;
  }

  static
  void update_neg(ca::HitPassI4& cell) {
    cell.passes += 1;
  }
};


class BinaryFloatUpdater {
public:
  static constexpr float UNKNOWN = 0.5f;
  static constexpr float OCCUPIED = 0.9f;
  static constexpr float FREE = 0.0f;
  static constexpr float UPDATE_POS = 0.08f;
  static constexpr float UPDATE_NEG = 0.01f;

public:

  static
  void init(float& cell) {
    cell = UNKNOWN;
  }

  static
  void update_pos(float& cell) {
    cell = std::min(cell + UPDATE_POS, OCCUPIED);
  }

  static
  void update_neg(float& cell) {
    cell = std::max(cell - UPDATE_NEG, FREE);
  }
};

class BinaryUint8Updater {
public:
  static constexpr uint8_t UNKNOWN = 127;
  static constexpr uint8_t OCCUPIED = 250;
  static constexpr uint8_t FREE = 10;
  static constexpr uint8_t UPDATE_POS = 8;
  static constexpr uint8_t UPDATE_NEG = 2;

public:

  static
  void init(uint8_t& cell) {
    cell = UNKNOWN;
  }

  static
  void update_pos(uint8_t& cell) {
    int32_t new_value = static_cast<int32_t>(cell) + static_cast<int32_t>(UPDATE_POS);
    cell = static_cast<uint8_t>(std::min(static_cast<int32_t>(OCCUPIED), new_value));
  }

  static
  void update_neg(uint8_t& cell) {
    int32_t new_value = static_cast<int32_t>(cell) - static_cast<int32_t>(UPDATE_NEG);
    cell = static_cast<uint8_t>(std::max(static_cast<int32_t>(FREE), new_value));
  }
};


template <class CellT,
          template <class, int> class GridT,
          template <class, int> class StorageT,
          class UpdaterT,
          int Dim>
class OccMap {
 public:
  typedef CellT CellType;
  typedef Eigen::Matrix<float, Dim, 1> Vecf;
  typedef Eigen::Matrix<float, Dim, Eigen::Dynamic> Matf;
  typedef Eigen::Matrix<grid_ix_t, Dim, 1> VecIx;

 public:
  OccMap() {

  }

  OccMap(const Vecf& center,
         const VecIx& dimension,
         float resolution) :
      grid_(center, dimension, resolution)
  {

    storage_.reset(dimension);

    CellT cell;
    UpdaterT::init(cell);
    storage_.fill(cell);
  }

  void reset() {
    CellT cell;
    UpdaterT::init(cell);
    storage_.fill(cell);
  }

  void update(const Vecf& originf,
              const Vecf& xyf) {
    bool hit = false, intersects = false;
    VecIx start_grid_ix(VecIx::Zero());
    VecIx end_grid_ix(VecIx::Zero());
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
    BresenhamIterator<Dim> bitr(start_grid_ix, end_grid_ix);
    while (!bitr.done()) {
      bitr.step();
      VecIx grid_ix(bitr.pos());
      if (!grid_.is_inside_grid(grid_ix)) {
        break;
      }

      // TODO considier other policies that map
      // grid_ix -> key
      // where key is not  ust memory. Or encapsulate
      // like map.update(grid_ix)
      mem_ix_t mem_ix = grid_.grid_to_mem(grid_ix);
      if (bitr.done() && hit) {
        CellT& cell(storage_[mem_ix]);
        UpdaterT::update_pos(cell);
      } else {
        CellT& cell(storage_[mem_ix]);
        UpdaterT::update_neg(cell);
      }
    }
  }

  void multi_update(const Matf& vp,
                    const Matf& p) {
    ROS_ASSERT(p.rows() == vp.rows());
    for (int i=0; i < p.cols(); ++i) {
      Vecf xyf(p.col(i));
      Vecf originf(vp.col(i));
      this->update(originf, xyf);
    }
  }

  StorageT<CellT, Dim>& get_storage() {
    return storage_;
  }

 private:
  void compute_start_end_grid_ix(const Vecf& origin,
                                 const Vecf& x,
                                 VecIx& start_grid_ix,
                                 VecIx& end_grid_ix,
                                 bool& hit,
                                 bool& intersects) {

    // inter1 and inter2 are clipped versions of origin, x
    Vecf inter1(Vecf::Zero());
    Vecf inter2(Vecf::Zero());

    auto box(grid_.box());
    intersects = box.clip_line(origin, x, inter1, inter2);

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
  GridT<float, Dim> grid_;
  StorageT<CellT, Dim> storage_;

};

} }

#endif /* end of include guard: FIXEDOCCMAP_HPP_I9FPFYWV */
