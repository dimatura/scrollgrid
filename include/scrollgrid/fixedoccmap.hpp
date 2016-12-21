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
    cell.hits = std::min(cell.hits+1, std::numeric_limits<int32_t>::max()-1);
  }

  static
  void update_neg(ca::HitPassI4& cell) {
    cell.passes = std::min(cell.passes+1, std::numeric_limits<int32_t>::max()-1);
  }

  static
  void update_decay(ca::HitPassI4& cell) {
    cell.hits = std::max(cell.hits-1, 0);
    cell.passes = std::max(cell.passes-1, 0);
  }

};


class BinaryFloatUpdater {
public:
  // constants from octomap paper
  // TODO: maybe constexpr log(x/(1.-x))
  static constexpr float UNKNOWN = 0.0; // 0.5
  static constexpr float OCCUPIED = 3.5; // 0.97
  static constexpr float FREE = -2.; // 0.12
  static constexpr float UPDATE_POS = 0.85; // 0.7
  static constexpr float UPDATE_NEG = -0.4; // 0.4
  static constexpr float UPDATE_DECAY = -0.05; // 0.4

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
    cell = std::max(cell + UPDATE_NEG, FREE);
  }

  static
  void update_decay(float& cell) {
    cell = std::max(cell + UPDATE_DECAY, UNKNOWN);
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
  typedef std::shared_ptr<OccMap> Ptr;
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
    if (!intersects) { return; }

    //std::cerr << "originf = " << originf.transpose() << std::endl;
    //std::cerr << "xyf = " << xyf.transpose() << std::endl;
    //std::cerr << "start_grid_ix = " << start_grid_ix.transpose() << std::endl;
    //std::cerr << "end_grid_ix = " << end_grid_ix.transpose() << std::endl;

    if (hit) {
      mem_ix_t mem_ix = grid_.grid_to_mem(end_grid_ix);
      CellT& cell(storage_[mem_ix]);
      UpdaterT::update_pos(cell);
    }

    BresenhamIterator<Dim> bitr(start_grid_ix, end_grid_ix);
    do {
      bitr.step();
      VecIx grid_ix(bitr.pos());

      if (!grid_.is_inside_grid(grid_ix)) {
        // sometimes there's edge effects with clipping
        continue;
      }

      // TODO consider other policies that map
      // grid_ix -> key
      // where key is not just memory. Or encapsulate
      // like map.update(grid_ix)
      if (!bitr.done()) {
        mem_ix_t mem_ix = grid_.grid_to_mem(grid_ix);
        CellT& cell(storage_[mem_ix]);
        UpdaterT::update_neg(cell);
      }
    } while (!bitr.done());
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

  void multi_update_single_vp(const Vecf& vp,
                              const Matf& p) {
    ROS_ASSERT(p.rows() == vp.rows());
    for (int i=0; i < p.cols(); ++i) {
      Vecf xyf(p.col(i));
      this->update(vp, xyf);
    }
  }


  void decay() {
    for (int i=0; i < grid_.num_cells(); ++i) {
      CellT& cell(storage_[i]);
      UpdaterT::update_decay(cell);
    }
  }


  StorageT<CellT, Dim>& get_storage() {
    return storage_;
  }

  GridT<float, Dim>& get_grid() {
    return grid_;
  }

 private:
  void compute_start_end_grid_ix(const Vecf& origin,
                                 const Vecf& x,
                                 VecIx& start_grid_ix,
                                 VecIx& end_grid_ix,
                                 bool& hit,
                                 bool& intersects) {

    // inter1 and inter2 are clipped versions of origin, x

    auto box(grid_.box());
#if 1
    Vecf inter1(Vecf::Zero());
    Vecf inter2(Vecf::Zero());
    intersects = box.clip_line(origin, x, inter1, inter2);
#else
    Vecf inter1 = origin;
    Vecf inter2 = x;
    intersects = true;
#endif

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

typedef OccMap<float, FixedGrid, ca::DenseArray, BinaryFloatUpdater, 2> FixedBinMap2f;
typedef OccMap<ca::HitPassI4, FixedGrid, ca::DenseArray, HitPassUpdater, 2> FixedHitPassMap2i;
typedef OccMap<float, FixedGrid, ca::DenseArray, BinaryFloatUpdater, 3> FixedBinMap3f;
typedef OccMap<ca::HitPassI4, FixedGrid, ca::DenseArray, HitPassUpdater, 3> FixedHitPassMap3i;


} }

#endif /* end of include guard: FIXEDOCCMAP_HPP_I9FPFYWV */
