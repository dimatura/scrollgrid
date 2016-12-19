/**
 * Copyright (c) 2015 Carnegie Mellon University, Daniel Maturana <dimatura@cmu.edu>
 *
 * For License information please see the LICENSE file in the root directory.
 *
 */

#include "scrollgrid/dense_array.hpp"
#include "scrollgrid/raycasting.hpp"
#include "scrollgrid/fixedgrid.hpp"
#include "scrollgrid/fixedoccmap.hpp"

int main(int argc, char *argv[]) {

  namespace csg = ca::scrollgrid;

  //csg::OccMap<ca::HitPass, csg::FixedGrid2f, csg::HitPassUpdater, 2> hpmap;
  //csg::OccMap<float, csg::FixedGrid2f, csg::BinaryFloatUpdater, 2> occmap;
  //ca::Vec2Ix dims(28, 28);
  //hpmap.init(dims);

  csg::OccMap<ca::HitPass, csg::FixedGrid, ca::DenseArray, csg::HitPassUpdater, 2> hpmap2;
  csg::OccMap<float, csg::FixedGrid, ca::DenseArray, csg::BinaryFloatUpdater, 2> occmap2;

  return 0;
}

#if 0
/**
 * I think I defined one already somewhere but forgot where
 */
class ClearFunctor : public ca::ClearCellsFun {
public:
  ClearFunctor(const ca::ScrollGrid3f* grid,
               ca::DenseArray<uint8_t, 3>* occ_array) {
    grid_ = grid;
    occ_array_ = occ_array;
  }

  virtual void operator()(const ca::Vec3Ix& start,
                          const ca::Vec3Ix& finish) const {
    if ((finish - start).prod() > 0) {
      ROS_ASSERT(grid_->is_inside_grid(start));
      ROS_ASSERT(grid_->is_inside_grid(finish-ca::Vec3Ix(1,1,1)));
      for (ca::grid_ix_t i = start[0]; i < finish[0]; ++i) {
        for (ca::grid_ix_t j = start[1]; j < finish[1]; ++j) {
          for (ca::grid_ix_t k = start[2]; k < finish[2]; ++k) {
            ca::mem_ix_t mem_ix = grid_->grid_to_mem(i, j, k);
            (*occ_array_)[mem_ix] = 0; // application-dependent
          }
        }
      }
    }
  }

  virtual ~ClearFunctor() { }
  const ca::ScrollGrid3f * grid_;
  ca::DenseArray<uint8_t, 3> * occ_array_;
};

int main(int argc, char *argv[]) {

  // scroll grid maps with origin in (0, 0, 0),
  // with dimensions (200, 200, 200), with
  // voxel size 0.5 (units are arbitrary but probably meters)
  // scrollgrid does not hold any data. it just maps coordinates.
  ca::ScrollGrid3f grid3(Eigen::Vector3f(0, 0, 0),
                         ca::Vec3Ix(200, 200, 200),
                         0.5);

  // a dense occupancy grid, type uint8. uses same dimensions as grid3.
  // doing otherwise would end badly.
  // a sparse array or any other backend could be used here.
  ca::DenseArray<uint8_t, 3> occ_array3(grid3.dimension());

  // fill with zeros
  occ_array3.fill(0);

  // manipulate the array directly.
  // if you go out of bounds, it will crash.
  occ_array3.get({100, 100, 100}) = 255;
  std::cout << "occarray3(100, 100, 100) " << int(occ_array3.get({100, 100, 100})) << "\n";

  // linear memory address also works, but does not know anything
  // about location in space (hence 'local')
  ca::mem_ix_t mix = occ_array3.local_grid_to_mem({100, 100, 100});
  std::cout << "occarray3[" << mix << "] " << int(occ_array3[mix]) << "\n";

  // use bresenham to draw a line. again, directly in voxel space
  ca::bresenham_trace3_increment(ca::Vec3Ix(20, 40, 10),
                                 ca::Vec3Ix(180, 170, 190),
                                 grid3,
                                 occ_array3);

  // now get an xyz point
  Eigen::Vector3f pt(10., 20., 30.);

  ca::Bresenham3Iterator b3itr(ca::Vec3Ix(20, 40, 10),
                               ca::Vec3Ix(180, 170, 190));
  while (!b3itr.done()) {
    b3itr.step();
    bool hit = b3itr.done();
    ca::mem_ix_t mix = grid3.grid_to_mem(b3itr.pos());
    if (hit) {
      occ_array3[mix] = occ_array3[mix] + 1;
    } else {
      occ_array3[mix] = occ_array3[mix] - 1;
    }
  }

  for (ca::Bresenham3Iterator b3itr(ca::Vec3Ix(20, 40, 10),
                                    ca::Vec3Ix(180, 170, 190));
      !b3itr.done();
      b3itr.step()) {
    bool hit = b3itr.done();
    ca::mem_ix_t mix = grid3.grid_to_mem(b3itr.pos());
  }

  // verify it is inside the volume
  if (grid3.is_inside_box(pt)) {
    std::cout << "inside\n";
  } else {
    std::cout << "outside\n";
  }

  // try with one outside
  if (grid3.is_inside_box(Eigen::Vector3f(1000., 0., 0.))) {
    std::cout << "inside\n";
  } else {
    std::cout << "outside\n";
  }

  // get discretized coordinates; this is *not* wrapped around,
  // so might be outside grid.
  ca::Vec3Ix ix = grid3.world_to_grid(pt);
  std::cout << ix.transpose() << "\n";

  // get memory coordinates. this is actually wrapped around
  // in space. world_to_grid -> grid_to_mem would also work.
  ca::mem_ix_t mix2 = grid3.world_to_mem(pt);
  std::cout << "mix2: " << mix2 << "\n";
  occ_array3[mix2] = 128;

  ClearFunctor clear_fn(&grid3, &occ_array3);

  // scroll by 10 cells, clear with clear_fn callback
  ca::Vec3Ix offset(10, 0, 0);
  grid3.scroll_and_clear(offset, clear_fn);

  // ijk never changes with scroll; mem_ix might or might not
  ca::Vec3Ix ix2 = grid3.world_to_grid(pt);
  std::cout << "ix2.T " << ix2.transpose() << "\n";

  // should retain same value
  std::cout << "occ_array3[grid3.world_to_mem(pt)]: " << int(occ_array3[grid3.world_to_mem(pt)]) << "\n";


  return 0;
}
#endif
