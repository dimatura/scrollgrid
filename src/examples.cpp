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

  csg::OccMap<ca::HitPassI4, csg::FixedGrid, ca::DenseArray, csg::HitPassUpdater, 2> hpmap2;
  csg::OccMap<float, csg::FixedGrid, ca::DenseArray, csg::BinaryFloatUpdater, 2> occmap2;

  return 0;
}
