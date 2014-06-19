
#include <iostream>

#include <boost/foreach.hpp>
#include <boost/assign/std/vector.hpp>

#include <Eigen/Core>

#include <scrollgrid/scrollgrid3.hpp>

#include <gtest/gtest.h>

using namespace boost::assign;
using namespace Eigen;
using namespace ca;

TEST(scrollgrid3, hash_roundtrip1) {

  Vector3f center(-10, 20, -30);
  Vec3Ix dim(200, 100, 400);
  float res = 0.15;
  ca::ScrollGrid3f grid3( center, dim, res);

  std::vector<grid_ix_t> ix;
  ix += -10, -3, 0, 4, 11;

  BOOST_FOREACH(grid_ix_t i, ix) {
    BOOST_FOREACH(grid_ix_t j, ix) {
      BOOST_FOREACH(grid_ix_t k, ix) {
        Vec3Ix gix(i, j, k);
        uint64_t hix = grid3.grid_to_hash(gix);
        Vec3Ix gix2 = grid3.hash_to_grid(hix);
        EXPECT_EQ(gix, gix2);
      }
    }
  }

}

TEST(scrollgrid3, hash_roundtrip2) {

  Vector3f center(-10, 20, -30);
  Vec3Ix dim(21, 31, 41);
  float res = 0.15;
  ca::ScrollGrid3f grid3( center, dim, res);

  std::vector<grid_ix_t> ix;
  ix += -10, -3, 0, 4, 11;

  BOOST_FOREACH(grid_ix_t i, ix) {
    BOOST_FOREACH(grid_ix_t j, ix) {
      BOOST_FOREACH(grid_ix_t k, ix) {
        Vec3Ix gix(i, j, k);
        uint64_t hix = grid3.grid_to_hash(gix);
        Vec3Ix gix2 = grid3.hash_to_grid(hix);
        EXPECT_EQ(gix, gix2);
      }
    }
  }

}

TEST(scrollgrid3, mem_roundtrip1) {

  Vector3f center(-10, 20, -30);
  Vec3Ix dim(200, 100, 400);
  float res = 0.15;
  ca::ScrollGrid3f grid3( center, dim, res);

  std::vector<grid_ix_t> ix;
  ix += -10, -3, 0, 4, 11;

  BOOST_FOREACH(grid_ix_t i, ix) {
    BOOST_FOREACH(grid_ix_t j, ix) {
      BOOST_FOREACH(grid_ix_t k, ix) {
        Vec3Ix gix(i, j, k);
        mem_ix_t mix = grid3.grid_to_mem_slow(gix);
        Vec3Ix gix2 = grid3.mem_to_grid(mix);
        std::cerr << "gix = " << gix.transpose() << "|" << " gix2 = " << gix2.transpose() << "\n";
        EXPECT_TRUE( gix.cwiseEqual(gix2).all() );
      }
    }
  }


}

TEST(scrollgrid3, mem_roundtrip2) {

  Vector3f center(-10, 20, -30);
  Vec3Ix dim(21, 31, 41);
  float res = 0.15;
  ca::ScrollGrid3f grid3( center, dim, res);

  std::vector<grid_ix_t> ix;
  ix += -10, -3, 0, 4, 11;

  BOOST_FOREACH(grid_ix_t i, ix) {
    BOOST_FOREACH(grid_ix_t j, ix) {
      BOOST_FOREACH(grid_ix_t k, ix) {
        Vec3Ix gix(i, j, k);
        mem_ix_t mix = grid3.grid_to_mem_slow(gix);
        Vec3Ix gix2 = grid3.mem_to_grid(mix);
        EXPECT_TRUE( gix.cwiseEqual(gix2).all() );
      }
    }
  }

}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
