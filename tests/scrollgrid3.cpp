#include <iostream>
#include <Eigen/Core>
#include <scrollgrid/scrollgrid3.hpp>

#include <gtest/gtest.h>

using namespace Eigen;
using namespace ca;

TEST(scrollgrid3, hash_roundtrip1) {

  Vector3f center(-10, 20, -30);
  Vec3Ix dim(200, 100, 400);
  float res = 0.15;
  ca::ScrollGrid3f grid3( center, dim, res);

  grid_ix_t ix[5] = { -10, -3, 0, 4, 11 };

  for (grid_ix_t ii=0; ii < 5; ++ii) {
    for (grid_ix_t jj=0; jj < 5; ++jj) {
      for (grid_ix_t kk=0; kk < 5; ++kk) {
        Vec3Ix gix(ix[ii], ix[jj], ix[kk]);
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

  grid_ix_t ix[5] = { -10, -3, 0, 4, 11 };
  for (grid_ix_t ii=0; ii < 5; ++ii) {
    for (grid_ix_t jj=0; jj < 5; ++jj) {
      for (grid_ix_t kk=0; kk < 5; ++kk) {
        Vec3Ix gix(ix[ii], ix[jj], ix[kk]);
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

  grid_ix_t ix[5] = { -10, -3, 0, 4, 11 };
  for (grid_ix_t ii=0; ii < 5; ++ii) {
    for (grid_ix_t jj=0; jj < 5; ++jj) {
      for (grid_ix_t kk=0; kk < 5; ++kk) {
        Vec3Ix gix(ix[ii], ix[jj], ix[kk]);
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

  grid_ix_t ix[5] = { -10, -3, 0, 4, 11 };
  for (grid_ix_t ii=0; ii < 5; ++ii) {
    for (grid_ix_t jj=0; jj < 5; ++jj) {
      for (grid_ix_t kk=0; kk < 5; ++kk) {
        Vec3Ix gix(ix[ii], ix[jj], ix[kk]);
        mem_ix_t mix = grid3.grid_to_mem(gix);
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
