#ifndef GRID_HASH_HPP_EXKMGXPF
#define GRID_HASH_HPP_EXKMGXPF

namespace ca { namespace scrollgrid

// TODO this is cut out from fixedgrid

// fixedgrid3
  /**
   * pack into grid_ix_t as [int16, int16, int16, 0]
   * note grid_ix_it is a signed 64-bit type
   * this gives range of [-32768, 32768] for each coordinate
   * so if voxel is 5cm, [-1638.4 m, 1638.4 m] relative to initial center.
   * this is useful as a unique hash.
   * this is better than mem_ix because mem_ix is ambiguous for absolute ijk.
   * if we use linear mem_ix as a hash,
   * then because of scrolling there may be collisions, and mem_ix become
   * invalidated once the corresponding voxel scrolls out.
   *
   */
  hash_ix_t grid_to_hash(const VecIx& grid_ix) const {
    //VecIx grid_positive(grid_ix - min_world_corner_ix_);
    Vec3Ix grid_positive(Vec3Ix::Zero());
    grid_positive.head<Dim>() = (grid_ix - min_world_corner_ix_);
    hash_ix_t hi = static_cast<hash_ix_t>(grid_positive[0]);
    hash_ix_t hj = static_cast<hash_ix_t>(grid_positive[1]);
    hash_ix_t hk = static_cast<hash_ix_t>(grid_positive[2]);
    hash_ix_t h = (hi << 48) | (hj << 32) | (hk << 16);
    return h;
  }

  Vec3Ix hash_to_grid(hash_ix_t hix) const {
    hash_ix_t hi = (hix & 0xffff000000000000) >> 48;
    hash_ix_t hj = (hix & 0x0000ffff00000000) >> 32;
    hash_ix_t hk = (hix & 0x00000000ffff0000) >> 16;
    Vec3Ix grid_ix(hi, hj, hk);
    grid_ix += min_world_corner_ix_;
    return grid_ix;
  }

  Mat3Ix multiple_hash_to_grid(const HashIxVector& hindices) const {
    Mat3Ix out(3, hindices.rows());
    for (size_t i=0; i < hindices.rows(); ++i) {
      Vec3Ix gix = this->hash_to_grid(hindices[i]);
      out.col(i) = gix;
    }
    return out;
  }

  /**
   * Like the above but does not offset by origin.
   * In the fixed case we assume the min ijk for origin is 0,0,0.
   */
  hash_ix_t local_grid_to_hash(const Vec3Ix& grid_ix) const {
    // TODO assumes grid_ix is inside box.
    hash_ix_t hi = static_cast<hash_ix_t>(grid_ix[0]);
    hash_ix_t hj = static_cast<hash_ix_t>(grid_ix[1]);
    hash_ix_t hk = static_cast<hash_ix_t>(grid_ix[2]);
    hash_ix_t h = (hi << 48) | (hj << 32) | (hk << 16);
    return h;
  }

  Vec3Ix hash_to_local_grid(hash_ix_t hix) const {
    hash_ix_t hi = (hix & 0xffff000000000000) >> 48;
    hash_ix_t hj = (hix & 0x0000ffff00000000) >> 32;
    hash_ix_t hk = (hix & 0x00000000ffff0000) >> 16;
    Vec3Ix grid_ix(hi, hj, hk);
    return grid_ix;
  }

  Mat3Ix multiple_local_hash_to_grid(const HashIxVector& hindices) const {
    Mat3Ix out(3, hindices.rows());
    for (size_t i=0; i < hindices.rows(); ++i) {
      Vec3Ix gix = this->local_hash_to_grid(hindices[i]);
      out.col(i) = gix;
    }
    return out;
  }


} }

#endif /* end of include guard: GRID_HASH_HPP_EXKMGXPF */

