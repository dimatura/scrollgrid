#include "scrollgrid/fixedoccmap3.h"
#include "scrollgrid/ray.hpp"
#include <scrollgrid/raycasting.hpp>

namespace ca { namespace scrollgrid {

void FixedOccMap3::Update(RowMatrixX3f& p, RowMatrixX3f& vp) {
  assert(p.rows() == vp.rows());

  const auto& box(_grid3.box());
  for (int i=0; i < p.rows(); ++i) {
    Eigen::Vector3f xyzf(p.row(i));
    Eigen::Vector3f originf(vp.row(i));

    Ray3<float> ray(originf, (xyzf-originf));
    if (!ca::aabb_ray_intersect(box, ray)) {
      // ray doesn't intersect box
      continue;
    }

    // TODO remove 1e-6 fudge factor
    Eigen::Vector3f start_ray = ray.point_at(ray.tmin+1e-6);
    Eigen::Vector3f end_ray = ray.point_at(ray.tmax-1e-6);
    if ((xyzf-originf).norm() < (start_ray-originf).norm()) {
      // ray ends before actually hitting box
      continue;
    }

    // adjust segment endpoints
    if (_grid3.is_inside_box(originf)) {
      start_ray = originf;
    }
    bool hit = false;
    if (grid3.is_inside_box(xyzf)) {
      hit = true;
      end_ray = xyzf;
    }

    ca::Vec3Ix start_grid_ix(_grid3.world_to_grid(start_ray));
    ca::Vec3Ix end_grid_ix(_grid3.world_to_grid(end_ray));

    // We check again because of border cases. see fudge factor above
    if (!grid3.is_inside_grid(end_grid_ix) ||
        !grid3.is_inside_grid(start_grid_ix)) {
      continue;
    }

    BresenhamOccupancyTrace(start_grid_ix, end_grid_ix, grid3, hit, occstats);
  }
}

} }
