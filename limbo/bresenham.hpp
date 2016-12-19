/**
 * Trace a straight line from start_pos to end_pos.
 * At each step fun(i, j, k) is called.
 *
 * NOTE start_pos and end_pos should be inside the grid.
 *
 * Reference: graphics gems article
 * TODO consider DDA-type raytracing.
 */
template<class TraceFunctor>
void bresenham_trace3(const Vec3Ix& start_pos,
                      const Vec3Ix& end_pos,
                      const TraceFunctor& fun) {
#define UPDATE_CELL_HIT
#define UPDATE_CELL_PASS \
  if (!fun(x, y, z)) { break; }
#include "bresenham3_macro.def"
#undef UPDATE_CELL_HIT
#undef UPDATE_CELL_PASS
}


template<class GridT, class ArrayT>
void bresenham_trace3_increment(const Vec3Ix& start_pos,
                                const Vec3Ix& end_pos,
                                const GridT& grid3,
                                ArrayT& array3) {
#define UPDATE_CELL_HIT
#define UPDATE_CELL_PASS \
  mem_ix_t mem_ix = grid3.grid_to_mem(x, y, z); \
  array3[mem_ix] += 1;
#include "bresenham3_macro.def"
#undef UPDATE_CELL_HIT
#undef UPDATE_CELL_PASS
}

template<class TraceFunctor>
void bresenham_trace2(const Vec2Ix& start_pos,
                      const Vec2Ix& end_pos,
                      const TraceFunctor& fun) {

  int x = start_pos[0], y = start_pos[1];

  int dx = end_pos[0] - start_pos[0],
      dy = end_pos[1] - start_pos[1];

  int sx, sy;

  //X
  if (dx>0) {
    sx = 1;
  } else if (dx<0) {
    sx = -1;
    dx = -dx;
  } else {
    sx = 0;
  }

  //Y
  if (dy>0) {
    sy = 1;
  } else if (dy<0) {
    sy = -1;
    dy = -dy;
  } else {
    sy = 0;
  }

  int ax = 2*dx,
      ay = 2*dy;

  if (dy <= dx){
    for (int decy=ay-dx;;
         x+=sx, decy+=ay) {
      bool end_cell = false;
      //Bresenham step
      if(!fun(x, y, end_cell)) {
        break;
      }
      if (x==end_pos[0]) {
        break;
      }
      if (decy>=0) {
        decy-=ax;
        y+=sy;
      }
    }
  } else if (dx <= dy){
    for (int decx=ax-dy;;
         y+=sy,decx+=ax) {
      bool end_cell = false;
      //Bresenham step
      if (!fun(x,y,end_cell)){
        break;
      }

      if (y==end_pos[1]) break;
      if (decx>=0) {
        decx-=ay;
        x+=sx;
      }
    }
  }
}


