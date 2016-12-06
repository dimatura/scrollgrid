# scrollgrid

** Author **: Daniel Maturana (dimatura@cmu.edu)

** Maintainer **: Daniel Maturana (dimatura@cmu.edu)

`scrollgrid` is best scrolling grid

`scrollgrid` is a 2D and a 3D scrolling grid to store volumetric data.

See `src/examples.cpp` and `tests/scrollgrid3.cpp` for examples.


## Concepts

```
"world_xyz" xyz world frame
"grid_xyz" xyz grid frame (shifted relative to world_xyz by origin)
"grid_ijk" scaled and discretized grid coordinates:
    i = floor( (x-origin_x-0.5)/resolution ) 
"local_ijk" grid_ijk shifted by scrolling and limited/wrapped to local extent
    li = (i - scroll_offset_i) modulo (dim_i)
"mem_ix" index into flat storage from local_ijk.
    mem_ix = local_ijk.dot(strides)
"hash_ix": bit-packed version of grid_ijk
```


### License ###
[This software is BSD licensed.](http://opensource.org/licenses/BSD-3-Clause)
 
Copyright (c) 2015, Carnegie Mellon University
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
