/*
 * This file is part of minmax2var - A two-variable linar min-max solver.
 * 
 * Copyright (C) 2012  Carmi Grushko, carmi.grushko@gmail.com
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 */
//  Created by Carmi Grushko on 6/29/12.

#ifndef minmax2var_solver2d_h
#define minmax2var_solver2d_h

#include "P3.h"
#include <utility>
#include "my_utilities.h"

namespace hough_min_max {
  namespace solver_2d {
    
    // Find minimum on corners: (0,0), (1,0), (0,1)
    void find_max_on_corners(P3 const *a, int n, Result *result_corner_x, Result *result_corner_y, Result *result_corner_0);
    
    // Find minimum on the lines x = 0 and y = 0;
    void solve_on_edges(const P3* ptArray, int n, Result *result_x, Result* result_y);
    
    // Find minimum just on the line x = 0;
    void solve_on_x(const P3* ptArray, int n, Result *result_x);
    
  } // namespace solver_2d
} // namespace hough_min_max

#endif
