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
//  Created by Carmi Grushko on 7/1/12.

#ifndef minmax2var_solver3d_h
#define minmax2var_solver3d_h

#include "P3.h"
#include <utility>
#include "my_utilities.h"

namespace hough_min_max {
  
  namespace solver_3d {
    
    bool solve(const P3* ptArray, int n, Result *result);
    
  } // namespace solver_3d
} // namespace hough_min_max


#endif
