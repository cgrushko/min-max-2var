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

#include "solver-main.h"
#include "my_utilities.h"
#include "solver2d.h"
#include "solver3d.h"

#include <limits>
#include <algorithm>

extern "C" void exactinit();

namespace hough_min_max {
  
  namespace {
    bool isValid(double u) {
      return u >= 0 && u <= 1;
    }
    
    void transformProblemToXY1(P3* points, int n) {
      for (int i=0; i<n; i++) {
        double b = points[i].y;
        points[i].x -= b;
        points[i].z -= b;
      }
    }
    
  } // anonymous namespace
  
  void init() {
    exactinit();
  }
  
  void solve(P3* points, int n, double* u, double* v, double* value) {
    Result result_unc;
    Result results[6];

    // Find unconstrained minimum
    solver_3d::solve(points, n, &result_unc);
    
    if (isValid(result_unc.u) && isValid(result_unc.v) && result_unc.u + result_unc.v <= 1) {
      *u = result_unc.u;
      *v = result_unc.v;
      *value = result_unc.value;
      return;
    }
    
    // Find minimum on x = 0 and y = 0 edges
    solver_2d::solve_on_edges(points, n, &results[3], &results[4]);
    if (!isValid(results[3].u))
      results[3].value = std::numeric_limits<double>::infinity();
    if (!isValid(results[4].v))
      results[4].value = std::numeric_limits<double>::infinity();
    
    // Find minimum on corners
    solver_2d::find_max_on_corners(points, n, &results[0], &results[1], &results[2]);

    // Find minimum on the line x+y = 1
    transformProblemToXY1(points, n);
    solver_2d::solve_on_x(points, n, &results[5]);
    if (!isValid(results[5].u))
      results[5].value = std::numeric_limits<double>::infinity();
    else {
      results[5].v = 1-results[5].u;
    }
    
    // Choose minimum
    Result* min = std::min_element(results, results + 6);
    *u = min->u;
    *v = min->v;
    *value = min->value;
  }
  
} // namespace hough_min_max
