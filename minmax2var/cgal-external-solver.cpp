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
//  Created by Carmi Grushko on 8/6/12.

#include "settings.h"

#ifdef EXTERNAL_SOLVER_CGAL
#include "cgal-external-solver.h"
#include <sys/time.h>

CGAL_Solver::CGAL_Solver() : n(-1) {
  c[0] = c[1] = 0;
  c[2] = 1;
}

void CGAL_Solver::init(int n) {
  if (n == this->n)
    return;

  for (int i=0; i<3; i++) {
    A_vectors[i].resize(n + 1); // +1 for the u+v<=1 constraint
    A[i] = A_vectors[i].data();
  }
  b.resize(n + 1);
  this->n = n;
  lp = new Program(3, n + 1, A, b.data(), CGAL::Const_oneset_iterator<CGAL::Comparison_result>(CGAL::SMALLER), c);
}

void CGAL_Solver::solve(const hough_min_max::P3 points[]) {
  for(int i=0; i<n; i++) {
    A[0][i] = points[i].x;
    A[1][i] = points[i].y;
    A[2][i] = -1;
    b[i] = points[i].z;
  }
  A[0][n] = 1;
  A[1][n] = 1;
  A[2][n] = 0;
  b[n] = 1;
  
  timeval t1, t2;
  gettimeofday(&t1, NULL);
  Solution s = CGAL::solve_nonnegative_linear_program(*lp, ET());
  gettimeofday(&t2, NULL);
  time += (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  time += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

  release_assert(!s.is_infeasible());
  bounded = !s.is_unbounded();
  
  value = CGAL::to_double(s.objective_value());
  u = CGAL::to_double(*s.variable_values_begin());
  v = CGAL::to_double(*(s.variable_values_begin()+1));
}


std::string CGAL_Solver::name() const {
  return "CGAL";
}

CGAL_Solver::~CGAL_Solver() {
  delete lp;
}
#endif
