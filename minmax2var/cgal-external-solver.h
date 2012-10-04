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

#ifndef __minmax2var__cgal_external_solver__
#define __minmax2var__cgal_external_solver__

#include "external_solver.h"
#include "settings.h"

#ifdef EXTERNAL_SOLVER_CGAL
#define CGAL_QP_NO_ASSERTIONS
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <vector>

using std::vector;

// choose exact integral type
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpq.h>
typedef CGAL::Gmpzf ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif

typedef CGAL::Nonnegative_linear_program_from_iterators
<double**,                                             // for A
double*,                                               // for b
CGAL::Const_oneset_iterator<CGAL::Comparison_result>,  // for r
double*>                                               // for c
Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

class CGAL_Solver : public ExternalSolver {
  Program* lp;
  int n;
  vector<double> A_vectors[3];
  double *A[3];
  vector<double> b;
  double c[3];
public:
  CGAL_Solver();
  virtual void init(int n);
  virtual void solve(const hough_min_max::P3 points[]);
  virtual std::string name() const;
  virtual ~CGAL_Solver();
};

#else
typedef NullExternalSolver CGAL_Solver;
#endif

#endif /* defined(__minmax2var__cgal_external_solver__) */
