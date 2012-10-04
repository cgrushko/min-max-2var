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

#ifdef EXTERNAL_SOLVER_GLPK

#include "glpk-external-solver.h"
#include <sys/time.h>

GLPK_Solver_base::GLPK_Solver_base() : n(-1) {
}

void GLPK_Solver_base::init(int n) {
  if (this->n == n)
    return;

  this->n = n;
  lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MIN);

  glp_add_rows(lp, n + 1); // +1 for u+v <= 1
  glp_add_cols(lp, 3);

  glp_set_col_bnds(lp, 1, GLP_DB, 0.0, 1.0);
  glp_set_col_bnds(lp, 1, GLP_DB, 0.0, 1.0);
  
  glp_set_obj_coef(lp, 3, 1);
}

void GLPK_Solver_base::solve(const hough_min_max::P3 points[]) {
  int indices[1+ 3] = {-1, 1, 2, 3};
  double values[1+ 3];
  values[3] = -1;
  
  glp_set_col_bnds(lp, 1, GLP_LO, 0, 1);
  glp_set_col_bnds(lp, 2, GLP_LO, 0, 1);
  glp_set_col_bnds(lp, 3, GLP_LO, 0, 1);
  
  for (int i=1; i<=n; i++) {
    values[1] = points[i-1].x;
    values[2] = points[i-1].y;
    glp_set_row_bnds(lp, i, GLP_UP, 0.0, points[i-1].z);
    glp_set_mat_row(lp, i, 3, indices, values);
  }
  
  // u+v <= 1
  values[1] = 1;
  values[2] = 1;
  values[3] = 0;
  glp_set_mat_row(lp, n+1, 3, indices, values);
  glp_set_row_bnds(lp, n+1, GLP_UP, 0.0, 1);

  // Solve
  timeval t1, t2;
  gettimeofday(&t1, NULL);
  glp_std_basis(lp);
  int res = solver(lp, params);
  if (0 != res) {
    printf("GLPK solver failed with error code %d. Trying exact arithmetic.\n", res);
    res = glp_exact(lp, (const glp_smcp*)params);
    if (0 != res) {
      printf("GLPK *exact* solver failed with error code %d. Giving up.\n", res);
      abort();
    }
  }
  gettimeofday(&t2, NULL);
  time += (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  time += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

  bounded = glp_get_status(lp) == GLP_OPT;
  
  value = glp_get_obj_val(lp);
  u = glp_get_col_prim(lp, 1);
  v = glp_get_col_prim(lp, 2);
}

GLPK_Solver_base::~GLPK_Solver_base() {
  glp_delete_prob(lp);
}

// ========================================================= //

GLPK_Solver_IPM::GLPK_Solver_IPM() {
  // Parameters
  params = new glp_iptcp;
  glp_init_iptcp((glp_iptcp*)params);
  
  // Solver
  solver = (glp_solver_function)glp_interior;
}

std::string GLPK_Solver_IPM::name() const {
  return "GLPK (IPM)";
}

GLPK_Solver_IPM::~GLPK_Solver_IPM() {
  delete (glp_iptcp*)params;
}

#endif
