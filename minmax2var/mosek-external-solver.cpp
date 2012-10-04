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
//  Created by Carmi Grushko on 9/9/12.

#include "settings.h"

#ifdef EXTERNAL_SOLVER_MOSEK

#include "mosek-external-solver.h"
#include <sys/time.h>

// Uncomment to log MOSEK operations to stdout.
// #define EXTERNAL_SOLVER_MOSEK_LOGGING

#ifdef EXTERNAL_SOLVER_MOSEK_LOGGING
static void MSKAPI printstr(void *handle,
                            char str[])
{
  printf("%s",str);
} /* printstr */
#endif

MOSEK_Environment_Singleton MOSEK_Environment_singleton;

MOSEK_Solver::MOSEK_Solver(MOSEK_Solver_Type::SolverType solverType) : n(-1) {
  release_assert(MSK_RES_OK == MSK_maketask (MOSEK_Environment_singleton.env, 0,0, &task));

  MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
  MSK_append(task, MSK_ACC_VAR, 3);

  // Variable bounds (0 <= u, v <= 1, 0 <= t)
  MSK_putbound(task, MSK_ACC_VAR, MOSEK_Solver_Variables::u, MSK_BK_RA, 0, 1);
  MSK_putbound(task, MSK_ACC_VAR, MOSEK_Solver_Variables::v, MSK_BK_RA, 0, 1);
  MSK_putbound(task, MSK_ACC_VAR, MOSEK_Solver_Variables::t, MSK_BK_LO, 0, MSK_INFINITY);

  // Objective: 0*u + 0*v + 1*t
  MSK_putcj(task, MOSEK_Solver_Variables::u, 0);
  MSK_putcj(task, MOSEK_Solver_Variables::v, 0);
  MSK_putcj(task, MOSEK_Solver_Variables::t, 1);

  // Turn off presolve steps (experiments hint it's faster like this)
  MSK_putintparam(task, MSK_IPAR_PRESOLVE_ELIMINATOR_USE, MSK_OFF);
  MSK_putintparam(task, MSK_IPAR_PRESOLVE_LINDEP_USE, MSK_OFF);
  
  // Select solver (IPM / Primal Simplex / Dual Simplex)
  mosek_solver = MSK_OPTIMIZER_BEGIN;
  switch (solverType) {
    case MOSEK_Solver_Type::IPM:
      mosek_solver = MSK_OPTIMIZER_INTPNT;
      break;
    case MOSEK_Solver_Type::PrimalSimplex:
      mosek_solver = MSK_OPTIMIZER_PRIMAL_SIMPLEX;
      break;
    case MOSEK_Solver_Type::DualSimplex:
      mosek_solver = MSK_OPTIMIZER_DUAL_SIMPLEX;
      break;
    case MOSEK_Solver_Type::PrimalDualSimplex:
      mosek_solver = MSK_OPTIMIZER_PRIMAL_DUAL_SIMPLEX;
      break;
  }
  MSK_putintparam(task, MSK_IPAR_OPTIMIZER, mosek_solver);
  
#ifdef EXTERNAL_SOLVER_MOSEK_LOGGING
  MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);
#endif
}

void MOSEK_Solver::init(int n) {
  if (this->n == n)
    return;
 
  this->n = n;
  
  // n+1 constraints, 3 variables (u, v, t)
  MSK_append(task, MSK_ACC_CON, n + 1); // n constraints + 1 for u+v <=1
  
  // Set u+v <= 1 constraint
  MSK_putaij(task, 0, MOSEK_Solver_Variables::u, 1);
  MSK_putaij(task, 0, MOSEK_Solver_Variables::v, 1);
  MSK_putaij(task, 0, MOSEK_Solver_Variables::t, 0);
  MSK_putbound(task, MSK_ACC_CON, 0, MSK_BK_RA, 0, 1);
}

void MOSEK_Solver::solve(const hough_min_max::P3 points[]) {
  for (int i=1; i<=n; i++) {
    MSK_putbound(task, MSK_ACC_CON, i, MSK_BK_UP, 0, points[i-1].z);
    MSK_putaij(task, i, MOSEK_Solver_Variables::u, points[i-1].x);
    MSK_putaij(task, i, MOSEK_Solver_Variables::v, points[i-1].y);
    MSK_putaij(task, i, MOSEK_Solver_Variables::t, -1);
  }
  
  // Solve
  timeval t1, t2;
  gettimeofday(&t1, NULL);
  MSK_optimize(task);
  gettimeofday(&t2, NULL);
  time += (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  time += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
  
  MSKprostae prosta;
  MSKsolstae solsta;
  MSK_getsolutionstatus(task, MSK_SOL_BAS, &prosta, &solsta);
  bounded = solsta == MSK_SOL_STA_OPTIMAL;

#ifdef EXTERNAL_SOLVER_MOSEK_LOGGING
  MSK_solutionsummary(task, MSK_STREAM_LOG);
#endif
  
  double vars[3];
  MSK_getsolutionslice(task, MSK_SOL_BAS, MSK_SOL_ITEM_XX, 0, 3, vars);
  u = vars[0];
  v = vars[1];
  value = vars[2];
}

MOSEK_Solver::~MOSEK_Solver() {
  MSK_deletetask(&task);
}

std::string solver_type_to_string(MSKoptimizertype_enum solverType) {
  switch (solverType) {
    case MSK_OPTIMIZER_INTPNT:
      return "IPM";
    case MSK_OPTIMIZER_PRIMAL_DUAL_SIMPLEX:
      return "Simplex, Prima-dual";
    case MSK_OPTIMIZER_DUAL_SIMPLEX:
      return "Simplex, Dual";
    case MSK_OPTIMIZER_PRIMAL_SIMPLEX:
      return "Simplex, Primal";
    default:
      abort();
  }
  abort();
}

std::string MOSEK_Solver::name() const {
  return "MOSEK (" + solver_type_to_string(mosek_solver) + ")";
}

#endif

