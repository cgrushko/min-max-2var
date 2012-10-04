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

#ifndef __minmax2var__mosek_external_solver__
#define __minmax2var__mosek_external_solver__

#include "external_solver.h"
#include "settings.h"

namespace MOSEK_Solver_Type {
  enum SolverType {IPM, PrimalSimplex, DualSimplex, PrimalDualSimplex};
}

#ifdef EXTERNAL_SOLVER_MOSEK

#include "mosek.h"

class MOSEK_Environment_Singleton {
public:
  MSKenv_t env;
  MOSEK_Environment_Singleton() {
    release_assert(MSK_RES_OK == MSK_makeenv(&env, NULL, NULL, NULL, NULL));
    MSK_initenv (env);
  }
  ~MOSEK_Environment_Singleton() {
    MSK_deleteenv(&env);
  }
};

extern MOSEK_Environment_Singleton MOSEK_Environment_singleton;

// ========================================================= //

namespace MOSEK_Solver_Variables {
  enum Variables {u, v, t};
}

class MOSEK_Solver : public ExternalSolver {
protected:
  MSKtask_t task;
  int n;
  MSKoptimizertype_enum mosek_solver;
public:
  MOSEK_Solver(MOSEK_Solver_Type::SolverType solverType);
  virtual void init(int n);
  virtual void solve(const hough_min_max::P3 points[]);
  virtual std::string name() const;
  virtual ~MOSEK_Solver();
};
#else
class MOSEK_Solver : public NullExternalSolver {
public:
  MOSEK_Solver(MOSEK_Solver_Type::SolverType solverType) {
    // Do nothing
  }
};
#endif

#endif /* defined(__minmax2var__mosek_external_solver__) */
