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

#ifndef __minmax2var__glpk_solver__
#define __minmax2var__glpk_solver__

#include "external_solver.h"
#include "settings.h"

#ifdef EXTERNAL_SOLVER_GLPK

#include <glpk.h>

typedef int (*glp_solver_function)(glp_prob*, const void*);

class GLPK_Solver_base : public ExternalSolver {
protected:
  glp_prob* lp;
  void* params;
  glp_solver_function solver;
  int n;
public:
  GLPK_Solver_base();
  virtual void init(int n);
  virtual void solve(const hough_min_max::P3 points[]);
  virtual ~GLPK_Solver_base();
};

template <bool Presolve>
class GLPK_Solver_primal : public GLPK_Solver_base {
public:
  GLPK_Solver_primal() {
    // Parameters
    params = new glp_smcp;
    glp_init_smcp((glp_smcp*)params);
    ((glp_smcp*)params)->msg_lev = GLP_MSG_OFF;
    ((glp_smcp*)params)->meth = GLP_PRIMAL;
    ((glp_smcp*)params)->presolve = Presolve ? GLP_ON : GLP_OFF;
    
    // Solver
    solver = (glp_solver_function)glp_simplex;
  }
  virtual std::string name() const {
    return std::string("GLPK (Simplex, primal") + (Presolve ? ", presolve" : "") + ")";
  }
  virtual ~GLPK_Solver_primal() {
    delete (glp_smcp*)params;
  }
};

template <bool Presolve>
class GLPK_Solver_dual : public GLPK_Solver_base {
public:
  GLPK_Solver_dual() {
    // Parameters
    params = new glp_smcp;
    glp_init_smcp((glp_smcp*)params);
    ((glp_smcp*)params)->msg_lev = GLP_MSG_OFF;
    ((glp_smcp*)params)->meth = GLP_DUALP;
    ((glp_smcp*)params)->presolve = Presolve ? GLP_ON : GLP_OFF;
    
    // Solver
    solver = (glp_solver_function)glp_simplex;
  }
  virtual std::string name() const {
    return std::string("GLPK (Simplex, dual") + (Presolve ? ", presolve" : "") + ")";
  }
  virtual ~GLPK_Solver_dual()  {
    delete (glp_smcp*)params;
  }
};

// Using the interior-point method is so slow I gave up.
// The optimality test code which is inherited from
// GLPK_Solver_base is wrong for IPM, so it will fail
// an assertion. Still, it's impractically slow,
// so I didn't bother fixing it.
class GLPK_Solver_IPM : public GLPK_Solver_base {
public:
  GLPK_Solver_IPM();
  virtual std::string name() const;
  virtual ~GLPK_Solver_IPM();
};

#else
template <bool Presolve>
class GLPK_Solver_primal : public NullExternalSolver {
};
template <bool Presolve>
class GLPK_Solver_dual : public NullExternalSolver {
};
#endif

#endif /* defined(__minmax2var__glpk_solver__) */
