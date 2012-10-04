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
//  Created by Carmi Grushko on 6/21/12.

#include "settings.h"

#ifdef EXTERNAL_SOLVER_LPSOLVE55

#include "lpsolve55-external-solver.h"
#include <sys/time.h>
#include <string>

lpsolve55_Solver::lpsolve55_Solver() : n(-1) {
}

void setConstraint(lprec* lp, int row, double a, double b, double c) {
  lp->orig_rhs[1+ row+1] = c;
  lp->matA->col_mat_value[1+ row] = a; // +1 because of the u1+u2<=1 constraint
  lp->matA->col_mat_value[1+ row + lp->matA->rows] = b; // +1 because of the u1+u2<=1 constraint  
}

static lprec* setup_lp(int N) {
  enum { u1 = 1, u2, t };
  int colno[3];
  double row[4] = {0, 3, 2, 1};
  
  lprec *lp = make_lp( N-1 +1 +1, 3 ); 
  
#ifdef _DEBUG
  set_col_name( lp, u1, "u1" );
  set_col_name( lp, u2, "u2" );
  set_col_name( lp, t, "t" );
#endif
  
  // Set objective function to be just t
  row[0] = 1;
  colno[0] = t;
  set_obj_fnex( lp, 1, row, colno );
  
  // u1+u2 <= 1
  row[0] = row[1] = 1;
  colno[0] = u1;
  colno[1] = u2;
  set_rh( lp, 1, 1 );
  set_constr_type( lp, 1, LE );
  set_rowex( lp, 1, 2, row, colno );
  
  row[0] = row[1] = 1;
  row[3] = -1; 
  
  // The rest of the constraints
  for( int i=2; i<=lp->rows; i++ )
  {
    set_row( lp, i, row );
    set_constr_type( lp, i, LE );
  }
  //  set_unbounded(lp, 1);
  //  set_bounds(lp, 1, 0, 1);
  //  set_unbounded(lp, 2);
  set_unbounded(lp, 3);
  
  set_verbose( lp, IMPORTANT );
  set_scaling( lp, SCALE_NONE );  // Experiments show a considerable speedup when disablding scaling
  
  return lp;
}

lprec* setup_unconstrained_lp(int N) {
  enum { u1 = 1, u2, t };
  int colno[3];
  double row[4] = {0, 3, 2, 1};
  
  lprec *lp = make_lp( N-1 +1, 3 ); 
 
  // Set objective function to be just t
  row[0] = 1;
  colno[0] = t;
  set_obj_fnex( lp, 1, row, colno );
  
  row[0] = row[1] = 1;
  row[3] = -1; 
  
  // The rest of the constraints
  for( int i=1; i<=lp->rows; i++ )
  {
    set_row( lp, i, row );
    set_constr_type( lp, i, LE );
  }
  set_unbounded(lp, 1);
  set_unbounded(lp, 2);
  set_unbounded(lp, 3);
  
  set_verbose( lp, IMPORTANT );
  set_scaling( lp, SCALE_NONE );  // Experiments show a considerable speedup when disablding scaling
  
  return lp;
}

void lpsolve55_Solver::init(int n) {
  if (n == this->n)
    return;
  this->n = n;
  lp = setup_lp(n);
}

void lpsolve55_Solver::solve(const hough_min_max::P3 points[]) {
  for (int i=0; i<n; i++) {
    setConstraint(lp, i, points[i].x, points[i].y, points[i].z);
  }

  timeval t1, t2;
  gettimeofday(&t1, NULL);
  default_basis(lp);
  bounded = ::solve(lp) == 0;
  gettimeofday(&t2, NULL);
  time += (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  time += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
 
  REAL var[3];
  get_variables(lp, var);
  u = var[0];
  v = var[1];
  value = var[2];
}

std::string lpsolve55_Solver::name() const {
  return "lpsolve55";
}

lpsolve55_Solver::~lpsolve55_Solver() {
  delete_lp(lp);
}

// ============================================================ //

#ifdef ETA_PFI_LIBRARY
static const char* etaPFI_lib = ETA_PFI_LIBRARY;

void lpsolve55_Solver_etaPFI::init(int n) {
  if (this->n == n)
    return;
  lpsolve55_Solver::init(n);
  release_assert(set_BFP(lp, (char*)etaPFI_lib));
}

void lpsolve55_Solver_etaPFI::solve(const hough_min_max::P3 points[]) {
  for (int i=0; i<n; i++) {
    setConstraint(lp, i, points[i].x, points[i].y, points[i].z);
  }
  
  timeval t1, t2;
  gettimeofday(&t1, NULL);
  default_basis(lp);
  bounded = true;
  int res = 0;
  if ((res = ::solve(lp)) != 0) {
    printf("etaPFI failed to solve problem; falling back to LUSOL.\n");
    delete_lp(lp);
    lp = setup_lp(n);

    for (int i=0; i<n; i++) {
      setConstraint(lp, i, points[i].x, points[i].y, points[i].z);
    }
    
    release_assert(::solve(lp) == 0);
    release_assert(set_BFP(lp, (char*)etaPFI_lib));
  }
  gettimeofday(&t2, NULL);
  time += (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  time += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
  
  REAL var[3];
  get_variables(lp, var);
  u = var[0];
  v = var[1];
  value = var[2];
}

std::string lpsolve55_Solver_etaPFI::name() const {
  return "lpsolve55 (etaPFI)";
}
#endif

#endif
