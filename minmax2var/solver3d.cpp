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

#include "P2.h"
#include "P3.h"
#include "performance-data.h"
#include "openmp-reductor.h"

#include "solver3d.h"
#include <limits>
#include "settings.h"
#include <stdlib.h>
#include <cmath>

typedef double REAL;

extern "C" double orient3d(const double *pa, const double *pb, const double *pc, const double *pd);

namespace hough_min_max {
  
  namespace solver_3d {
    
    V2 project_to_XY(const P3& p) {
      return V2(p.x, p.y);
    }
    
    bool colinear(const P3& p0, const P3& p1, const P3& p2) {
      return filtered_orient2d(p0.coord, p1.coord, p2.coord) == 0;
    }
    
    double getIntersect(const P3& p0, const P3& p1, const P3& p2) {
      V3 N = cross(subtract(p1, p0), subtract(p2, p0));
      return dot(N, p0) / N.z;
    }
    
    // Last coordinate of result, Nz, is guarenteed to be positive             
    V3 planeNormalFrom3Points(const P3& p0, const P3& p1, const P3& p2) {
      V3 N = cross(subtract(p1, p0), subtract(p2, p0));
      if (N.z < 0) {
        N.x = -N.x;
        N.y = -N.y;
        N.z = -N.z;
      }
      
      return N;
    }
    
    bool inDaZone(const P3& p0, const P3& p1, const P3& p2);
      
    bool assertDualSolution(const P3* ptArray, int n, const P3& p0, const P3& p1, const P3& p2) {
      release_assert(inDaZone(p0, p1, p2));
      V3 N = planeNormalFrom3Points(p0, p1, p2);
      
      for (int i=0; i<n; i++) {
        const P3& p = ptArray[i];
        release_assert(!PtIsBelowPlane(p0, p1, p2, N, p));
      }
      
      return true; // dummy return, so this function can be used in an assert
    }
    
    V2 lineNormalThatPassesThroughPointAndOrigin(const P2& p) {
      return rotate90(p);
    }
    
    bool findFeasiblePt(const P3* ptArray, int n, const P3& p0, const P3& p1, const P3** p2) {
      for (int i=0; i<n; i++) {
        if (inDaZone(p0, p1, ptArray[i]) && !colinear(p0, p1, ptArray[i])) {
          *p2 = &ptArray[i];
          return true;
        }
      }
      return false;
    }
    
    //template <int Factor>  // Factor should be either 1 or -1, for negative y, or positive y, respectively.
    bool findMostCW(const P3* ptArray, int n, P3 const ** mostCW, const P3** mostCCW, int Factor) {
      // Find most-CW and most-CCW with negative y
      if (n < 3) {
        return false;
      }
      int feasiblePoints = 0;
      P3 const* pMostCW = NULL;
      P3 const* pMostCCW = NULL;
      for (int i=0; i<n; i++) {
        const P3& p = ptArray[i];
        if (Factor * p.y <= 0) {
          feasiblePoints++;
          if (feasiblePoints == 1) {
            pMostCW = &p;
            pMostCCW = &p;
            continue;
          }
          if (cross(project_to_XY(p), project_to_XY(*pMostCW)) > 0) {
            pMostCW = &p;
          }
          if (cross(project_to_XY(p), project_to_XY(*pMostCCW)) < 0) {
            pMostCCW = &p;
          }
        }
      }
      
      if (feasiblePoints < 2) {
        return false;
      }
      
      if (Factor * pMostCW->y <= 0 && Factor * pMostCCW->y <= 0) {
        *mostCW = pMostCW;
        *mostCCW = pMostCCW;
        return true;
      }
      return false;
    }
    
    int find_quadrant(const P3& p) {
      if (p.x > 0 && p.y > 0) return 0;
      if (p.x > 0 && p.y < 0) return 3;
      if (p.x < 0 && p.y > 0) return 1;
      if (p.x < 0 && p.y < 0) return 2;
      return -1;
    }

    bool set_points_if_inDaZone(const P3& a, const P3& b, const P3& c, const P3** p0, const P3** p1, const P3** p2) {
        if (inDaZone(a, b, c)) {
          *p0 = &a;
          *p1 = &b;
          *p2 = &c;
          return true;
        }
        return false;    
    }

    bool experimental_findInitialGuess(const P3* points, int n, const P3** p0, const P3** p1, const P3** p2) {
      bool found_in_quadrant[4] = {false, false, false, false};
      const P3* quadrant_points[4] = {NULL, NULL, NULL, NULL};
      int quadrants_found = 0;
      for (int i=0; i<n && quadrants_found < 4; i++) {
        int q = find_quadrant(points[i]);
        if (q != -1 && !found_in_quadrant[q]) {
          found_in_quadrant[q] = true;
          quadrant_points[q] = &points[i];
          quadrants_found++;
        }
      }

      if (quadrants_found == 4) {
        if (set_points_if_inDaZone(*quadrant_points[0], *quadrant_points[1], *quadrant_points[2], p0, p1, p2))
          return true;
        if (set_points_if_inDaZone(*quadrant_points[0], *quadrant_points[1], *quadrant_points[3], p0, p1, p2))
          return true;
        if (set_points_if_inDaZone(*quadrant_points[0], *quadrant_points[2], *quadrant_points[3], p0, p1, p2))
          return true;
        *p0 = quadrant_points[1];
        *p1 = quadrant_points[2];
        *p2 = quadrant_points[3];
        return true;
      }
      return false;
    }

    bool findInitialGuess(const P3* ptArray, int n, const P3* * p0, const P3** p1, const P3** p2) {
      if (experimental_findInitialGuess(ptArray, n, p0, p1, p2)) {
        return true;
      }

      bool success = findMostCW(ptArray, n, p0, p1, 1);
      if (success) {
        success = findFeasiblePt(ptArray, n, **p0, **p1, p2);
        if (success) {
          internal_verification(inDaZone(**p0, **p1, **p2));
          return true;
        }
      }
      
      success = findMostCW(ptArray, n, p0, p1, -1);
      if (success) {
        success = findFeasiblePt(ptArray, n, **p0, **p1, p2);
        if (success) {
          internal_verification(inDaZone(**p0, **p1, **p2));
          return true;
        }
      }
      
      return false;
    }
    
    bool inDaZone(const P3& p0, const P3& p1, const P3& p2) {
      double a = filtered_det2(p0.coord, p1.coord);
      double b = filtered_det2(p1.coord, p2.coord);
      double c = filtered_det2(p2.coord, p0.coord);
      return ((a >= 0 && b >= 0 && c >= 0) || (a <= 0 && b <= 0 && c <= 0)) && !colinear(p0, p1, p2);
    }
    
    void letItFall_worker(const P3& new_point, const P3** p0, const P3** p1, const P3** p2, bool *changed, V3 *N) {
      if (PtIsBelowPlane(**p0, **p1, **p2, *N, new_point)) {
        bool z0 = inDaZone(**p1, **p2, new_point);
        bool z1 = inDaZone(**p0, **p2, new_point);
        bool z2 = inDaZone(**p0, **p1, new_point);
        
        if (z2 && !colinear(**p0, **p1, P3::zero)) {
          *p2 = &new_point;
        } else if (z1 && !colinear(**p0, **p2, P3::zero)) {
          *p1 = &new_point;
        } else if (z0 && !colinear(**p1, **p2, P3::zero)) {
          *p0 = &new_point;
        } else {
          if (z2) {
            *p2 = &new_point;
          } else if (z1) {
            *p1 = &new_point;
          } else if (z0) {
            *p0 = &new_point;
          } else {
            return;
          }
        }
        *changed = true;
        *N = planeNormalFrom3Points(**p0, **p1, **p2);
      }
    }
    
    bool letItFall(const P3* ptArray, int n, const P3** p0, const P3** p1, const P3** p2) {
      V3 N = planeNormalFrom3Points(**p0, **p1, **p2);
      bool changed = false;
      internal_verification(inDaZone(**p0, **p1, **p2));

      for (int i=0; i<n; i++) {
        letItFall_worker(ptArray[i], p0, p1, p2, &changed, &N);
      }
      return changed;
    }
    
    struct LetItFallBody {
      V3 N;
      bool changed;
      P3 const* p0;
      P3 const* p1;
      P3 const* p2;

      void init(const P3* p0, const P3* p1, const P3* p2) {
        this->p0 = p0; this->p1 = p1; this->p2 = p2;
        changed = false;
        N = planeNormalFrom3Points(*p0, *p1, *p2);
      }
      
      void operator()(const P3* iterator) {
        letItFall_worker(*iterator, &p0, &p1, &p2, &changed, &N);
      }
      
      void join(const LetItFallBody& rhs) {
        letItFall_worker(*rhs.p0, &p0, &p1, &p2, &changed, &N);
        letItFall_worker(*rhs.p1, &p0, &p1, &p2, &changed, &N);
        letItFall_worker(*rhs.p2, &p0, &p1, &p2, &changed, &N);
      }
    };
    
    void dualSolutionToMinMaxSolution(const P3& p0, const P3& p1, const P3& p2, double* u, double* v, double* val) {
      V3 N = cross(subtract(p1, p0), subtract(p2, p0));
      internal_verification(N.z != 0);
      *u = -N.x / N.z;
      *v = -N.y / N.z;
      *val = - dot(N, p0) / N.z;
    }

    bool solve(const P3* ptArray, int n, Result *result) {
      P3 const *p0 = NULL;
      P3 const *p1 = NULL;
      P3 const *p2 = NULL;
      if (!findInitialGuess(ptArray, n, &p0, &p1, &p2)) {
        result->u = std::numeric_limits<double>::infinity();
        result->v = std::numeric_limits<double>::infinity();
        return false;
      }

#ifdef PARALLEL
      static OpenMPReductor<const P3*, LetItFallBody> reductor;
      LetItFallBody body;
      body.init(p0, p1, p2);
#endif
      int counter = 0;
      bool c1;
      do {
        counter++;
        release_assert(counter < 1000);
#ifndef PARALLEL
        c1 = letItFall(ptArray, n, &p0, &p1, &p2);
#else
        body.changed = false;
        reductor.reduce(ptArray, ptArray + n, &body);
        c1 = body.changed;
#endif
      } while(c1);
#ifdef PARALLEL
      p0 = body.p0;
      p1 = body.p1;
      p2 = body.p2;
#endif
      
      performanceDataSingleton.add_scans_count_3d(counter);
      
      dualSolutionToMinMaxSolution(*p0, *p1, *p2, &result->u, &result->v, &result->value);
      
      internal_verification(assertDualSolution(ptArray, n, *p0, *p1, *p2));
      
      return true;
    }
    
  } // namespace solver_3d
} // namespace hough_min_max
