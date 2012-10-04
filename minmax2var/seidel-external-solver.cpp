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
//  Created by Carmi Grushko on 9/6/12.

#include "settings.h"

#ifdef EXTERNAL_SOLVER_SEIDEL
#include "seidel-external-solver.h"
#include <sys/time.h>
#include <assert.h>

Seidel_Solver::Seidel_Solver() : n(-1), halves_count(-1) {
  n_vec[0] = 0;
  n_vec[1] = 0;
  n_vec[2] = 1;
  n_vec[3] = 0;
  d_vec[0] = 0;
  d_vec[1] = 0;
  d_vec[2] = 0;
  d_vec[3] = 1;
}

void Seidel_Solver::init(int n) {
  if (n == this->n)
    return;

  this->n = n;
  d = 3;
  halves_count = n + numberOfExtraHalves();
  perm = new int[halves_count-1];
  next = new int[halves_count];
  prev = new int[halves_count];
  work = new FLOAT[(halves_count+3)*(d+2)*(d-1)/2];
  
  halves = new FLOAT[halves_count*(d+1)];
  input_halves = halves;
  addConstraint(0, 0, 0, 1);
  addConstraint(1, 0, 0, 0);
  addConstraint(0, 1, 0, 0);
  addConstraint(-1, -1, 0, 1);
  assert(input_halves == halves + numberOfExtraHalves() * (d+1));
}

int Seidel_Solver::numberOfExtraHalves() {
  // w>0, x>0, y>0, x+y<1
  return 4;
}

void Seidel_Solver::addConstraint(FLOAT a, FLOAT b, FLOAT c, FLOAT d) {
  input_halves[0] = a;
  input_halves[1] = b;
  input_halves[2] = c;
  input_halves[3] = d;
  input_halves += 4;
}

void Seidel_Solver::preparePermutationArrays() {
  randperm(halves_count-1,perm);

  /* previous to 0 should never be used */
  prev[0] = 1234567890;
  /* link the zero position in at the beginning */
  next[0] = perm[0]+1;
  prev[perm[0]+1] = 0;
  /* link the other planes */
  for(int i=0; i<halves_count-2; i++) {
    next[perm[i]+1] = perm[i+1]+1;
    prev[perm[i+1]+1] = perm[i]+1;
  }
  /* flag the last plane */
  next[perm[halves_count-2]+1] = halves_count;
}

void Seidel_Solver::fillHalvesArray(const hough_min_max::P3 points[]) {
  FLOAT *it = input_halves;
  for (int i=0; i<n; i++) {
    *it++ = -points[i].x;
    *it++ = -points[i].y;
    *it++ = 1;
    *it++ = -(-points[i].z); // points[] is in dual space, negate z coordinate
  }
}

void Seidel_Solver::solve(const hough_min_max::P3 points[]) {
  timeval t1, t2;
  gettimeofday(&t1, NULL);
  fillHalvesArray(points);
  preparePermutationArrays();
  bounded = linprog(halves, 0, halves_count, n_vec, d_vec, d,
                    opt, work, next, prev, halves_count) == MINIMUM;
  gettimeofday(&t2, NULL);
  time += (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  time += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

  u = opt[0] / opt[3];
  v = opt[1] / opt[3];
  value = opt[2] / opt[3];
}

std::string Seidel_Solver::name() const {
  return "Seidel";
}

Seidel_Solver::~Seidel_Solver() {
  delete[] perm;
  delete[] next;
  delete[] prev;
  delete[] work;
}

#endif
