
/*
    Modularized simplex inverse modules - w/interface for lp_solve v5.0+
   ----------------------------------------------------------------------------------
    Author:        Michel Berkelaar (to lp_solve v3.2),
                   Kjell Eikland
    Contact:       kjell.eikland@broadpark.no
    License terms: LGPL.

    Requires:      lp_etaPFI.h, lp_lib.h, lp_colamdMDO.h

    Release notes:
    v1.0    1 September 2003    Implementation of the original lp_solve product form
                                of the inverse using a sparse eta matrix format with
                                a significant zero'th index containing the objective.
                                The new implementation includes optimal column
                                ordering at reinversion using the colamd library.
                                For lp_solve B-inverse details confer the text at the
                                beginning of lp_lib.cpp.
    v2.0    1 April 2004        Repackaged and streamlined for both internal and
                                external usage, using the new inverse/factorization
                                interface library.
    v2.0.1  23 May 2004         Moved mustrefact() function into the BFP structure.
    v2.1.0  20 July 2004        Reworked with flexible lp_solve matrix storage model.
    v2.2.0  20 February 2005    Changed to explicit OF vector mode.
                                Added updating of maximum |rhs|.
    v2.3.0  12 May 2005         Improved speed by eliminating redundant code, using
                                local variables and pointer operations.  Thanks to 
                                Henri Gourvest for suggesting several of these.
    v2.4.0  18 June 2005        Made changes to allow for "pure" factorization;
                                i.e. without the objective function included.

   ----------------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <string.h>

#include "lp_lib.h"
#include "commonlib.h"
#include "lp_etaPFI.h"

#ifdef FORTIFY
# include "lp_fortify.h"
#endif


/* Include routines common to inverse implementations; unfortunately,
   etaPFI requires the optional shared routines to be unique. */
#include "lp_BFP1.c"


/* MUST MODIFY */
char * BFP_CALLMODEL bfp_name(void)
{
#if INVERSE_ACTIVE == INVERSE_LEGACY
  return( "etaPFI v1.4" );
#else
  return( "etaPFI v2.4" );
#endif
}

/* MUST MODIFY */
MYBOOL BFP_CALLMODEL bfp_init(lprec *lp, int size, int delta, char *options)
{
  INVrec *eta;

  eta = (INVrec *) calloc(1, sizeof(*lp->invB));
  lp->invB = eta;

  if(!lp->bfp_resize(lp, ETA_START_SIZE))
     return( FALSE);

  if(delta <= 0)
    delta = lp->bfp_pivotmax(lp);
  if(!lp->bfp_pivotalloc(lp, size+delta) ||
     !lp->bfp_restart(lp))
    return( FALSE );

  lp->bfp_preparefactorization(lp);
  eta->num_refact = 0;
  eta->num_timed_refact = 0;
  eta->num_dense_refact = 0;
  eta->timed_refact = DEF_TIMEDREFACT;

  eta->statistic1 = 0;
  eta->statistic2 = delta;

  return( TRUE );
}


MYBOOL BFP_CALLMODEL bfp_restart(lprec *lp)
{
  INVrec *eta;

  eta = lp->invB;
  if(eta == NULL)
    return( FALSE );

  eta->status = BFP_STATUS_SUCCESS;
  eta->max_Bsize = 0;
  eta->max_colcount = 0;
  eta->max_etasize = 0;
  eta->num_refact = 0;
  eta->num_timed_refact = 0;
  eta->num_dense_refact = 0;
  eta->pcol = NULL;
/*  eta->set_Bidentity = FALSE; */
  eta->num_pivots = 0;

  eta->last_colcount = 0;
  eta->statistic1 = 0;
  eta->statistic2 = lp->bfp_pivotmax(lp);

  return( TRUE );
}

MYBOOL BFP_CALLMODEL bfp_implicitslack(lprec *lp)
{
  return( TRUE );
}

MYBOOL BFP_CALLMODEL bfp_pivotalloc(lprec *lp, int newsize)
{
  INVrec *eta;
  MYBOOL isfirst;

  eta = lp->invB;
  isfirst = (MYBOOL) (eta->eta_col_end == NULL);
  newsize += 1;
  if(!allocINT(lp, &eta->eta_col_end, newsize, AUTOMATIC) ||
     !allocINT(lp, &eta->eta_col_nr, newsize, AUTOMATIC))
    return( FALSE );
  if(isfirst)
    eta->eta_col_end[0] = 0;
  return( TRUE );
}


MYBOOL BFP_CALLMODEL bfp_resize(lprec *lp, int newsize)
{
  INVrec *eta;

  eta = lp->invB;
  if(eta->eta_matalloc == 0)
    eta->eta_matalloc = newsize;
  else {
    while(eta->eta_matalloc <= newsize)
      eta->eta_matalloc += eta->eta_matalloc / RESIZEFACTOR;
  }
  return( allocREAL(lp, &eta->eta_value, eta->eta_matalloc + 1, AUTOMATIC) &&
           allocINT(lp, &eta->eta_row_nr, eta->eta_matalloc + 1, AUTOMATIC) );
}


void BFP_CALLMODEL bfp_free(lprec *lp)
{
  INVrec *eta;

  eta = lp->invB;
  if(eta == NULL)
    return;

  FREE(eta->eta_value);
  FREE(eta->eta_row_nr);
  FREE(eta->eta_col_end);
  FREE(eta->eta_col_nr);
  FREE(eta);
  lp->invB = NULL;
}


/* MUST MODIFY */
MYBOOL BFP_CALLMODEL bfp_canresetbasis(lprec *lp)
{
  return( TRUE );
}

int BFP_CALLMODEL bfp_preparefactorization(lprec *lp)
{
  INVrec *eta = lp->invB;

  /* Finish any outstanding business */
  if(eta->is_dirty == AUTOMATIC)
    lp->bfp_finishfactorization(lp);

  /* Reset additional indicators */
  lp->bfp_updaterefactstats(lp);
  
  return( 0 );
}


void BFP_CALLMODEL bfp_finishfactorization(lprec *lp)
{
  INVrec *eta = lp->invB;

  /* Collect and optionally report statistics */
  if((lp->bb_totalnodes <= 1) && (lp->verbose & MSG_PERFORMANCE)) {
    REAL hold;
    int  stepsize = 5;

    /* Report size of eta */
    hold = lp->bfp_efficiency(lp);
    if((hold >= 10) && (hold >= stepsize*(1 + (int) eta->statistic1 / stepsize)))
      lp->report(lp, NORMAL, "Reduced speed with inverse density at %.1fx basis matrix.\n",
                             hold);
    SETMAX(eta->statistic1, hold);
    /* Report numeric stability */
    hold = (REAL) (lp->total_iter-lp->total_bswap)/
                  (lp->bfp_refactcount(lp, BFP_STAT_REFACT_TOTAL)+1);
    if((hold <= 30) && (hold <= stepsize*(1 + (int) eta->statistic2 / stepsize)))
      lp->report(lp, NORMAL, "Reduced numeric accuracy with %.1f pivots/refactorization.\n",
                             hold);
    SETMIN(eta->statistic2, hold);
  }
  eta->last_colcount = lp->bfp_colcount(lp);
  SETMAX(eta->max_colcount, eta->last_colcount);
  SETMAX(eta->max_etasize, lp->bfp_nonzeros(lp, FALSE));

  /* Signal that we done reinverting */
  lp->invB->is_dirty = FALSE;
  lp->clear_action(&lp->spx_action, ACTION_REINVERT | ACTION_TIMEDREINVERT);
  lp->invB->force_refact = FALSE;

  /* Store information about the current inverse */
  eta->num_pivots = 0;

}


int BFP_CALLMODEL bfp_colcount(lprec *lp)
{
  return(lp->invB->user_colcount);
}


int BFP_CALLMODEL bfp_nonzeros(lprec *lp, MYBOOL maximum)
{
  if(maximum == TRUE)
    return(lp->invB->max_etasize);
  else if(maximum == AUTOMATIC)
    return(lp->invB->max_Bsize);
  else
    return(lp->invB->eta_col_end[lp->invB->user_colcount]);
}


int BFP_CALLMODEL bfp_memallocated(lprec *lp)
{
  return(lp->invB->eta_matalloc);
}


int bfp_ETApreparepivot(lprec *lp, int row_nr, int col_nr, REAL *pcol, MYBOOL *frow)
/* Find optimal pivot row and do any column preprocessing before
   being added to the eta file -- Only used in bfp_refactorize() */
{
  int  i, rowalt, refnr;
  REAL hold, test;

#if 0
  fsolve(lp, col_nr, pcol, NULL, lp->epsmachine, 1.0, TRUE);
#else
  i = lp->get_lpcolumn(lp, col_nr, pcol, NULL, NULL);
  lp->bfp_ftran_prepare(lp, pcol, NULL);
#endif

  if(frow == NULL) {
    if(fabs(pcol[row_nr]) < lp->epsmachine)
      return( -1 );
    else
      return( row_nr );
  }

  /* Find largest coefficient row available to be pivoted */
  refnr = row_nr;
  row_nr  = lp->rows + 1;
  rowalt = 0;
  hold  = -lp->infinite;
  /* Loop over rows that currently contain slacks;
     (i.e. the positions that are not already taken by another pivot) */
#ifdef UseMarkowitzStatistic
  i = 0;
  while((i = nextActiveLink((LLrec *) frow, i)) != 0) {
#else
  for(i = 1; i <= lp->rows; i++)
    if(frow[i] == FALSE) {
#endif
      /* Get largest absolute value */
      test = fabs(pcol[i]);
      if(test > hold) {
        hold = test;
        row_nr = i;
#ifdef LegacyEtaPivotChoice
        if(hold > lp->epspivot)
          break;
#endif
      }
      if((rowalt == 0) && (test == 1))   /* Unit row */
        rowalt = i;
    }

  if((row_nr > lp->rows) || (hold < lp->epsmachine))
    row_nr = -1;

  return( row_nr );

}

int bfp_ETApivotrow(lprec *lp, int datacolumn)
/* Get pivot row for a particular data column */
{
  INVrec *eta = lp->invB;
  int    n = -1;

  if(datacolumn <= eta->last_colcount) {
    n = eta->eta_col_end[datacolumn] - 1;
    n = eta->eta_row_nr[n];
    if(n < 1) {
      lp->report(lp, CRITICAL, "bfp_ETApivotrow: Invalid pivot row identified");
      n = -1;
    }
  }
  return( n );
}


void bfp_ETAupdateCounts(lprec *lp, MYBOOL *usedpos, int rownr, int colnr)
{
  usedpos[rownr] = AUTOMATIC;
  usedpos[colnr] = AUTOMATIC;
}
void bfp_ETAreduceCounts(lprec *lp, MYBOOL *usedpos,
                         int *rownum, int *colnum, int rownr, int colnr, MYBOOL defer)
{
  int    i, j, k;
  REAL   hold;
  MATrec *mat = lp->matA;

  /* Reduce row NZ counts */
  if(rownum != NULL) {
    int *rowidx;

    hold = lp->get_OF_active(lp,  lp->rows+colnr, 1.0);
    if((usedpos[0] == FALSE) && (hold != 0))
      rownum[0]--;

    j = mat->col_end[colnr - 1];
    k = mat->col_end[colnr];
    rowidx = &COL_MAT_ROWNR(j);
    for(; j < k; j++, rowidx += matRowColStep) {
      if(usedpos[*rowidx] == FALSE)
        rownum[*rowidx]--;
    }
    rownum[rownr] = -1;
  }
  usedpos[rownr] = AUTOMATIC+defer; /* To distinguish it from basic slacks! */

  /* Reduce column NZ counts */
  if(colnum != NULL) {

    if(rownr == 0) {
      for(j = 1; j <= lp->columns; j++) {
        hold = lp->get_OF_active(lp, lp->rows + j, 1.0);
        if((usedpos[lp->rows + j] == TRUE) && (hold != 0))
          colnum[j]--;
      }
    }
    else {
      k = mat->row_end[rownr];
      for(j = mat->row_end[rownr - 1]; j < k; j++) {
        i = ROW_MAT_COLNR(j);
        if(usedpos[lp->rows + i] == TRUE)
          colnum[i]--;
      }
    }
    colnum[colnr] = -1;
  }
  usedpos[lp->rows+colnr] = AUTOMATIC+defer;

}

void bfp_ETAsimpleiteration(lprec *lp, int row_nr, int col_nr, REAL *pcol)
/* Called from bfp_refactorize() with two cases;
    pcol == NULL : A singleton row or column that allows easy elimination,
    pcol != NULL : A dense row/column combination presolved with ftran    */
{
  lp->bfp_prepareupdate(lp, row_nr, col_nr, pcol);
  lp->set_basisvar(lp, row_nr, col_nr);
  lp->bfp_finishupdate(lp, FALSE);
} /* bfp_ETAsimpleiteration */


int BFP_CALLMODEL bfp_factorize(lprec *lp, int uservars, int Bsize, MYBOOL *usedpos, MYBOOL final)
{
  REAL    *pcol, hold;
  int     *colnum, *rownum, *col, *row;
  int     *mdo = NULL;
#ifdef UseMarkowitzStatistic
  REAL    absval, testval;
  LLrec   *nextrow;
  int     ii;
#endif
  int     k, kk, i, j, jj, numit, rownr, colnr;
  int     singularities = 0;
  MATrec  *mat = lp->matA;
  REAL    *value;
  int     *rowidx;

 /* Check if there is anyting to do */
  SETMAX(lp->invB->max_Bsize, Bsize+(1+lp->rows-uservars));
  if(uservars == 0)
    return(singularities);

 /* Allocate other necessary working arrays */
  allocINT(lp,  &col,    lp->rows + 1, TRUE);
  allocINT(lp,  &row,    lp->rows + 1, TRUE);
  allocREAL(lp, &pcol,   lp->rows + 1, TRUE);
  allocINT(lp,  &rownum, lp->rows + 1, TRUE);
  allocINT(lp,  &colnum, lp->columns + 1, TRUE);

 /* Get (net/un-eliminated) row and column entry counts for basic user variables */
#ifdef ExcludeCountOrderOF
  usedpos[0] = TRUE;
#else
  if(bfp_rowoffset(lp) > 0)
    usedpos[0] = FALSE;
  else
    usedpos[0] = TRUE;
#endif
  for(j = 1; j <= lp->columns; j++) {

    /* If it is a basic user variable...*/
    if(usedpos[lp->rows+j]) {

      numit = 0;
      if(!usedpos[0]) {
        hold = lp->get_OF_active(lp, lp->rows+j, 1.0);
        if(hold != 0) {
          numit++;
          colnum[j]++;
          rownum[0]++;
        }
      }

      i = mat->col_end[j - 1];
      kk = mat->col_end[j];
      rowidx = &COL_MAT_ROWNR(i);
      value = &COL_MAT_VALUE(i);

      /* Count relevant non-zero values */
      for(; i < kk; 
          i++, numit++, rowidx += matRowColStep, value += matValueStep) {

       /* Exclude pre-eliminated rows due to basic slacks;
          this is a characteristic of the eta-model that presumes an
          initial basis with all slacks; i.e. an identity matrix */
        if(!usedpos[*rowidx]) {
          numit++;
          colnum[j]++;
          rownum[*rowidx]++;
        }
      }
    }
  }

 /* Initialize counters */
  numit = 0;
  k = 0;

#ifdef UseMarkowitzStatistic
 /* Create a linked list for the available pivot rows */
  createLink(lp->rows, &nextrow, usedpos);
#endif

 /* Loop over constraint rows, hunting for ROW singletons */
#ifdef ReprocessSingletons
Restart:
#endif
  kk = 0;
  for(i = 1; i <= lp->rows; i++) {

   /* Only process if the corresponding slack is non-basic (basis slot is available) */
    if((usedpos[i] == FALSE) && (rownum[i] == 1)) {

     /* Find first basic user column available to be pivoted */
      j = mat->row_end[i - 1];
      jj = mat->row_end[i];
      while((j < jj) && (usedpos[lp->rows + ROW_MAT_COLNR(j)] != TRUE))
        j++;
#ifdef Paranoia
      if(j >= jj)
        lp->report(lp, SEVERE, "bfp_factorize: No column to pivot in due to internal error.\n");
#endif
      colnr = ROW_MAT_COLNR(j);

     /* Reduce item counts for the selected pivot column/row */
#ifdef UseMarkowitzStatistic
      removeLink(nextrow, i);
#endif
      bfp_ETAreduceCounts(lp, usedpos, rownum, colnum, i, colnr, TRUE);

     /* Perform the pivot */
      bfp_ETAsimpleiteration(lp, i, lp->rows+colnr, NULL);
      k++;
      kk++;
    }
  }

 /* Loop over columns, hunting for COLUMN singletons;
    (only store row and column indexes for pivoting in at the end of refactorization) */
  if(k < lp->rows)
  for(colnr = 1; (k < lp->rows) && (colnr <= lp->columns); colnr++) {

   /* Only accept basic user columns not already pivoted in */
    if((usedpos[lp->rows + colnr] == TRUE) && (colnum[colnr] == 1)) {

     /* Find first available basis column to be pivoted out */
      j = mat->col_end[colnr - 1];
      jj = mat->col_end[colnr];
      rowidx = &COL_MAT_ROWNR(j);
      for(; (j < jj) && (usedpos[*rowidx] != FALSE);
          j++, rowidx += matRowColStep);
#ifdef Paranoia
      if(j >= jj)
        lp->report(lp, SEVERE, "bfp_factorize: No column to pivot out due to internal error.\n");
#endif

     /* Reduce item counts for the selected pivot column/row */
#ifdef UseMarkowitzStatistic
      removeLink(nextrow, *rowidx);
#endif
      bfp_ETAreduceCounts(lp, usedpos, rownum, colnum, *rowidx, colnr, FALSE);

      /* Store pivot information and update counters */
      col[numit] = colnr;
      row[numit] = *rowidx;
      numit++;
      k++;
      kk++;
    }
  }

 /* Check timeout and user abort again */
  if(lp->userabort(lp, -1))
    goto Cleanup;
  if(k >= lp->rows)
    goto Process;

 /* Reprocess the singleton elimination loop until supply exhausted */
#ifdef ReprocessSingletons
  if(kk > 0)
    goto Restart;
#endif

 /* Determine the number of remaining pivots and fill a minimum degree ordering column index array */
  k = lp->rows - k;
  mdo = bfp_createMDO(lp, usedpos, k, 
#ifdef UseLegacyOrdering
                      FALSE);
#else
                      TRUE);
#endif
  if(mdo == NULL)
    goto Cleanup;
  else if(mdo[0] == 0)
    goto Process;  
  kk = mdo[0];  

 /* Loop over all unprocessed basic user columns, finding an appropriate
    unused basis column position to pivot the user column into */
  for(i = 1; i <= kk; i++) {

  /* Get the entering variable */
    colnr = mdo[i];

  /* Solve for and eliminate the entering column / variable */
#ifdef UseMarkowitzStatistic
    rownr = bfp_ETApreparepivot(lp, -i, colnr, pcol, (MYBOOL*) nextrow);
#else
    rownr = bfp_ETApreparepivot(lp, -i, colnr, pcol, usedpos);
#endif
    if(rownr < 0) {
     /* This column is singular; Just let it leave the basis, making one of the
        slack variables basic in its place. (Source: Geosteiner changes!) */
      if(lp->spx_trace)
        lp->report(lp, DETAILED, "bfp_factorize: Skipped singular column %d\n", colnr);
      singularities++;
      continue;
    }
    else {
#ifdef UseMarkowitzStatistic
     /* Do a simple "local" Markowitz-based pivot selection;
        this generally reduces eta NZ count, but does not always give
        a speed improvement due to weaker numerics and added overhead
        (the numeric pivot selection tolerance limit for application
         of the Markowitz metric is 0.1, which is a typical value) */
      k = rownr;
      absval = fabs(pcol[rownr]);
      hold = 0.1*absval;
      j = 0;
      while((j = nextActiveLink(nextrow, j)) != 0) {

        /* Skip previously selected positions and focus row */
        if(j == rownr)
          continue;
        /* Pick a pivot row that is shorter than the previous best */
        jj = rownum[j];
        ii = rownum[k];
        if((jj > 0) && (jj <= ii)) {
          /* Make sure that we preserve numeric accuracy */
          testval = fabs(pcol[j]);
          if(((jj == ii) && (testval > absval)) ||
             ((jj < ii)  && (testval > hold))) {
            k = j;
            absval = testval;
          }
        }
      }
      rownr = k;
      bfp_ETAsimpleiteration(lp, rownr, colnr, pcol);
    }

   /* Reduce item counts for the selected pivot column/row */
    if(rownr > 0)
      removeLink(nextrow, rownr);
    bfp_ETAreduceCounts(lp, usedpos, rownum, NULL, rownr, colnr-lp->rows, FALSE);

#else
      bfp_ETAsimpleiteration(lp, rownr, colnr, pcol);
    }

   /* Update occupancy states for pivot column/row */
    bfp_ETAupdateCounts(lp, usedpos, rownr, colnr);
#endif

   /* Check timeout and user abort again */
    if(lp->userabort(lp, -1))
      goto Cleanup;
  }

 /* Perform pivoting of the singleton columns stored above */
Process:
  for(i = numit - 1; i >= 0; i--) {
    colnr = col[i];
    rownr = row[i];
    bfp_ETAsimpleiteration(lp, rownr, lp->rows+colnr, NULL);
  }

 /* Finally, wrap up the refactorization */
Cleanup:
  FREE(mdo);
  if(pcol != NULL) {
    FREE(pcol);
    FREE(col);
    FREE(row);
    FREE(rownum);
    FREE(colnum);
  }
#ifdef UseMarkowitzStatistic
  freeLink(&nextrow);
#endif

  lp->invB->num_singular += singularities;
  return(singularities);
}


LREAL BFP_CALLMODEL bfp_prepareupdate(lprec *lp, int row_nr, int col_nr, REAL *pcol)
/* Was condensecol() in versions of lp_solve before 4.0.1.8 - KE */
{
  int    i, j, jj, k, colusr, elnr, min_size;
  LREAL  pivValue;
  INVrec *eta = lp->invB;
  MATrec *mat = lp->matA;

  eta->pcol = pcol;

  elnr = eta->eta_col_end[eta->user_colcount];
  pivValue = 0;

#ifdef Paranoia
  if(row_nr < 1 || col_nr < 1)
    lp->report(lp, CRITICAL, "bfp_prepareupdate: Invalid row/column combination specified");
#endif

  min_size = elnr + lp->rows + 2;
  if(min_size >= eta->eta_matalloc) /* maximum local growth of Eta */
    lp->bfp_resize(lp, min_size);

  /* Fill the eta-column from the A matrix */
  if(pcol == NULL) {
    REAL   *value;
    int    *rowidx;

    k = 0;
    colusr = col_nr - lp->rows;

    /* Handle phase 1 objective function adjustments */
    if(bfp_rowoffset(lp) > 0) {
      pivValue = lp->get_OF_active(lp, col_nr, 1.0);
      if(pivValue != 0) {
        eta->eta_row_nr[elnr] = 0;
        eta->eta_value[elnr] = pivValue;
        elnr++;
      }
    }

    j = mat->col_end[colusr - 1];
    jj = mat->col_end[colusr];
    rowidx = &COL_MAT_ROWNR(j);
    value  = &COL_MAT_VALUE(j);
    
    for(; j < jj; 
        j++, rowidx += matRowColStep, value += matValueStep) {

      /* The pivot row is stored at the end of the column; store its index */
      if(*rowidx == row_nr) {
        k = j;
        continue;
      }

      /* Append other non-zero column values */
      eta->eta_row_nr[elnr] = *rowidx;
      eta->eta_value[elnr] = *value;
      elnr++;
    }
    eta->eta_row_nr[elnr] = row_nr;
    eta->eta_value[elnr] = COL_MAT_VALUE(k);
    elnr++;
  }

  /* Fill the eta-colum from a dense, ftran-preprocessed column
     where the data has been retrieved with obtain_column().
     Note that phase 1 objective function adjustments are done in
     obtain_column/expand_column, called in fsolve() */
  else {
    int  r = lp->rows,                  *etarownr = &eta->eta_row_nr[elnr];
    REAL *p = pcol+1-bfp_rowoffset(lp), *etavalue = &eta->eta_value[elnr];

    for(i = 1 - bfp_rowoffset(lp); i <= r; i++, p++) {
      if((i != row_nr) && (*p != 0)) {
        *etarownr = i;
        etarownr++;
        *etavalue = *p;
        etavalue++;
        elnr++;
      }
    }
    elnr++;
    *etarownr = row_nr;
    *etavalue = pcol[row_nr];
    pivValue = *etavalue;
  }
  i = eta->user_colcount;
  eta->eta_col_nr[i] = col_nr;
  eta->eta_col_end[i+1] = elnr;
  SETMAX(eta->last_colcount, i+1);

  /* Set completion status; but hold if we are reinverting */
  if(eta->is_dirty != AUTOMATIC)
    lp->invB->is_dirty = TRUE;

  return(pivValue);
}


REAL BFP_CALLMODEL bfp_pivotRHS(lprec *lp, LREAL theta, REAL *pcol)
/* Was rhsmincol(), ie. "rhs minus column" in versions of lp_solve before 4.0.1.8 - KE */
{
  int    i;
  LREAL  f = 0;

 /* Round net RHS values (should not be smaller than the factor used in recompute_solution) */
  register LREAL  roundzero = lp->epsvalue;
  register LREAL  rhsmax = 0;

  if(pcol == NULL) {
    int    j, k, *nv;
    REAL   *vv;
    INVrec *eta = lp->invB;

    j = eta->eta_col_end[eta->user_colcount];
    k = eta->eta_col_end[eta->user_colcount + 1];
    for(i = j, nv = eta->eta_row_nr+j, vv = eta->eta_value+j;
        i < k; i++, nv++, vv++) {
      f = lp->rhs[*nv] - theta * (*vv);
      my_roundzero(f, roundzero);
      SETMAX(rhsmax, fabs(f));
      lp->rhs[*nv] = f;
    }
    f = eta->eta_value[k - 1];
  }
  else {
    int       ibase = 0; /*1 - bfp_rowoffset(lp);*/
    register LREAL *rhs;

    for(i = ibase, rhs = lp->rhs + ibase, pcol += ibase; 
        i<= lp->rows; i++, rhs++, pcol++) {
      if(*pcol == 0)
        continue;
      *rhs -= theta * (*pcol);
      my_roundzero((*rhs), roundzero);
      SETMAX(rhsmax, fabs(*rhs));
    }
    f = 0;
  }
  lp->rhsmax = rhsmax;

  return( f );

}


MYBOOL BFP_CALLMODEL bfp_finishupdate(lprec *lp, MYBOOL changesign)
/* Was addetacol() in versions of lp_solve before 4.0.1.8 - KE */
{
  int    i, j, k;
  REAL   theta, *value;
  INVrec *eta = lp->invB;

  /* Check if a data column has been loaded by bfp_putcolumn */
  if(!eta->is_dirty)
    return( FALSE );

 /* Do fast eta transformation (formally "u") */
  j = eta->eta_col_end[eta->user_colcount];
  eta->user_colcount++;
  k = eta->eta_col_end[eta->user_colcount] - 1;

 /* Handle following cases: 1) changesign == FALSE, theta == -1    -> Drop straight through
                            2)            == FALSE        != -1
                            3)            == TRUE         ==  1    -> Sign change only
                            4)            == TRUE         !=  1 */
 /* Amended old style */
  if(changesign)
    for(i = j, value = eta->eta_value+j; i <= k; i++, value++)
      *value = -(*value);

  value  = eta->eta_value+k;
  theta  = 1.0 / (*value);
  *value = theta;
  if(fabs(theta+1) > EPS_ETAMACHINE) {
    theta = -theta;
    while(k > j) {
      k--;
      value--;
      *value *= theta;
    }
  }

  eta->num_pivots++;

  /* Reset indicators; treat reinversion specially */
  if(eta->is_dirty != AUTOMATIC) {
    eta->is_dirty = FALSE;
  }

  return( TRUE );

} /* bfp_finishupdate */


void BFP_CALLMODEL bfp_ftran_normal(lprec *lp, REAL *pcol, int *nzidx)
/* Note that ftran does not "expand" B indexes to the actual basis
   column, and that both the input and output range is [0..rows] */
{
  INVrec *eta = lp->invB;
  LREAL  theta, *vcol;
  REAL   *valuep, *values = eta->eta_value;
  int    i, j, k, r, *rowp, *rownr = eta->eta_row_nr,
         ucc = eta->user_colcount, *col_end = eta->eta_col_end;

  /* Initialize case where long doubles may be used */
  if(sizeof(LREAL) == sizeof(REAL))
    vcol = pcol;
  else {
    r = lp->rows;
    if(!allocLREAL(lp, &vcol, r + 1, FALSE))
      return;
    for(i = 1 - bfp_rowoffset(lp); i <= r; i++)
      vcol[i] = pcol[i];
  }

  /* Loop forward over the eta matrices */
  for(i = 1; i <= ucc; i++) {
    j = *col_end;
    col_end++;
    k = (*col_end) - 1;
    r = rownr[k];
    theta = vcol[r];
    if(theta != 0) {

      /* CPU intensive loop, let's do pointer arithmetic */
      for(rowp = rownr + j, valuep = values + j;
          j < k; j++, rowp++, valuep++) {
        vcol[*rowp] += theta * (*valuep);
      }
      vcol[r] *= (*valuep);
    }
  }

  /* Finalize case where long doubles may be used */
  if(sizeof(LREAL) != sizeof(REAL)) {
    r = lp->rows;
    for(i = 1 - bfp_rowoffset(lp); i <= r; i++)
      pcol[i] = vcol[i];
    FREE(vcol);
  }

} /* ftran */

void BFP_CALLMODEL bfp_ftran_prepare(lprec *lp, REAL *pcol, int *nzidx)
/* Does nothing particular in the etaPFI version of the inverse;
   in other versions it is used to additionally store necessary data
   to prepare for pivoting / update of the inverse */
{
  lp->bfp_ftran_normal(lp, pcol, nzidx);
}


void BFP_CALLMODEL bfp_btran_normal(lprec *lp, REAL *prow, int *nzidx)
/* Note that btran does not "expand" B indexes to the actual basis
   column, and that both the input and output range is [0..rows] */
{
  int      i, jb, je, k, *rowp;
  LREAL    fmax;
  REGISTER LREAL f;
  REAL     *valuep;
  INVrec   *eta = lp->invB;

  /* Loop backward over the eta matrices */
  k = 0;
  fmax = 0;
  for(i = eta->user_colcount; i >= 1; i--) {
    f = 0;
    jb = eta->eta_col_end[i-1];
    je = eta->eta_col_end[i] - 1;
    for(rowp = eta->eta_row_nr + jb, valuep = eta->eta_value + jb;
        jb <= je; jb++, rowp++, valuep++) {
      f += prow[(*rowp)] * (*valuep);
    }
    SETMAX(fmax, fabs(f));
#ifndef EtaBtranRoundRelative
    my_roundzero(f, EPS_ETAPIVOT);
#endif
    jb = eta->eta_row_nr[je];
    prow[jb] = f;

  }

  /* Help numeric accuracy in case we have not scaled to equilibration */
#ifdef EtaBtranRoundRelative
  if(fmax > 0) {
    int ie = 1 - bfp_rowoffset(lp);
    fmax *= EPS_ETAPIVOT;
    for(i = lp->rows, valuep = prow; i >= ie; i--, valuep++)
      if(fabs(*valuep) < fmax)
        *valuep = 0;
  }
#endif

} /* btran */


void BFP_CALLMODEL bfp_btran_double(lprec *lp, REAL *prow, int *pnzidx, REAL *drow, int *dnzidx)
{
  int      i, j, k, *rowp;
  LREAL    dmax, fmax;
  REGISTER LREAL  d, f;
  REAL     *valuep;
  INVrec   *eta = lp->invB;

  fmax = 0;
  dmax = 0;
  for(i = eta->user_colcount; i >= 1; i--) {
    d = 0;
    f = 0;
    k = eta->eta_col_end[i] - 1;
    j = eta->eta_col_end[i - 1];

  /* This is one of the loops where the program consumes a lot of CPU time;
     let's help the compiler by doing some pointer arithmetic instead of array indexing */
    for(rowp = eta->eta_row_nr + j, valuep = eta->eta_value + j;
        j <= k; j++, rowp++, valuep++) {
      f += prow[(*rowp)] * (*valuep);
      d += drow[(*rowp)] * (*valuep);
    }
    SETMAX(fmax, fabs(f));
    SETMAX(dmax, fabs(d));

    j = eta->eta_row_nr[k];

#ifndef EtaBtranRoundRelative
    my_roundzero(f, EPS_ETAPIVOT);
    my_roundzero(d, EPS_ETAPIVOT);
#endif
    prow[j] = f;
    drow[j] = d;
  }

  /* Help numeric accuracy in case we have not scaled to equilibration */
#ifdef EtaBtranRoundRelative
  if((fmax > 0) && (dmax > 0)) {
    int  ie = 1 - bfp_rowoffset(lp);
    REAL *valued;
    fmax *= EPS_ETAPIVOT;
    dmax *= EPS_ETAPIVOT;
    for(i = lp->rows, valuep = prow, valued = drow; i >= ie; i--, valuep++, valued++) {
      if(fabs(*valuep) < fmax)
        *valuep = 0;
      if(fabs(*valued) < dmax)
        *valued = 0;
    }
  }
#endif

}

/* MUST MODIFY - Routine to find maximum rank of equality constraints */
int BFP_CALLMODEL bfp_findredundant(lprec *lp, int items, getcolumnex_func cb, int *maprow, int *mapcol)
{
  /* Are we capable of finding redundancy with this BFP? */
  if((maprow == NULL) && (mapcol == NULL))
    return( 0 );

  /* If so, process */
  return( 0 );
}
