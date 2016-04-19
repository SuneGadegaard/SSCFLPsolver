/********************************************************************************

FILE      : vico.c

VERSION   : 2.9.0
CHANGE_LOG: DATE              VER.-No.  CHANGES MADE
            ---------------------------------------------------------------------
            May  12-19, 2003    1.0.0   first implementation
	    July     7, 2003    1.1.0   bug in VICuflohi removed
	    Aug      6, 2003    2.0.0   routine VIClci and VICkconf added
	    Aug      7, 2003    2.1.0   routine VICgapfn added
	    Aug     15, 2003    2.2.0   routine VICuflcov added
	    Aug     22, 2003    2.3.0   VIClci, VICkconf modified to ease handling
	    Aug     23, 2003    2.3.1   bugs in VICuflohi removed
	                                (shortest path computation now based on
					 FIFO, error in cut existence check removed)
	    Aug     26, 2003    2.3.2   minor changes, some "statics" inserted
	    Aug     27, 2003    2.3.3   minor changes: routines VICsearchcut,
	                                VICclearlst added
	    Feb     20, 2004    2.4.0	routine VICeci added
	    Feb     24, 2004    2.4.1   routine VICsort added
	    Feb     25, 2004    2.5.0   routine VICkpreduce added
	    March   11, 2004    2.5.1   small bug in VIClci removed
	    May     10, 2005    2.5.2   bug in VICkconf "(if card==n ) ... return" removed
	    August   1, 2005    2.5.3   bug in kirestore() removed (all coefficients of
	                                a generated cut were modified if at least one
					variable was inverted!)
            May 30    , 2007    2.6.0   Started with implementing Fenchel cuts based on
	                                single-node flow structures
            June 5    , 2007    2.6.1   Different subgradient strategies for Fenchel cut
	                                generation implemented.
            June 14   , 2007    2.6.2   Single-node flow problems arising in Fenchel cut
	                                generation may optionally be solved using different
					methods for this problem
            June 15   , 2007    2.6.3   Additional parameter "Freq" for Fenchel
                                        cut generation introduced
            June 22   , 2007    2.7.0   Additional parameter "Reduce" for Fenchel cut generation
                                        introduced. If Reduce=0 all "zero arcs" are removed and
                                        the inequality is generated for the reduced polytope
            Aug 20    , 2007    2.7.0   small bug in vicsnffenchel removed: capacities are now
	                                allowed to be zero (variable can then be ignored)
            Jan  9    , 2008    2.8.0   start to implement routine for exact knapsack separation
	    Jun 18    , 2008    2.8.1   Some bugs removed in VICkplift and VICkpsep
	    Jun 19    , 2008    2.8.2   Bug remove in VICkplift ( (t-pi) could be negative)
            Aug 20    , 2012    2.8.3   Adjusted the uplifting in VICkplift and included possibility
	                                to exclusively fix variables of zero LP value when
					defining the reduced knapsack polytope in the exact knapsack
					separation procedure
            Dec 18    , 2012    2.9.0   Inclusion of a procedure suggested by Kaparis and Letchford (2010)
                                        to separate extended cover inequalities
LANGUAGE : c
AUTHOR   : Andreas Klose

SUBJECT  : Valid Inequalities for selected Combinatorial Optimization problems
           (separation routines)

EXTERNALS: (in addition to routines provided by compiler libraries )

  ----------------------------------------------------------
  Function         File       Description
  ----------------------------------------------------------
  combo            combo.h    solves binary knapsack problems
  CPXcreateprob    cplex.h    create LP problem
  CPXcopylp        cplex.h    copy LP data to CPLEX
  CPXfreeprob      cplex.h    free LP object
  CPXgetmipobjval  cplex.h    get MIP objective value
  CPXgetobjval     cplex.h    get LP objective value
  CPXmipopt        cplex.h    optimize MIP
  CPXoptimize      cplex.h    optimize LP
  ssfctp_mt1       ssfctp.h   solves single-node flow problem
  ----------------------------------------------------------

  ----------------------------------------------------------
  Typdefs          File       Description
  ----------------------------------------------------------
  item             combo.h    item of a knapsack
  CPXENVptr        cplex.h    pointer to CPLEX's environment
  CPXLPptr         cplex.h    pointer to an LP object

FUNCTIONS PROVIDED BY THIS FILE:

A. Functions which will be exported ("SCOPE: Export")
   VICaddtolst( )
   VICallocut( )
   VICcflfc( )
   VICcflsmi( )
   VICclearlst( )
   VICeci( )
   VICecikl( )
   VICfreecut( )
   VICfreelst( )
   VICgapfn( )
   VICkconf( )
   VICkpreduce( )
   VIClci( )
   VICsearchcut( )
   VICsort( )
   VICuflcov( )
   VICuflohi( )
   VICsnffenchel( )
   VICkpsep( )
   VICsetdefaults( )

B. Functions used for internal purposes only ("SCOPE: Intern")
   VICallocut( )
   VICcovsep( )
   VICcovsepheu( )
   VICdcompare( )
   VICfreecut( )
   VICicompare( )
   VICkidefault( )
   VICkplift( )
   VICkpsepaddcol( )
   VICkpsepsub( )
   VICinitsepmaster( )
   VICintcoeff( )
   VICkirestore( )
   VICisintvec( )
   VICliftkc( )
   VICshrinkcut( )
   VICswap( )

********************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "combo.h"
#include "ssfctp.h"
#include <ilcplex/cplex.h>
#include "vico.h"

#define VICzero 1.0E-9
#define VICtol 1.0E-4     /* real number x assumed to be zero if abs(x)<VICtol*/
#define VICone 0.9999     /* real number x <= 1 assumed to be one if x>VICone */
#define VICinf 1.0E+31    /* approximates infinity                            */
#define VICbig 2147483647 /* largest integer                                  */

#define MAX( x, y ) ( (x) > (y) ? (x) : (y) )
#define MIN( x, y ) ( (x) < (y) ? (x) : (y) )
#define GETEDGE( i, j ) ( (i) > (j) ? ((i*(i+1))/2+j) : ((j*(j+1))/2+i) )

int*    VICiptr   = NULL; /* pointer to an array of integers; used for sorting */
double* VICdptr   = NULL; /* pointer to an array of doubles; used for sorting  */
int     VICsrtsgn = 1;    /* sorting direction: 1->ascending, -1->descending   */

/* structure to store a depot subset. Used in function VICcflsmi() */
typedef struct TJSet {
  int card;
  char* inset;
  struct TJSet* succ;
} TJSet;

TFENCHELopt FCHELopt;

/*------------------------------------------------------------------------------*/

void VICsetdefaults( void ) {
/* Set parameters to default values */
  FCHELopt.maxit    = 100;
  FCHELopt.sg_strat = SG_SMOOTH;
  FCHELopt.alpha    = 2.0;
  FCHELopt.H        = 5;
  FCHELopt.CHK      = 0;
  FCHELopt.Algo     = 0;
  FCHELopt.Freq     = 1;
}

/*------------------------------------------------------------------------------*/

void VICaddtolst( VICcut** first, VICcut* firstnew ) {
/*
   Description:
   -------------
   Adds a linked list of new cuts to a linked list of already existing cuts.
   The list of new cuts is added add the top of the list of existing cuts

   Scope: Export
   -------------

   Parameters:
   -------------
   first    : pointer to the pointer to the first cut in the existing list of
              cuts. If *first is NULL the list of existing cuts is assumed to
	      be empty
   firstnew : pointer to the first cut in the list of new cuts to be inserted
              in the list of existing cuts
*/
  VICcut *tmp=*first;
  VICcut *last=firstnew;

  if ( firstnew != NULL ) {
    *first = firstnew;
    while ( last->nextcut != NULL ) last = last->nextcut;
    last->nextcut=tmp;
  }

}

/*------------------------------------------------------------------------------*/

void VICfreecut( VICcut** cut ) {
/*
  Desciption: free memory allocated for the cut *cut and return NULL pointer
  Scope     : Export
  Parameters: pointer to the pointer to the cut
*/
  if ( (*cut)->nzval ) free ( (*cut)->nzval );
  if ( (*cut)->nzind ) free ( (*cut)->nzind );
  if ( (*cut)->UsrDatPtr ) free ( (*cut)->UsrDatPtr );
  free( (*cut) );
  *cut = NULL;
}

/*------------------------------------------------------------------------------*/

void VICfreelst ( VICcut** first ) {
/*
   Description:
   -------------
   Free the memory allocated by a linked list of cuts

   Scope: Export
   -------------

   Parameters:
   -------------
   - first : pointer to the pointer to the first cut in the linked list of cuts
*/
   VICcut* delcut;
   while ( *first != NULL ) {
     delcut = *first;
     *first = (*first)->nextcut;
     VICfreecut( &delcut );
   }
}

/*------------------------------------------------------------------------------*/

char VICallocut( int numnz, VICcut** cut ) {
/*
   Description: Allocate memory required by a cut
   Scope      : Export
   Parameters : numnz is number of nonzeros in the cut.
                If not successfull the NULL pointer is returned in *cut
*/

   VICcut *newcut= NULL;

   *cut = NULL;
   if ( numnz <= 0 ) return( 0 );
   newcut = (VICcut*) malloc( sizeof(VICcut) );
   if ( newcut == NULL ) return( 0 );
   newcut->sense = 'L';
   newcut->rhs = 0.0;
   newcut->nzcnt = numnz;
   newcut->nzind = NULL;
   newcut->nzval = NULL;
   newcut->UsrRVal = 0.0;
   newcut->UsrIVal = 0;
   newcut->UsrDatPtr = NULL;
   newcut->nextcut = NULL;
   newcut->nzval = (double*) calloc( numnz, sizeof(double) );
   newcut->nzind = (int*) calloc( numnz, sizeof(int) );
   if ( (newcut->nzval==NULL) || (newcut->nzind==NULL) ) VICfreecut( &newcut );
   *cut = newcut;
   return( ( cut != NULL ) );
}

/*------------------------------------------------------------------------------*/

VICcut* VICsearchcut ( VICcut* cutlst, VICcut* cut ) {
/*
   Description:
   -------------
   Searches for a cut in given list of cuts

   Scope: Export
   -------------

   Parameters:
   -------------
   - cutlst : pointer to the first cut in the list which has to be searched
   - cut    : pointer to the cut which is search in the list of cuts

   Return value:
   ------------
   NULL if the cut *cut is not containted in the list, and
   the pointer to the cut in this list otherwise
*/
  VICcut *cur_cut;
  char   found;
  int    nz;

  if ( cut==NULL ) return( NULL );
  for ( cur_cut=cutlst; cur_cut != NULL; cur_cut = cur_cut->nextcut ) if (cut != cur_cut) {
    found = ( (cut->nzcnt) == (cur_cut->nzcnt) );
    if ( found ) found = ( fabs(cut->rhs - cur_cut->rhs) < VICtol );
    for ( nz=0; (found) && (nz < cut->nzcnt); nz++ )
      found = ( cut->nzind[nz] == cur_cut->nzind[nz] );
    for ( nz=0; (found) && (nz < cut->nzcnt); nz++ )
      found = ( fabs( cut->nzval[nz] - cur_cut->nzval[nz] ) < VICtol );
    if (found) return( cur_cut );
  }
  return( NULL );
}

/*------------------------------------------------------------------------------*/

void VICclearlst ( VICcut** first ) {
/*
  Description:
  ------------
  Removes duplicated cuts in a list of cuts

  Scope: Export
  -------------

  Parameters:
  - first : pointer to the pointer to the first cut in the list of cuts
*/
  VICcut *cut = *first, *foundcut, *predcut=NULL, *delcut;

  while ( cut != NULL ) {
    foundcut = VICsearchcut( cut->nextcut, cut );
    if ( foundcut ) { /* cut is duplicate: remove the cut */
      delcut = cut;
      cut = cut->nextcut;
      VICfreecut( &delcut );
      if ( predcut == NULL ) *first = cut;
      else predcut->nextcut=cut;
    } else {
      predcut = cut;
      cut = cut->nextcut;
    }
  }

}

/*------------------------------------------------------------------------------*/

static
char VICshrinkcut( int numnz, VICcut* cut ) {
/*
   Description: reduce the space allocated for the cut *cut to space for
                numnz nonzeros
   Scope      : Intern
   Parameters : numnz is new number of nonzeros, "cut" to the cut
                return 1 if successful and 0 otherwise
*/

  int* indices=NULL;
  double* values= (double*) realloc( cut->nzval, sizeof(double)*numnz );

  if ( values == NULL ) {
    VICfreecut(&cut);
    return( 0 );
  }
  cut->nzval = values;
  indices = (int*) realloc( cut->nzind, sizeof(int)*numnz );
  if ( indices==NULL ) {
    VICfreecut(&cut);
    return(0);
  }
  cut->nzind = indices;
  cut->nzcnt = numnz;
  return( 1 );
}

/*------------------------------------------------------------------------------*/

static
void VICswap( int* ptr1, int* ptr2 ) {
/*
  Description: swap contents of addresses to which ptr1 and ptr2 point
  Scope      : Intern
  Parameters : ptr1, ptr2 pointer to integers
*/
  int tmp;
  tmp=*ptr1, *ptr1=*ptr2, *ptr2=tmp;
}

/*------------------------------------------------------------------------------*/

#ifdef SUNPAS
short int VICdcompare ( int *Order1, int* Order2 )
#else
int VICdcompare ( const int *Order1, const int *Order2 )
#endif
{
/*
  Description : comparison function used by qsort (stdlib.h) for sorting
  --------------

  Scope: Intern
  --------------

  Parameters:
  --------------
  - Order1, Order2 : The two elements to be sorted by access to VICdptr
*/
  if ( VICdptr[*Order1] < VICdptr[*Order2] )
    return( -VICsrtsgn );
  else if ( VICdptr[*Order1] > VICdptr[*Order2] )
    return( VICsrtsgn );
  else return( 0 );

}

/*------------------------------------------------------------------------------*/

#ifdef SUNPAS
short int VICicompare ( int *Order1, int* Order2 )
#else
int VICicompare ( const int *Order1, const int *Order2 )
#endif
{
/*
  Description : comparison function used by qsort (stdlib.h) for sorting
  --------------

  Scope: Intern
  --------------

  Parameters:
  --------------
  - Order1, Order2 : The two elements to be sorted by access to VICdptr
*/
  if ( VICiptr[*Order1] < VICiptr[*Order2] )
    return( -VICsrtsgn );
  else if ( VICiptr[*Order1] > VICiptr[*Order2] )
    return( VICsrtsgn );
  else return( 0 );

}

/*------------------------------------------------------------------------------*/

void VICsort( int n, int ascending, int doinit, int size, void* numbers, int* order ) {
/* Description:
   ------------
   sorts an array of integers

   Parameters:
   ------------
   - n         : number of integer numbers to be sorted
   - ascending : array is sorted in ascending order if ascending=1
   - doinit    : doinit=1 indicates that the array order is not initialised
   - numbers   : pointer to an array of integers of dimension of at least n
                 containing the numbers that have to be sorted
   - order     : pointer to an array of integers of dimension of at least n.
                 On output order determines the sorted array as
		 numbers[order[0], ..., order[n-1]]
*/
  int i;
  if (doinit) for ( i=0; i < n; i++ ) order[i]=i;
  if (ascending) VICsrtsgn = 1; else VICsrtsgn = -1;
  if (size <= sizeof(int)) {
    VICiptr = (int*)numbers;
    qsort( order, n, sizeof(int), (void*) VICicompare );
  } else {
    VICdptr = (double*)numbers;
    qsort( order, n, sizeof(int), (void*) VICdcompare );
  }
}

/*------------------------------------------------------------------------------*/

static
void VICcovsepheu ( char useall, int n, int cap, int* weight, double* xlp,
                    double* redcost, int* cover, int* card ) {
/* Decription:
   -------------
   Determines a minimal cover from a knapsack inequality

                   \sum_{j=1}^n weight[j]* x[j]<= cap

   by means of Dantzig's greedy method. The separation problem is to solve

                     \min \sum_{j=1}^n ( 1 - xlp[j] )*z(j)
    s.t.:            s.t.:\sum_{j=1}^n weight[j]*z(j) >= cap+1
                          z(j) \in {0,1} for j=1,...,n

  Scope: Intern
  -------------

  Parameters:
  -------------
  - useall  : if equal to 1 even variables x[j] with LP value xlp[j] of 0 may
              be put into the cover; otherwise only variables with LP value
	      larger than zero are included in the cover
  - n       : number of variables in the knapsack inequality
  - cap     : capacity of the knapsack
  - weight  : pointer to an array of integers of size of at least n containing
              the weights of variables/items in the knapsack inequality
  - xlp     : pointer to an array of doubles of size of at least n containing
              the LP solution
  - redcost : pointer to an array of doubles of size of at least n containing the
              reduced costs of the variables.
	      Reduced cost values are used to decide which variables to
	      remove from a cover if the initial cover is not minimal.
	      redcost may be NULL. Reduced cost values of nonbasic variables
	      at their lower bound are nonnegative, and reduced cost values
	      of nonbasis variables at their upper bound are nonpositive.
  - cover   : pointer to an array of integers of dimension of at least n.
              On output, the first *card entries form the cover and the
	      last n-*card entries are items/variables not in the cover
  - card    : pointer to an integer containing the cardinality of the cover.
              Card is 0 if no cover is found, which is the case if all (free)
	      variables/items j fit into the knapsack (useall=1) or if all
	      (free) items with LP value larger than zero fit into the
	      knapsack (useall=0)
  */

  int     crd, idx, j, jrem, last, notmin, npos, wsum;
  double  Hirco, Lowrco, Lowx, rco;
  double  *ratio = (double*) calloc(n, sizeof(double));

  *card=0;

  /* Determine ratio of "cost coefficients" to weights */
  for ( j=0, npos=0, last=n; j < n; j++ ) {
    ratio[j] = ( 1.0-xlp[j] )/weight[j];
    if ( (useall) || (xlp[j] > VICtol) )
      cover[npos++] = j;
    else last--, cover[last]=j;
  }

  /* sort items in non-descending order of ratios */
  VICdptr   = ratio;
  VICsrtsgn = 1;
  qsort( cover, npos, sizeof(int), (void*) VICdcompare );
  free(ratio);

  /* Use above ordering to put items in the cover as long there is
     enough capacity */
  *card=n, crd=0, wsum=0;
  while ( (crd < npos) && (wsum <= cap) ) wsum += weight[cover[crd++]];
  if ( wsum <= cap ) return;

  /* Transform cover into a minimal one. In case that an item can be
     removed from the cover, remove the item j with minimal LP value
     xlp[j]; in case of a ty and xlp[j] < 1 prefer the item with a
     high magnitude of reduced costs; in case of a ty and xlp[j]=1
     prefer the item with small magnitude of reduced cost */
  notmin = 1;
  Hirco  = -1.0;
  Lowrco = VICinf;
  Lowx   = VICinf;
  while ( notmin ) {
    jrem = -1;
    j    = crd-1;
    while ( j >= 0 ) {
      idx = cover[j];
      if ( wsum - weight[idx] > cap ) {
        if ( xlp[idx] < Lowx )
	  jrem = j;
        else if ( ( redcost ) && ( ! ( xlp[idx] > Lowx ) ) ) {
	  rco = fabs(redcost[idx]);
	  if ( xlp[idx] > VICone ) {
	    if ( rco < Lowrco ) jrem = j;
	  }
	  else if ( rco > Hirco ) jrem = j;
	}
      }
      if ( j==jrem ) {
        idx = cover[jrem];
        Lowx = xlp[idx];
	if ( redcost ) {
	  if ( Lowx < VICone ) Hirco = fabs(redcost[idx]); else Lowrco = fabs(redcost[idx]);
	}
      }
      j--;
    }
    notmin = ( jrem >= 0 );
    if ( notmin ) { /* remove item cover[jrem] from the cover */
      crd--;
      wsum -= weight[cover[jrem]];
      VICswap( cover+crd, cover+jrem );
    }
  }
  *card=crd;
}

/*------------------------------------------------------------------------------*/

static
char VICcovsep ( int n, int cap, int* weight, double* xlp, int* cover, int* card ) {
/* Description:
   -------------
   Exact separation of a violated cover inequality. See function VICcovsepheu
   for a description of the separation problem to be solved.

  Scope: Intern
  -------------

  Parameters:
  -------------
  - n       : number of variables in the knapsack inequality
  - cap     : capacity of the knapsack
  - weight  : pointer to an array of integers of size of at least n
              containing the weights of variables/items in the knapsack
	      inequality
  - xlp     : pointer to an array of doubles of size of at least n
              containing the LP solution
  - cover   : pointer to an array of integers of size of at least n.
              On output, the first *card entries form the cover and the
	      last n-*card entries are items/variables not in the cover.
  - card    : pointer to an integer containing the cardinality of the cover.
              Card is 0 if no cover is found, which is the case if all
	      variables/items j fit into the knapsack (useall=1) or if all
	      items with LP value larger than zero fit into the knapsack
	      (useall=0)
*/
  item   *items = ( item* ) calloc( n, sizeof( item ) );
  int    idx, j, IObj, icoeff, kpcap, nn, last, ones, ub, wone, wsum;
  char   Viol;
  double SPObj;

  /* Solve the separation problem, which is a binary knapsack problem,
     by means of the combo algorithm of Martello, Toth, Pisinger */
  SPObj = (double) n;
  kpcap = -cap-1;
  for ( j=0,nn=0,ones=0,wone=0,wsum=0,ub=0; j < n; j++ ) {
    SPObj -= xlp[j];
    kpcap += weight[j];
    icoeff = (1.0-xlp[j])*10000.0;
    if ( icoeff > 0 ) {
      wsum += weight[j];
      ub += icoeff;
      items[nn].p = icoeff;
      items[nn].w = weight[j];
      items[nn].x = 1;
      items[nn++].num = j;
    }
    else wone += weight[j], cover[ones++] = j;
  }
  *card = ones;
  if ( wsum > kpcap ) {
    if ( nn > 1 )
      IObj = combo( items, items+nn-1, kpcap, 0, ub, 1, 0 );
    else items->x=0;
  }
  for ( j=0,wsum=wone,last=n; j < nn; j++ ) {
    idx=items[j].num;
    if ( items[j].x ) {
      SPObj -= (1.0-xlp[idx]);
      cover[--last] = idx;
    }
    else {
      cover[(*card)++] = idx;
      wsum += weight[idx];
    }
  }
  Viol = ( SPObj < VICone );

  /* remove variables with LP value of 1 from the cover until
     the cover is minimal */
  for ( j=ones-1; j >= 0; j-- ) if ( wsum - weight[cover[j]] > cap ) {
    (*card)--;
    wsum -= weight[cover[j]];
    VICswap( cover+j, cover+(*card) );
  }

  free( items );
  return( Viol );
}

/*------------------------------------------------------------------------------*/

static
void VICkidefault( int n, char sense, char* invers, int* cap, int* weight, double* xlp ) {
/*
  Description:
  -------------
  Transform an inequality encompassing binary variables in a default knapsack inequality

  Scope: Intern
  -------------

  Parameters:
  -------------
  - n      : number of variables in inequality
  - sense  : sense (that is <= or >=) of original inequality
  - invers : on output invers[j]=1 (0) if variable j was complemented or not
  - cap    : right hand side of original and transformed inequality
  - weight : coefficients in original and transformed inequality
  - xlp    : fractional solution values of original and transformed variables
*/
  int j;

  if ( sense == 'G' ) {
    *cap *= -1;
    for ( j=0; j < n; j++ ) weight[j] *= -1;
  }
  for ( j=0; j < n; j++ ) {
    invers[j] = ( weight[j] < 0 );
    if ( invers[j] ) {
      *cap -= weight[j];
      weight[j] *= -1;
      if (xlp) xlp[j] = 1.0 - xlp[j];
    }
  }

}

/*------------------------------------------------------------------------------*/

static
void  VICkirestore( int n, char sense, char* invers, int* cap, int* weight,
                    double* xlp) {
/*
  Description:
  --------------
  Restore original data and express generated cuts in original variables in
  case that a knapsack inequality was transformed to default form

  Scope: Intern
  --------------

  Parameters:
  --------------
  - n      : number of variables in the cut
  - sense  : sense (i.e. <= or >=) of original inequality
  - invers : invers[j]=1 if binary variable was complemented (i.e.
             substituted by 1-x)
  - weight : coefficients of variables in transformed knapsack inequality
  - xlp    : fractional solution of (transformed) variables
*/
  int     j, rhs=*cap;

  for ( j=0; j < n; j++ ) if ( invers[j] ) {
    rhs -= weight[j];
    weight[j] *= -1;
    if (xlp) xlp[j] = 1.0-xlp[j];
  }
  *cap = rhs;

  if ( sense == 'G' ) {
    for ( j=0; j < n; j++ ) weight[j] *= -1;
    *cap = -rhs;
  }

}

/*------------------------------------------------------------------------------*/

void VICkpreduce ( int WHATRED, int TRYCLI, int n, int* cap, char sense,
                   int* weight, int* order, int* indx, VICcut** clique ) {
/*
  Description:
  ------------
  Applies coefficient reduction to a knapsack problem/inequality given by
       w_1 x_1 + ... + w_n x_n <= c, where x_j = 0,1 for j = 1,...,n
  The strengthening of the inequality is based on the following steps:
  1) if x_j = 1 is feasible and implies x_k=0 for all k \ne j then set w_j = c
  2) Try to find a clique by first finding a cover of 2 items and obtaining
     the cover's extension
  3) Reduce coefficients of variables in the clique by replacing the original
     inequality by
       \sum_{j\in C} (w_j - d)x_j + \sum_{j\notin C} w_j x_j <= c-d
     where
       d = max   \sum_{j \notin C} w_j x_j
           s.t.: \sum_{j \notin C} w_j x_j <= c
	          x_j = 0,1 for all j \notin C

  Parameters:
  -----------
   - WHATRED:  Determines which coeffiecient improvement scheme is applied:
               WHAT_RED = 0 -> no coefficient improvement at all
	       WHAT_RED = 1 -> just simple coefficient improvement
	                      (increase w_j to c if x_j=1 implies all other
		               x_k = 0)
               WHAT_RED = 2 -> only apply reduction of coefficients of
	                       variables in the clique
               WHAT_RED = 3 -> do both types of improvement if possible
   - TRYCLI  : If TRYCLI=0 no cliques are generated, otherwise cliques
               are derived if possible
   - n       : number of variables in the knapsack inequality
   - cap     : right hand side of the knapsack inequality
   - sense   : sense (direction) of the knapsack inequality.
               In case of sense='L' the inequality is given by
               weight[0]*x[indx[0]] +....+ weight[n-1]*x[indx[n-1]] <= cap
	       Otherwise it is assumed that the inequality is
               weight[0]*x[indx[0]] +....+ weight[n-1]*x[indx[n-1]] >= cap
   - weight  : pointer to an array of integers of size of at least n
               containing the coefficients of variables/items in the knapsack
	       inequality. The coefficients are not restricted to be
	       nonnegative!
   - order   : pointer to an array of integers of size of at least n
               determining an ordering according to increasing weights
	       of the items in the knapsack. order may be NULL.
   - indx    : pointer to an array of integers of size of at least n
               containing the indices of the variables arising in the
	       knapsack inequality.
	       If indx is NULL, it is assumed that these variables are
	       numbered from 0 to n-1
   - clique  : pointer to a pointer to a cut structure VICcut defined in
               vico.h. *clique points to the generated clique inequality
	       if one is found. If no clique was found, then *clique=NULL
*/
  char   *invers = calloc(n, sizeof(char) );
  char   dosort = (order==NULL), dorestore=0;
  int    delta, k, j, jj, nn, numnz, rcap, wmin, wsum;
  item   *items = NULL;

  *clique = NULL;
  if (dosort) order = (int*) calloc( n, sizeof(int) );
  if ( (invers==NULL)||(order==NULL) ) goto RETURN;

  /* transform knapsack inequality to the default form if required */
  dorestore = 1;
  VICkidefault( n, sense, invers, cap, weight, NULL );

  /* sort items according to increasing weights */
  if ( dosort ) {
    VICiptr   = weight;
    VICsrtsgn = 1;
    for ( j=0; j < n; j++ ) order[j] = j;
    qsort( order, n, sizeof(int), (void*) VICicompare );
  }

  /* If x_j=1 implies x_k=0 forall k\ne j increase coeff of x_j to rhs */
  wmin = weight[order[0]];
  rcap = *cap - wmin;
  if ( (WHATRED==1) || (WHATRED==3) ) for ( j=1; j < n; j++ ) {
    k = order[j];
    if ( (weight[k] < *cap ) && ( weight[k] > rcap ) ) weight[k] = *cap;
  }

  /* Try to find a cover of 2 items */
  if ( !(TRYCLI) ) goto RETURN;
  for ( k=1, wsum=0; k < n; k++ ) {
    wsum = weight[order[k-1]]+weight[order[k]];
    if (wsum > *cap) break;
  }
  if (wsum <= *cap) goto RETURN;

  /* The clique consists of variables order[k-1], order[k], ...,order[n-1] */
  numnz = n-k+1;
  if ( ! (VICallocut( numnz, clique )) ) goto RETURN;
  (*clique)->sense = 'L';
  (*clique)->rhs = 1.0;
  (*clique)->nzcnt = numnz;
  for ( j=k-1,numnz=0; j < n; j++ ) {
    jj = order[j];
    (*clique)->nzval[numnz] = (invers[jj]) ? -1.0 : 1.0;
    (*clique)->nzind[numnz++] = (indx) ? indx[jj] : jj;
    if ( invers[jj] ) (*clique)->rhs -= 1.0;
  }

  /* Now try do reduce the coefficients of clique members */
  if ( WHATRED < 2 ) goto RETURN;
  nn    = k-1;
  delta = *cap;
  if ( nn > 1 ) {
    items = (item*) calloc( nn, sizeof(item ) );
    if ( items ) {
      for ( j=0,wsum=0; j < nn; j++ ) {
        wmin = weight[order[j]];
        wsum += wmin;
        items[j].p = wmin;
        items[j].w = wmin;
      }
      delta = wsum;
      if (wsum > *cap) delta = combo( items, items+nn-1, *cap, 0, *cap, 0, 0 );
    }
  } else if ( nn==1 ) delta = weight[order[0]];
  /* Reduce right hand side and clique coeffs by *cap - delta */
  delta = *cap - delta;
  if ( delta > 0 ) {
    (*cap) -= delta;
    for (j=k-1; j < n; j++ ) weight[order[j]] -= delta;
  }

RETURN:
  if ( dorestore ) VICkirestore( n, sense, invers, cap, weight, NULL );
  if (invers) free(invers);
  if (items) free(items);
  if ( (dosort) && (order) ) free(order);
}
/*------------------------------------------------------------------------------*/

void VIClci ( int n, int cap, char sense, int* weight, int* indx, double* xlp,
              double* rco, VICcut** lci ) {
/* Description:
   -------------
   Heuristic of Savelsbergh et al. for determining a lifted cover inequality
   which cuts off a given fractional solution. For more information on this
   procedure see:
   Gu Z, Nemhauser GL, Savelsbergh MWP (1998). Lifted cover inequalities
   for 0-1 linear programs: Computation. INFORMS J. on Computing 10:427-437

   Scope: Export
   -------------

   Parameters:
   ------------
   - n       : number of variables in the knapsack inequality
   - cap     : right hand side of the knapsack inequality
   - sense   : sense (direction) of the knapsack inequality. In case of
               sense='L' the inequality is given by
               weight[0]*x[indx[0]] + .... + weight[n-1]*x[indx[n-1]] <= cap
	       Otherwise it is assumed that the inequality is
               weight[0]*x[indx[0]] + .... + weight[n-1]*x[indx[n-1]] >= cap
   - weight  : pointer to an array of integers of size of at least n
               containing the coefficients of variables/items in the knapsack
	       inequality. The coefficients are not restricted to be
	       nonnegative!
   - indx    : pointer to an array of integers of size of at least n
               containing the indices of the variables arising in the
	       knapsack inequality.
	       If indx is NULL, it is assumed that these variables are
	       numbered from 0 to n-1
   - xlp     : pointer to an array of doubles of size of at least n
               containing the LP solution.
   - rco     : pointer to an array of doubles of size of at least n
               containing the absolute values of reduced costs
   - lci     : pointer to a pointer to a cut structure VICcut defined in
               vico.h. If no violated lifted cover is found, the null
	       pointer is returned; otherwise **lci contains the cut.
*/

  double lhs;
  double *crb = (double*) calloc(n, sizeof(double) );
  char  *down = (char*) calloc(n, sizeof(char) );
  char  *invers = (char*) calloc(n, sizeof(char) );
  int   *cover = (int*) calloc(n, sizeof(int) );
  int   *coeff = (int*) calloc(n, sizeof(int) );
  int   i, j, jj, k, last, ncov, npos, numnz, ones, r, t, wsum;
  int   *partialsum=NULL, psize=1000, *lseq=NULL, rcap, rhs;

  *lci = NULL;
  /* Transform Knapsack inequality to default */
  VICkidefault( n, sense, invers, &cap, weight, xlp );

  /* Sort variables with positive LP-value according to nonincreasing LP values */
  for ( j=0,npos=0,last=n; j < n; j++ ) {
    down[j] = 0;
    if ( weight[j] > cap ) rco[j] = CPX_INFBOUND;
    if ( xlp[j] > VICtol ) cover[npos++]=j; else cover[--last]=j;
  }
  VICdptr   = xlp;
  VICsrtsgn = -1;
  qsort( cover, npos, sizeof(int), (void*) VICdcompare );

  /* put variables into cover until capacity is exceeded */
  ncov=0, wsum=0;
  while ( (ncov < npos) && (wsum <= cap) ) wsum += weight[cover[ncov++]];
  if ( wsum <= cap ) goto RETURN;

  /* variables in cover with LP value of 1 will be down-lifted */
  ones = 0;
  while ( (ones<n)&&(xlp[cover[ones]]>VICone) ) down[cover[ones++]]=1;

  /* remove items from cover until initial cover is minimal */
  j=ncov-1;
  while ( j >= ones) {
    if ( wsum - weight[cover[j]] > cap ) {
      ncov--;
      wsum -= weight[cover[j]];
      VICswap(cover+j,cover+ncov);
    } else j--;
  }

  /* Determine the lifting sequence as follows:
     (1) variables with positive LP value not in the cover:
         lifting sequence according to decreasing value
	 of LP values times weight
     (2) variable in the cover with LP value of 1:
         lifting sequence according to ascending absolute
	 values of reduced costs
     (3) variables not in the cover with LP value of 0:
         lifting sequence according to ascending absolute
	 values of reduced costs
  */
  lseq = (int*) calloc(n, sizeof(int) );
  /* 1. variables in the cover with LP value less than 1 */
  for ( j=ones; j < ncov; j++ ) lseq[j-ones] = cover[j];
  /* 2. variables not in the cover with positive LP value */
  for ( j=ncov; j < npos; j++ ) {
    k=cover[j];
    crb[k] = weight[k]*xlp[k];
    lseq[j-ones] = k;
  }
  VICdptr = crb;
  VICsrtsgn = -1;
  last = ncov-ones;
  qsort( lseq+last, npos-ncov, sizeof(int), (void*) VICdcompare );
  /* 3. variables in the cover with LP value of one */
  for ( j=0; j < ones; j++ ) lseq[npos-ones+j] = cover[j];
  VICdptr = rco;
  VICsrtsgn=1;
  qsort( lseq+npos-ones, ones, sizeof(int), (void*) VICdcompare );
  /* 4. variables not in the cover with LP value of zero */
  for ( j=npos; j < n; j++ ) lseq[j]=cover[j];
  qsort( lseq+npos, n-npos, sizeof(int), (void*) VICdcompare );

  free(cover);
  cover = lseq;

  /* Lift the inequality by means of dynamic programming : */

  if ( last > psize-1 ) psize = last+1;
  partialsum = (int*) calloc(psize,sizeof(int) );

  /* a. sort items in cover with LP value less than 1 according
        to increasing weights */
  VICiptr = weight;
  VICsrtsgn = 1;
  qsort( cover, last, sizeof(int), (void*) VICicompare );
  r   = last;
  rhs = r-1;
  partialsum[0]=0;
  for ( j=0,lhs=0; j < last; j++ ) {
    k = cover[j], lhs += xlp[k], coeff[k]=1, partialsum[j+1] = weight[k]+partialsum[j];
  }
  for ( j=last,rcap=cap; j < n; j++ ) if ( down[cover[j]] ) rcap -= weight[cover[j]];

  /* b. Apply the dynamic programming algorithm */
  for ( i=r,numnz=r; i < n; i++ ) {
    k=cover[i], t = r;
    if ( weight[k] > cap ) break; /* variables cover[k]...cover[n-1] must be zero */
    if ( down[k] ) rcap += weight[k]; else rcap -= weight[k];
    while ( ( partialsum[t] > rcap ) && ( t > 0 ) ) t--;
    if ( down[k] ) {
      coeff[k] = t - rhs;
      if ( coeff[k] > 0 ) rhs += coeff[k]; else coeff[k] = 0;
    }
    else {
      coeff[k] = MAX(0, rhs - t);
      rcap += weight[k];
    }
    if ( coeff[k] > 0 ) {
      numnz++;
      lhs += coeff[k]*xlp[k];
      if ( r+coeff[k]+1 > psize ) {
        psize = MAX( psize+psize/2, r+coeff[k]+1 );
        partialsum = realloc( partialsum, psize*sizeof(int) );
      }
      for ( t=r+coeff[k]; t > r; t-- ) partialsum[t]=partialsum[t-coeff[k]]+weight[k];
      for ( t=r; t >= coeff[k]; t-- )
        partialsum[t] = MIN( partialsum[t], weight[k]+partialsum[t-coeff[k]] );
      for ( t=coeff[k]-1; t >= 0; t-- ) partialsum[t]=MIN(partialsum[t],coeff[k]);
      r += coeff[k];
    }
    if ( ( xlp[k] < VICtol ) && (lhs < rhs) ) break;
  }

  /* If inequality is violated, return the cut in structure lci */
  if ( lhs*VICone <= rhs ) goto RETURN;
  if ( ! (VICallocut( numnz, lci )) ) goto RETURN;
  (*lci)->sense = 'L';
  (*lci)->rhs = (double)rhs;
  (*lci)->nzcnt = numnz;
  for ( j=0,numnz=0; j < n; j++ ) {
    jj = cover[j];
    if ( coeff[jj] > 0 ) {
      if ( invers[jj] ) {
        coeff[jj] *= -1;
	(*lci)->rhs += (double)coeff[jj];
      }
      (*lci)->nzval[numnz] = (double)coeff[jj];
      (*lci)->nzind[numnz++] = (indx) ? indx[jj] : jj;
    }
  }


RETURN:
  /* transform data and inequality to original form */
  VICkirestore( n, sense, invers, &cap, weight, xlp );
  free( down );
  free( invers );
  free( cover );
  free( coeff );
  free( crb );
  if ( partialsum ) free( partialsum );

}

/*------------------------------------------------------------------------------*/

static
int VICliftkc( int n, int rcap, int* W, double* xlp, int* order, int* coeff ){
/*
  Description:
  --------------
  performs uplifting of a (1,k)-configuration inequality given by

       x[order[0]] +... + x[order[r-1]] + coeff[order[r]]*x[order[r]] <= r

  where r=coeff[n]. The procedure is based on Zemel's dynamic programming approach.

  Scope: Intern
  --------------

  Parameters:
  ------------
  - n    : number of variables in knapsack constraint
  - rcap : capacity of knapsack
  - W    : W[j] is the weight of the variable j in the knapsack constraint
  - xlp  : xlp[j] = LP value of variable j
  - order: variables order[0], ..., order[r] have nonzero coefficient on
           input. Furthermore W[order[0]] <= ...<= W[order[r-1]].
           order[r+1], ..., order[n-1] gives the lifting sequence, where
           variables j with positive LP value xlp[j] have to appear first.
  - coeff: On input coeff[order[j]] = 1 for j=0, ..., r-1 and
            0 < coeff[order[r]] <= r
           On output coeff[j] is the lifted coefficient of a variable
           which had a zero coefficient on input.

  Return value: 1 if lifted inequality is violated by xlp and 0 otherwise
  -------------
*/
int    i, k, t, T, r;
int    wsize, nwsize, isviol, tmp;
double lhs;
int*   wsums = NULL;

  /* Check if a lifted inequality can be violated */
  r   = coeff[n];
  lhs = coeff[order[r]]*xlp[order[r]];
  for ( i=0; i < r; i++ ) lhs += xlp[order[i]];
  isviol = (lhs > r);
  for ( i=r+1; (i < n) && (isviol==0); i++ ) isviol = ( xlp[order[i]] > VICtol );
  if ( isviol==0 ) return( isviol );

  /* Allocate memory for partial sums of weights */
  wsize = r+coeff[order[r]] < 1000 ? 1000 : r+coeff[order[r]]+1;
  wsums = (int*)calloc(wsize,sizeof(int));
  if ( wsums==NULL ) return( 0 );
  for ( i=0, wsums[0]=0; i < r; i++ ) wsums[i+1] = wsums[i]+W[order[i]];
  k = order[r];
  for ( t = r+coeff[k]; t > r; t-- ) wsums[t] = W[k] + wsums[t-coeff[k]];
  for ( t = r; t >= coeff[k]; t-- ) {
    tmp = W[k]+wsums[t-coeff[k]];
    if ( tmp < wsums[t] ) wsums[t] = tmp;
  }
  for ( t = coeff[k]-1; t > 0; t-- ) if ( W[k] < wsums[t] ) wsums[t] = W[k];

  /* Apply Zemel's dynamic programming approach to uplift a variable's coefficient */
  for ( i=r+1,T=r; (i < n)&&(isviol); i++ )  {
    T += coeff[k];
    k  = order[i];
    t  = T;
    while ( (wsums[t] > rcap - W[k] ) && ( t > 0 ) ) t--;;
    coeff[k] = coeff[n]-t;
    if ( coeff[k] > 0 ) {
      nwsize = T + coeff[k] + 1;
      if ( nwsize > wsize ) {
        wsize = nwsize + nwsize/2;
        wsums = (int*) realloc( wsums, wsize*sizeof(int) );
        if ( wsums==NULL ) return( 0 );
      }
      lhs += coeff[k]*xlp[k];
      for ( t=T+coeff[k]; t > T; t-- ) wsums[t] = W[k] + wsums[t-coeff[k]];
      for ( t=T; t >= coeff[k]; t-- ) {
        tmp = W[k]+wsums[t-coeff[k]];
        if ( tmp < wsums[t] ) wsums[t] = tmp;
      }
      for ( t=coeff[k]-1; t > 0; t-- ) if ( W[k] < wsums[t] ) wsums[t] = W[k];
    }
    else coeff[k]=0;
    if ( xlp[k] < VICtol ) isviol = ( lhs > coeff[n] ) ;
  }
  lhs = lhs*(1.0-VICtol);
  isviol = (lhs > coeff[n]);
  free(wsums);
  return ( isviol );
}

/*------------------------------------------------------------------------------*/

void VICkconf( int n, int cap, char sense, int* weight, int* indx, double* xlp,
               double* rco, VICcut** kconf ) {
/*
  Description:
  --------------
  Tries to get a (1,k)-configuration inequality from a cover using  the separation
  heuristic of Crowder, Johnson, Padberg in Oper. Res. 31 (1983). For an alternative
  separation heuristic for (1,k)-configurations see Carlos E. Ferreira (1997). On
  Combinatorial Optimization Problems arising in Computer Systems Design. Phd Thesis,
  Technische Universit\E4t Berlin.

  Given the Knapsack-Polytop
            X = { x : \sum_{j\in N} w[j]*x[j] <= c, x_j = 0,1 }
  a (1,k)-configuration is a set NP \cup {z}, where NP \subset N, such that

  (i)  \sum_{j \in NP} w[j] <= c
  (ii) The set K \cup {z} is a cover with respect to N for all
       subsets K of NP with cardinality k

  The corresponding (1,k)-configuration inequality is given by

          (r - k + 1)x[z] + \sum_{j\in NP} x[j] <= r, where r=|NP|

  Crowder, Johnson, Padberg propose the following separation heuristic:

  1. Let S \subset N be the cover, and z\in S the item with
     maximum weight. Set NP = S-{z} and k=|NP|.
  2. For all j \in N-S with \sum_{l \in N-S} w[l] <= c do:
     a. Check if K \cup \{z} is a cover for any K \subseteq NP with |K|=k
     b. If this is the case build the corresponding (1,k)-configuration
        inequality and lift it.
     c. If the lifted inequality is violated, add the inequality


  Scope: Export
  --------------

  Parameters:
  -------------
   - n       : number of variables in the knapsack inequality
   - cap     : right hand side of the knapsack inequality
   - sense   : sense (direction) of the knapsack inequality. In case of
               sense='L' the inequality is given by
               weight[0]*x[indx[0]] + .... + weight[n-1]*x[indx[n-1]] <= cap
	       Otherwise it is assumed that the inequality is
               weight[0]*x[indx[0]] + .... + weight[n-1]*x[indx[n-1]] >= cap
   - weight  : pointer to an array of integers of size of at least n
               containing the coefficients of variables/items in the knapsack
	       inequality. The coefficients are not restricted to be
	       nonnegative!
   - indx    : pointer to an array of integers of size of at least n
               containing the indices of the variables arising in the
	       knapsack inequality.
	       If indx is NULL, it is assumed that these variables are
	       numbered from 0 to n-1
   - xlp     : pointer to an array of doubles of size of at least n
               containing the LP solution.
   - rco     : pointer to an array of doubles of size of at least n
               containing the absolute values of reduced costs
   - kconf   : pointer to a pointer to a cut structure VICcut defined in
               vico.h. If no violated inequality is found, the null
	       pointer is returned; otherwise **kconf contains the cut.


   Meaning of some important local variables used here:
   - cover    cover[0, ..., card-1] are the items in the cover, that is the
              variable indices indx[cover[0]], ..., indx[cover[card-1]].
              indx[cover[card-1]] is the item in the cover with maximum weight
              indx[cover[0]], ..., indx[cover[num-1]] is the set NP
              indx[cover[card]], ..., indx[cover[n-1]] is the Set N-S
   - wsum     sum of the weights of all items in NP
   - wksum    sum of the weights of the first k items in NP if sorted in ascending order
              of weights
   - wkmax    weight of the k-th ordered item from NP
   - minw     minimum value required for w[j] for NP \cup {j} \cup {z} beeing a
              (1,k)-configuration
   - maxw     maximum value of w[j] such that sum of weights of items in NP \cup {j}
              does not exceed the capacity c
   - nwksum   new vallue of "WKSum" if j is added to the set NP
   - num      number of items in set NP
   - NPXsum   sum of LP values of the variables in set NP
   - jpos     position at which item j \in N-S has to be inserted in set NP
*/

  int    *cover = (int*) calloc( n, sizeof(int) );
  int    *coeff = (int*) calloc( n+1, sizeof(int) );
  int    wsum=0, wksum, wkmax, card=0, minw, maxw,nwksum, num=0;
  int    k, l, j, jj, jpos, z, numnz;
  char   *invers = calloc(n, sizeof(char) );
  double NPXsum=0.0;
  VICcut *curcut;

  *kconf=NULL;
  VICkidefault( n, sense, invers, &cap, weight, xlp );

  /* determine a cover using Danztig's heuristics */
  VICcovsepheu( 0, n, cap, weight, xlp, rco, cover, &card );
  if ( card==n ) goto RETURN;

  card = 0;
  while ( wsum <= cap ) NPXsum += xlp[cover[card]], wsum += weight[cover[card++]];

  /* sort items in the cover in ascending order of weights */
  VICiptr = weight;
  VICsrtsgn = 1;
  qsort( cover, card, sizeof(int), (void*)VICicompare );
  k      = card-1;
  z      = cover[k];
  wsum  -= weight[z];
  wksum  = wsum;
  wkmax  = weight[cover[k-1]];
  NPXsum-= xlp[z];
  num    = k;

  /* sort items not in the cover in descending order of LP values */
  VICdptr = xlp;
  VICsrtsgn = -1;
  qsort( (cover+card), n-card, sizeof(int), (void*)VICdcompare );

  /* sort items not in the cover with LP values of 0 in ascending
     order of magnitude of reduced cost values */
  VICsrtsgn=1;
  VICdptr = rco;
  l = n;
  do l--; while ( (l>=card) && (xlp[cover[l]] < VICtol) );
  l++;
  if ( (l >= card) && (l<n) ) qsort( (cover+l), n-l, sizeof( int), (void*)VICdcompare);

  /* For every item cover[j] not in the cover check if the set
     NP={cover[0],...,cover[k-1]} \cup cover[j] \cup {z} is a (1,k)-configuration */
  maxw = cap - wsum;
  minw = cap - weight[z] - wksum + wkmax;
  for (j=card; (j < n)&&(minw < maxw); j++) {
    jj = cover[j];
    if ( (weight[jj] > minw) && (weight[jj] <= maxw) ){
      nwksum = weight[jj] < wkmax ? nwksum = wksum - wkmax + weight[jj] : wksum;
      jpos   = nwksum + weight[z] > cap ? num : -1;
      if ( jpos >= 0 ) {
        wksum = nwksum;
        wsum += weight[jj];
        NPXsum += xlp[jj];
        /* insert item cover[j] at right position */
        while ( (jpos > 0) && (weight[jj] < weight[cover[jpos-1]]) ) jpos--;
        for ( l=j; l > jpos; l-- ) cover[l] = cover[l-1];
        cover[jpos] = jj;
        wkmax = weight[cover[k-1]];
        maxw  = cap - wsum;
        minw  = cap - weight[z] - wksum + wkmax;
        num++;
        /* set coefficients of the inequality */
        for ( l=0; l < num; l++ ) coeff[cover[l]] = 1;
        for ( l=num; l < n; l++ ) coeff[cover[l]] = 0;
        coeff[n] = num;
        coeff[z] = num-k+1;
        /* lift the inequality and add it to LP if violated */
        if ( VICliftkc( n, cap, weight, xlp, cover, coeff ) ) {
	  for (l=0,numnz=0; l < n; l++ ) if ( coeff[l] > 0 ) numnz++;
	  if ( VICallocut( numnz, &curcut ) ) {
	    curcut->rhs = (double) coeff[n];
	    curcut->nzcnt = numnz;
	    curcut->sense = 'L';
	    for (l=0,numnz=0; l < n; l++ ) if ( coeff[l] > 0 ) {
	      if ( invers[l] ) {
	        coeff[l] *= -1;
		curcut->rhs += (double)coeff[l];
	      }
	      curcut->nzval[numnz] = (double)coeff[l];
	      curcut->nzind[numnz++] = ( indx ) ? indx[l] : l;
	    }
	    VICaddtolst( kconf, curcut );
	  }
        }
      }
    }
  }

RETURN:
  VICkirestore( n, sense, invers, &cap, weight, xlp );
  free( cover );
  free( coeff );
  free( invers );

}

/*------------------------------------------------------------------------------*/

void VICeci( int n, int cap, char sense, int* weight, int* order, int* indx,
             double* xlp, VICcut** eci ) {
/* Description:
   -------------
   Algorithm of Grabel & Minoux for separating a most violated extended
   cover inequality. For more information on this procedure see:
   - V. Gabrel, M. Minoux (2002). A scheme for exact separation of
     extended cover inequalities and application to multidimensional
     knapsack problems. Oper. Res. Lett. 30:252-264

   Scope: Export
   -------------

   Parameters:
   ------------
   - n       : number of variables in the knapsack inequality
   - cap     : right hand side of the knapsack inequality
   - sense   : sense (direction) of the knapsack inequality. In case of
               sense='L' the inequality is given by
               weight[0]*x[indx[0]] + .... + weight[n-1]*x[indx[n-1]] <= cap
	       Otherwise it is assumed that the inequality is
               weight[0]*x[indx[0]] + .... + weight[n-1]*x[indx[n-1]] >= cap
   - weight  : pointer to an array of integers of size of at least n
               containing the coefficients of variables/items in the knapsack
	       inequality. The coefficients are not restricted to be
	       nonnegative!
   - order   : if not NULL order must contain an ordering of the items
               {0,...,n-1} in the knapsack such that
                 weight[order[i]] <= weight[order[i+1]] for i=0,...n-2
   - indx    : pointer to an array of integers of size of at least n
               containing the indices of the variables arising in the
	       knapsack inequality.
	       If indx is NULL, it is assumed that these variables are
	       numbered from 0 to n-1!!
   - xlp     : pointer to an array of doubles of size of at least n
               containing the LP solution.
   - eci     : pointer to a pointer to a cut structure VICcut defined in
               vico.h. If no violated lifted cover is found, the null
	       pointer is returned; otherwise **eci contains the cut.
*/

  char   *invers = calloc(n, sizeof(char) );
  char   *in_ec  = calloc(n, sizeof(char) );
  int    *Wsums  = calloc(n, sizeof(int) );
  item   *items = ( item* ) calloc( n, sizeof( item ) );
  char   dosort, goon;
  double f_y, g_y, lambda, max_obj, max_rho, new_lam, profit;
  long   scale;
  int    ccard, combo_obj, h, i, j, k, kmin, nn, numnz, wsum, rcap;

  *eci = NULL;
  if ( (invers==NULL)||(in_ec==NULL)||(Wsums==NULL)||(items==NULL) ) goto RETURN;
  dosort = ( order==NULL );
  if ( dosort ) order = (int*) calloc( n, sizeof( int ) );
  if ( order == NULL ) goto RETURN;

  /* transform knapsack inequality to the default form if required */
  VICkidefault( n, sense, invers, &cap, weight, xlp );

  /* sort items according to increasing weights */
  if ( dosort ) {
    VICiptr   = weight;
    VICsrtsgn = 1;
    for ( j=0; j < n; j++ ) order[j] = j;
    qsort( order, n, sizeof(int), (void*) VICicompare );
  }

  /* Determine smallest index kmin such that sum of weights of first
     kmin items exceeds knapsack's capacity, etc */
  for ( j=0, wsum=0, kmin=0; j < n; j++ ) {
    wsum += weight[order[j]];
    Wsums[j] = wsum;
    if ( wsum <= cap ) kmin++;
  }
  if ( kmin == n ) goto RETURN;

  /* If xlp[j] > 0 for any item j with weight[j] > cap then return
     the trivial extended cover inequality with right hand side of
     zero */
  ccard = 1;
  max_rho = 2.0;
  for ( k=kmin; k < n; k++ ) {
    j = order[k];
    if ( (weight[j] > cap) && (xlp[j] > VICtol) ) {
      for ( i=j; i < n; i++ ) in_ec[order[i]]=1;
      goto STORECUT;
    }
  }

  /* Apply algorithm of Gabrel & Minoux for finding most violated ECI */
  max_rho = 0.0;
  for ( k=kmin; k < n; k++ ) { /* Solve separation problem SP^k */
    /* Calculate starting solution */
    wsum = cap - weight[order[k]];
    for ( h=0; h < n; h++ ) if ( Wsums[h] > wsum ) break;
    for ( j=0, f_y=0.0; j <= h; j++ ) f_y += xlp[order[j]];
    for ( j=k; j < n; j++ ) f_y += xlp[order[j]];
    lambda = f_y/( (double)(h+1) );
    do { /* solve knapsack problems by means of combo */
      max_obj = lambda*((double)k);
      scale = 2;
      while ( scale <= max_obj ) scale *= 2;
      scale = (long)(0x80000000UL/scale);
      f_y   = g_y = 0.0;
      rcap  = Wsums[k]-cap-1;
      for ( j=0, nn=0, wsum=0.0; j < k; j++ ) {
        i = order[j];
        profit = lambda - xlp[i];
        if ( (profit < VICtol) || (weight[i]>rcap) ) {
          f_y += xlp[i];
	  g_y += 1;
        } else {
	  items[nn].p   = (long)(profit*scale);
  	  items[nn].w   = weight[i];
	  items[nn].num = i;
          items[nn].x   = ( items[nn].w <= rcap );
	  wsum += items[nn++].w;
        }
      }
      if ( (wsum > rcap) && (nn > 1) ) {
        combo_obj = combo( items, items+nn-1, rcap, 0, (int)(scale*max_obj), 1, 0 );
      }
      for ( j=0; j < nn; j++ ) if ( !(items[j].x) ) {
        i = items[j].num, f_y += xlp[i], g_y += 1.0;
      }
      for ( j=k; j < n; j++ ) f_y += xlp[order[j]];
      new_lam = f_y/g_y;
      goon = ( new_lam > lambda + VICtol );
      if ( goon ) lambda = new_lam;
    } while ( goon );
    /* save current best extended cover */
    if ( lambda > max_rho ) {
      max_rho= lambda;
      ccard  = k+1;
      for ( j=0; j < n; j++ ) in_ec[j] = 1;
      for ( j=0; j < nn; j++ ) if ( items[j].x ) {
        i = items[j].num, in_ec[i] = 0, ccard--;
      }
    }
  }

STORECUT:
  /* if violated then store the extended cover inequality in structure eci* */
    if ( max_rho > 1.0+VICtol ) {
    for ( j=0, numnz=0; j < n; j++ ) numnz += (int)in_ec[j];
    if ( VICallocut( numnz, eci ) ) {
      (*eci)->sense = 'L';
      (*eci)->rhs = (double)(ccard-1);
      for ( j=0, numnz=0; j < n; j++ ) if ( in_ec[j] ) {
        (*eci)->nzval[numnz] = (invers[j]) ? -1.0 : 1.0;
        (*eci)->nzind[numnz++] = (indx) ? indx[j] : j;
	if ( invers[j] ) (*eci)->rhs -= 1.0;
      }
      (*eci)->nzcnt = numnz;
    }
  }

  VICkirestore( n, sense, invers, &cap, weight, xlp );
RETURN:
  if ( invers) free( invers );
  if ( in_ec ) free( in_ec );
  if ( Wsums ) free( Wsums );
  if ( items ) free( items );
  if ( ( dosort ) && ( order ) ) free( order );

}

/*------------------------------------------------------------------------------*/

void VICecikl( int n, int cap, int do_exact, char sense, int* weight, int* indx,
               double* xlp, VICcut** eci ) {
/* Procedure for generating extended cover inequalities as suggested in
   K. Kaparis, A.N. Letchford (2010). Separation algorithms for 0-1 knapsack
   polytopes, Math. Prog. 124:69-91
*/

  *eci = NULL;
  char *invers = calloc(n, sizeof(char) );
  double* xpos = (double*) calloc( n, sizeof(double) );
  int* wpos = (int*) calloc( n, sizeof(int) );
  int* order = (int*) calloc( n, sizeof(int) );
  int* cover = (int*) calloc( n, sizeof(int) );
  int  card, ecard, k, kk, itmp, j, jj, nn, wk, wsum;
  char viol;
  double lhs, rtmp;

  if ( (invers==NULL) || (xpos==NULL) || (wpos==NULL) || (order==NULL)
       || (cover==NULL) ) goto RETURN;

  /* Transform knapsack inequality to the default form if required */
  VICkidefault( n, sense, invers, &cap, weight, xlp );

  /* Remove variables showing zero LP value from consideration */
  for ( j=nn=wsum=0; j < n; j++ ) if ( xlp[j] > VICtol ) {
    wsum += weight[j];
    xpos[nn] = xlp[j];
    wpos[nn] = weight[j];
    order[nn++] = j;
  }

  /* Main loop of algorithm (cf. p. 75 in Kaparis and Letchford (2010) */
  while ( wsum > cap ) {
    // Find min. cover of items in set \tilde{N}
    if ( do_exact )
      viol = VICcovsep( nn, cap, wpos, xpos, cover, &card );
    else
     VICcovsepheu( 1, nn, cap, wpos, xpos, NULL, cover, &card );
    // Let S* be the set of items in the cover of largest weight. Find item k*
    // from S* that shows smallest LP value
    for ( j=1, k=0, kk=cover[0], lhs=xpos[kk]; j < card; j++ ) {
      jj = cover[j];
      lhs += xpos[jj];
      if ( wpos[jj] > wpos[kk] ) k=j, kk=jj;
      else if ((wpos[jj] == wpos[kk]) && (xpos[jj] < xpos[kk])) k=j, kk=jj;
    }
    // Check if the extended cover inequality is violated
    for ( j=ecard=card, wk=wpos[kk]; j < nn; j++ ) if ( wpos[cover[j]] >= wk )
      lhs += xpos[cover[j]];
    if ( lhs > (double)(card-1) + VICtol ) {
      // violated cover inequality found
      if ( VICallocut( ecard, eci ) ) {
        (*eci)->sense = 'L';
        (*eci)->rhs = (double)(card-1);
        for ( j=0; j < card; j++ ) {
          jj = order[cover[j]];
          (*eci)->nzval[j] = (invers[jj]) ? -1.0 : 1.0;
          (*eci)->nzind[j] = (indx) ? indx[jj] : jj;
          if ( invers[jj] ) (*eci)->rhs -= 1.0;
        }
        for ( j=card; j < ecard; j++ ) if ( wpos[cover[j]] >= wk ) {
          jj = order[cover[j]];
          (*eci)->nzval[j] = (invers[jj]) ? -1.0 : 1.0;
          (*eci)->nzind[j] = (indx) ? indx[jj] : jj;
          if ( invers[jj] ) (*eci)->rhs -= 1.0;
        }
        (*eci)->nzcnt = ecard;
        goto RETURN;
      }
    }
    else {
      // Remove the k-th item as well as all items showing larger weight than
      // the k-th from the set \tilde{N}
      j = 0;
      kk = order[kk];
      while ( j < nn )
      {
        jj = order[j];
        if ( (weight[jj] > wk) || ( jj==kk ) )
        {
          wsum-= weight[jj];
          order[j] = order[--nn], order[nn] = jj;
        }
        else j++;
      }
      for ( j=0; j < nn; j++ ) xpos[j]=xlp[order[j]], wpos[j]=weight[order[j]];
    }
  }

  VICkirestore( n, sense, invers, &cap, weight, xlp );

RETURN:
  if ( invers ) free( invers );
  if ( xpos ) free( xpos );
  if ( wpos ) free( wpos );
  if ( order ) free( order );
  if ( cover ) free( cover );
}


/*------------------------------------------------------------------------------*/


void VICcflfc( int m, int n, int* demand, int* capaci, double* x, double* y,
               VICcut** fc_cut ) {
/* Description:
   --------------
   Heuristic of Van Roy and Wolsey for finding a violated extended flow cover
   inequality for the capacitated facility location problem/single node flow problem.
   Let D = \sum_{i=1}^m demand[i] denote the total demand and N be the set of
   potential depot sites. Let C be a subset of N, which covers total demand,
   that is E = \sum_{j \in C} capaci[j] - D > 0. An extended flow cover
   inequality is then given by:

      \sum_{j\in C} ( x[j] + max{0,capaci[j] - E} ) * (1-y[j])
     +\sum_{j\in L} ( x[j] - (capaci*[j] - E)* y[j] ) <= D,
   where
      L is a subset of N\C
      capaci*[j] = max{ capaci[j], max{capaci[j] : j \in C} }
      y[j] is a location variable corresponding to depot j
      x[j] is the flow from depot j to the customers.

   The separation heuristic proceeds as follows: Let xlp, ylp denote the
   (fractional) solution values for the flow and location variables:
   (1) Determine a (minimal) cover C
   (2) Select L as L = { j \in N\C : xlp[j] > capaci*[j] - E }
   (3) Check if the resulting inequality is violated.

   For details see:
   -  Padberg M, Van Roy TJ, Wolsey LA (1985). Valid linear inequalities
      for fixed charge problems. Operations Research 33, 842-861
   -  Van Roy TJ, Wolsey LA (1987). Solving mixed integer programming
      problems using automatic reformulation. Operations Research 35,
      45-57


   Scope: Export
   -------------

   Parameters:
   -------------
   - n      : number of potential depot sites
   - m      : number of customers
   - demand : pointer to an array of integers of size of at least m containing
              the customers' demands
   - capaci : pointer to an array of integers of size of at least n containing
              the depot capacities
   - x      : pointer to an array of doubles of size of at leat m*n containing
              the allocation part of the fractional solution, which should be
	      separated by a flow cover inequality.
	      Let i=0,...,m-1 and j=0,...n-1 be the indices of customers and
	      depot sites, resp. Then x[i*n + j] is the solution value of
	      the allocation variable x(i,j), where 0 <= x(i,j) <= 1. The
	      variable x(i,j) denotes the fraction of customer i's demand met
	      from facility j.
   - y      : pointer to an array of doubles of size of at least n containing
              the location part of the fractional solution, which should be
	      separated by a flow cover inequality. 0 <= y[j] <= 1
   - fc_cut : pointer to a pointer to a cut structure VICcut defined in
              vico.h. If no violated flow cover is found, the null
	      pointer is returned; otherwise **fc_cut contains the cut.
*/
  int     mn=m*n, TotDem;
  int     bigcap, cap, card, lastl, depidx, excess, i, j, k, numnz;
  int     *cover = (int*) calloc( n, sizeof(int) );
  double  *flow = (double*) calloc( n, sizeof(double) );
  double  lhs, tmp;
  char    covviol;
  VICcut* cut=NULL;

  if ( (cover==NULL) || (flow==NULL) ) goto RETURN;
  for ( i=0,TotDem=0; i < m; i++ ) TotDem += demand[i];

  /* Determine a (violated) cover */
  covviol = VICcovsep( n, TotDem, capaci, y, cover, &card );
  if ( card == 0 ) goto RETURN;

  /* Determine flows of depots to customers in solution x */
  for ( j=0; j < n; j++ )
    for ( i=0, flow[j]=0.0; i < m; i++ ) flow[j] += x[i*n+j]*demand[i];


  /* Get the excess and largest capacity of depots in the cover */
  for ( j=0, bigcap=0, excess=-TotDem; j < card; j++ ) {
    depidx = cover[j];
    excess += capaci[depidx];
    if ( capaci[depidx] > bigcap ) bigcap = capaci[depidx];
  }

  /* Select the set L. Store the set in cover[card, ..., lastl-1].
     Start accumulating the left-hand side value and number of
     non-zeros of the resulting inequality */
  lastl = n, j=card, lhs=0.0;
  while ( j < lastl) {
    depidx = cover[j];
    cap = MAX( bigcap, capaci[depidx] );
    tmp = flow[depidx] - (cap - excess)*y[depidx];
    if ( tmp > VICtol ) lhs += tmp, j++;
    else cover[j] = cover[--lastl], cover[lastl] = depidx;
  }

  /* Count non-zeros and lhs for variables in cover C */
  numnz = m*lastl + lastl - card;
  for ( j=0; j < card; j++ ) {
    depidx = cover[j];
    lhs += flow[depidx];
    if ( capaci[depidx] > excess ) {
      lhs += (capaci[depidx]-excess)*(1.0-y[depidx]);
      numnz++;
    }
  }
  covviol = lhs*(1.0-VICtol) > (double) TotDem;

  /* If violated, return the inequality in structure **fc_cut */
  if ( covviol ) covviol = VICallocut( numnz, &cut );
  if ( covviol ) {
    cut->rhs = (double)TotDem;
    for ( i=0, k=0, numnz=0; i < m; i++ ) {
      for ( j=0; j < lastl; j++ ) {
        depidx = cover[j];
	cut->nzval[numnz] = demand[i];
	cut->nzind[numnz++] = k+depidx;
      }
      k += n;
    }
    for ( j=0; j < card; j++ ) {
      depidx=cover[j];
      if ( capaci[depidx] > excess ) {
        tmp = (double)( excess - capaci[depidx] );
	cut->rhs += tmp;
        cut->nzval[numnz] = tmp;
	cut->nzind[numnz++] = mn+depidx;
      }
    }
    for ( j=card; j < lastl; j++ ) {
      depidx = cover[j];
      cap = MAX( bigcap, capaci[depidx] );
      cut->nzval[numnz] = (double)(excess-cap);
      cut->nzind[numnz++] = mn+depidx;
    }
  }

RETURN:
  if ( cover ) free( cover );
  if ( flow ) free( flow );
  *fc_cut = cut;
}

/*------------------------------------------------------------------------------*/

void VICcflsmi( int m, int n, int* demand, int* capaci, double* x, double* y,
                VICcut** first_smi ) {
/*
   Description:
   -------------
   Aardal's separation heuristic for finding violated submodular inequalties
   for the capacitated facility location problem. See

   - Aardal K (1998). Capacitated facility location: Separation algorithms and
     computational experience. Mathematical Programming 81, 149-175
   - Aardal K, Pochet Y, Wolsey LA (1995). Capacitated facility location: Valid
     inequalities and facets. Mathematics of Operations Research 20, 552-582

   Aardal, Pochet, Wolsey propose the special cases of submodular
   valid inequalities for the CFLP. These inequalities, called
   effective capacity and single-depot inequalities, generalize flow
   cover inequalities for the CFLP.

   Scope: Export
   -------------

   Parameters:
   -------------
   - n         : number of potential depot sites
   - m         : number of customers
   - demand    : pointer to an array of integers of size of at least m containing
                 the customers' demands
   - capaci    : pointer to an array of integers of size of at least n containing
                 the depot capacities
   - x         : pointer to an array of doubles of size of at leat m*n containing
                 the allocation part of the fractional solution, which should be
	         separated by a submodular inequality.
	         Let i=0,...,m-1 and j=0,...n-1 be the indices of customers and
	         depot sites, resp. Then x[i*n + j] is the solution value of
	         the allocation variable x(i,j), where 0 <= x(i,j) <= 1. The
	         variable x(i,j) denotes the fraction of customer i's demand met
	         from facility j.
   - y         : pointer to an array of doubles of size of at least n containing
                 the location part of the fractional solution, which should be
	         separated by a submodular inequality. 0 <= y[j] <= 1
   - first_smi : pointer to the first cut (structure VICcut defined in vico.h)
                 of a linked list of violated submodular inequalities found
	         by this procedure. The NULL pointer is returned if no violated
	         submodular inequality was found.
*/

  int lfrac,    /* last fractional location variable used to start heuristic   */
      nFSet,    /* number of fractional location varaibles                     */
      npos,     /* number of positive location variables                       */
      nJ0Set,   /* number of depots in set J during loop                       */
      nJSet,    /* number of depots in set J after reduction                   */
      nK0Set,   /* number of customers in set K during loop                    */
      nKSet,    /* number of customers in set K after reduction                */
      nQSet,    /* number of depots in set Q                                   */
      nPSet,    /* number of depots in set P                                   */
      K0dem,    /* demand of set K during loop                                 */
      Kdem,     /* demand of set K after reduction                             */
      J0cap,    /* capacity of depot set J during loop                         */
      Jcap,     /* capacity of depot set J after reduction                     */
      excess;   /* excess capacity of depot set J over demand of K             */

  int *FSet=NULL,    /* set of depots with fractional LP value                 */
      *J0Set=NULL,   /* depot set J during loop                                */
      *JSet=NULL,    /* depot set J after reduction                            */
      *PDep=NULL,    /* PDep[k]=p if k uniquely supplied by a "p-depot"        */
      *K0Set=NULL,   /* customer set K during loop                             */
      *KSet=NULL,    /* customer set K after reduction                         */
      *Pflow=NULL,   /* unique demands to be met by depots in set J0           */
      *J0flow=NULL,  /* demands to be met by depots in set J0                  */
      *Jflow=NULL;   /* demands to be met by depots in set J                   */

  double *flow0=NULL,/* flow through depots                                    */
         *flow=NULL; /* flow through depots after reduction                    */

  char *IsPDep=NULL; /* IsPDep[j] > 0 if j is a "p-depot"                      */

  TJSet *JRoot=NULL, /* pointer to first generated depot set J                 */
        *JEnd=NULL;  /* pointer to last generated depot set J                  */

  int    col, fail, i, ii, j, jj, jsel, k, ncard, num, nserved;
  int    Rdem, Rcap, nPKSet, minflow, newexc;
  double coeff, lhs, sigma, minsigma;
  char   equal;
  TJSet  *AJSet=NULL;
  VICcut *cur_smi=NULL;

  *first_smi = NULL;

  /* Allocate required memory */
  FSet  = (int*) calloc( n, sizeof(int) ); if (FSet==NULL) goto RETURN;
  J0Set = (int*) calloc( n, sizeof(int) ); if (J0Set==NULL)goto RETURN;
  JSet  = (int*) calloc( n, sizeof(int) ); if (JSet==NULL) goto RETURN;
  PDep  = (int*) calloc( m, sizeof(int) ); if (PDep==NULL) goto RETURN;
  K0Set = (int*) calloc( m, sizeof(int) ); if (K0Set==NULL) goto RETURN;
  KSet  = (int*) calloc( m, sizeof(int) ); if (KSet==NULL) goto RETURN;
  Pflow = (int*) calloc( n, sizeof(int) ); if (Pflow==NULL) goto RETURN;
  J0flow= (int*) calloc( n, sizeof(int) ); if (J0flow==NULL) goto RETURN;
  Jflow = (int*) calloc( n, sizeof(int) ); if (Jflow==NULL) goto RETURN;
  flow0 = (double*) calloc( n, sizeof(double) ); if (flow0==NULL) goto RETURN;
  flow  = (double*) calloc( n, sizeof(double) ); if (flow==NULL) goto RETURN;
  IsPDep= (char*) calloc( n, sizeof(char) ); if (IsPDep==NULL) goto RETURN;

  /* Get flows through depots, set of fractional location variables, etc */
  for ( j=0, npos=0, nFSet=0; j < n; j++ ) {
    J0flow[j] = 0, flow0[j] = 0.0;
    if ( y[j] > VICtol ) {
      J0Set[npos++] = j, col=j;
      for ( i=0; i < m; i++ ){
         flow0[j] += x[col]*demand[i];
	 if ( x[col] > VICtol ) J0flow[j] += demand[i];
	 col += n;
      }
      if ( y[j] < VICone ) FSet[nFSet++] = j;
    }
  }

  /* Apply Aardal's heuristic for identifying set J and K */
  lfrac=0;

NEWJSET:
  if ( lfrac >= nFSet ) goto RETURN;
  fail = 0;
  /* Initialize subset J with next depot of fractional value
     and customer subset K with customers served by this depot */
  for ( j=0; j < n; j++ ) if ( FSet[lfrac] == J0Set[j] ) break;
  VICswap( J0Set, J0Set+j );
  j = J0Set[0], nJ0Set = 1, J0cap = MIN( capaci[j], J0flow[j] );
  for ( i=0, nK0Set=0, col=j, ii=m; i < m; i++ ) {
    if ( x[col] > VICtol ) K0Set[nK0Set++] = i;
    else  K0Set[--ii] = i;
    col += n;
  }
  K0dem = J0flow[j], lfrac++;

SELNEXT:
  /* Augment set J with depots j of positive LP-value and small value of
     min{ d(Kj), capaci[j] } - flow[j], where Kj is the set of customers
     served by j and d(Kj) its demand. Augment set K by sets Kj.
     Repeat until J covers K */
  do {
    for ( j=nJ0Set, jsel=-1, minsigma=VICinf; j < npos; j++ ) {
      jj = J0Set[j];
      sigma = MIN( J0flow[jj], capaci[jj] ) - flow0[jj];
      if ( sigma < minsigma )
	minsigma = sigma, jsel=j;
      else if ( (!(sigma > minsigma)) && (flow0[jj] > flow0[J0Set[jsel]]) )
	minsigma = sigma, jsel=j;
    }
    if ( jsel >= 0 ) {
      VICswap( J0Set+jsel, J0Set+nJ0Set );
      jsel  = J0Set[nJ0Set++];
      J0cap+= MIN( capaci[jsel], J0flow[jsel] );
      ii = m;
      while ( nK0Set < ii ) {
        k = K0Set[nK0Set];
	if ( x[k*n+jsel] > VICtol ) K0dem += demand[k], nK0Set++;
	else ii--, VICswap( K0Set+nK0Set, K0Set+ii );
      }
    }
  } while ( ( J0cap <= K0dem ) && (jsel >= 0 ) );

  /* If set J could not be enlarged try a new set J of depots */
  if ( jsel < 0 ) goto NEWJSET;

  /* Check if set J=J0 has already been investigated */
  equal=0, AJSet=JRoot;
  while ( (AJSet != NULL) && ( !(equal) ) ) {
    equal = ( AJSet->card == nJ0Set );
    for (j=0; (j<nJ0Set)&&(equal); j++ ) equal=(AJSet->inset[J0Set[j]]);
    AJSet = AJSet->succ;
  }
  if ( equal ) {
    if ( nK0Set < m ) goto SELNEXT; else goto NEWJSET;
  }

  /* Copy sets J0 and K0 found so far into sets J and K for
     cut generation. Store depot set J in linked list */
  AJSet = (TJSet*) malloc( sizeof( TJSet ) );
  if ( AJSet == NULL ) goto RETURN;
  AJSet->card = nJ0Set;
  AJSet->succ = NULL;
  AJSet->inset= (char*) calloc(n, sizeof(char) );
  if (JRoot==NULL) JRoot=JEnd=AJSet; else JEnd->succ = AJSet; JEnd=AJSet;
  if ( AJSet->inset == NULL ) goto RETURN;
  Jcap=J0cap, Kdem= K0dem, nJSet=nJ0Set, nKSet=nK0Set;
  for ( i=0; i < nKSet; i++ ) KSet[i]=K0Set[i];
  for ( j=0; j < n; j++ ) AJSet->inset[j] = 0;
  for ( j=0; j < nJSet; j++ ) {
    jj = J0Set[j];
    AJSet->inset[jj] = 1;
    JSet[j]=jj, Jflow[jj]=J0flow[jj], flow[jj]=flow0[jj];
  }
  goto GETQSET;


REMCUST:
  /* Remove customers KSet[ncard, ... ,nKSet-1] from set K */
  for ( i=ncard; i < nKSet; i++ ) {
     ii = KSet[i], col = ii*n;
     Kdem -= demand[ii];
     for ( j=0; j < nJSet; j++ ) if ( x[JSet[j]+col] > VICtol ) {
       jj = JSet[j];
       Jflow[jj] -= demand[ii];
       flow[jj] -= x[jj+col]*demand[ii];
     }
  }
  nKSet = ncard;

  /* Move depots with effective capacity smaller than capacity
     to the set Q of depot */
GETQSET:
  nQSet=j=Jcap=0, ii=nJSet;
  while ( j < ii ) {
    jj = JSet[j];
    if ( Jflow[jj] < capaci[jj] ) {
      ii--, VICswap( JSet+j, JSet+ii );
      if ( Jflow[jj] > 0 ) nQSet++, Jcap += Jflow[jj];
      else nJSet--, VICswap( JSet+ii, JSet+nJSet );
    }
    else Jcap += capaci[jj], j++;
  }
  excess = Jcap - Kdem;
  if ( excess <= 0 ) goto SELNEXT;

  /* Remove customers served by more than one q-depot from
     set K and obtain new set Q */
  i=0, ncard = nKSet;
  if ( nQSet > 1 ) {
    while ( i < ncard ) {
      for ( j=nJSet-nQSet, nserved=0; (j<nJSet)&&(nserved<2); j++ )
	if ( x[KSet[i]*n+JSet[j]] > VICtol ) nserved++;
     if ( nserved > 1 ) ncard--, VICswap( KSet+i, KSet+ncard );
     else i++;
    }
    if ( ncard < nKSet ) goto REMCUST;
  }

  /* Remove depots with effective capacity less than the excess
     from the set Q and J. Remove customers (solely) served by
     a removed q-depot from the set K */
  num = nJSet, j=nJSet-nQSet;
  while ( j < num ) {
    if ( Jflow[JSet[j]] <= excess ) num--, VICswap(JSet+j,JSet+num);
    else j++;
  }
  if ( num < nJSet ) {
    i=0, ncard = nKSet;
    while ( i < ncard ) {
      col = KSet[i]*n;
      for ( j=num; j < nJSet; j++ ) if ( x[col+JSet[j]] > VICtol ) break;
      if ( j < nJSet ) ncard--, VICswap( KSet+i, KSet+ncard );
      else i++;
    }
    nJSet = num;
    if ( ncard < nKSet ) goto REMCUST;
  }

  /* Apply Aardal's heuristic for identifying a single-depot structure:*/
  for ( i=0; i < nKSet; i++ ) PDep[KSet[i]] = -1;
  for ( j=0; j < nJSet; j++ ) IsPDep[JSet[j]] = 0, Pflow[JSet[j]] = 0;
  nPSet = 0;
  if ( nJSet > 2 ) {  /* get candidate set of p-depots */
    Rdem = Kdem, Rcap = Jcap, num = nQSet, nPKSet = 0;
    for ( i=0; i < nKSet; i++ ) {
      ii=KSet[i], col=ii*n;
      for ( j=0, nserved=0; (j < nJSet)&&(nserved < 2); j++ )
        if ( x[col+JSet[j]] > VICtol ) PDep[ii]=JSet[j], nserved++;
      if ( nserved > 1 )
        PDep[ii] = -1;
      else {
        jj=PDep[ii], nPKSet++, Pflow[jj] += demand[ii];
	IsPDep[jj] = ( capaci[jj] < Kdem );
      }
    }
    if ( nPKSet > 0 ) for ( j=0; j < nJSet; j++ ) {
      jj = JSet[j];
      if ( IsPDep[jj] ) IsPDep[jj] = ( Pflow[jj] < capaci[jj] );
      if ( IsPDep[jj] ) {
        Rdem -= Pflow[jj], Rcap -= MIN( capaci[jj], Jflow[jj] ), nPSet++;
	if ( capaci[jj] > Jflow[jj] ) num--;
      }
    }
    while ( (Rcap <= Rdem) || (nJSet < nPSet+num) ) {
      /* reduce candidate set of p-depots until single-depot structure */
      for ( j=0, minflow=VICbig; j < nJSet; j++ ) if ( IsPDep[JSet[j]] ) {
	if ( Pflow[JSet[j]] < minflow ) jj=JSet[j], minflow = Pflow[jj];
      }
      Rcap += MIN(capaci[jj],Jflow[jj]), Rdem += Pflow[jj], nPSet--;
      if ( capaci[jj] > Jflow[jj] ) num--;
      IsPDep[jj] = 0;
    }
    for ( i=0; i < nKSet; i++ ) {
      ii = KSet[i], jj = PDep[ii];
      if ( (jj >= 0) && ( ! (IsPDep[jj]) ) ) PDep[ii] = -1;
    }
    newexc = -Kdem;
    if ( nPSet > 0 ) { /* calculate new excess of capacity over demand */
      for ( j=0; j < nJSet; j++ ) {
        jj = JSet[j];
	if ( IsPDep[jj] ) newexc += MIN(capaci[jj],Rdem+Pflow[jj]);
	else newexc += MIN(capaci[jj],Jflow[jj]);
      }
    }
    if ( (newexc > excess)&& (num > 0) ) {
      /* check if q-depots have noncero coefficients */
      for ( j=nJSet - nQSet; j < nJSet; j++ ) if ( ! (IsPDep[JSet[j]]) ) {
        if ( Jflow[JSet[j]] <= newexc ) break;
      }
      if ( j < nJSet ) {
        for ( i=0; i < nKSet; i++ ) PDep[KSet[i]] = -1;
	for ( jj=0; jj < nJSet; jj++ ) IsPDep[JSet[jj]] = 0, Pflow[JSet[jj]] = 0;
	nPSet = 0;
      }
    }
    if ( nPSet > 0 ) excess = newexc;
  }

  /* Check if the inequality corresponding to the constructed effective
     capacity/single-depot structure is violated */
  for ( j=0, lhs=0.0; j < nJSet; j++ ) {
    jj = JSet[j];
    if ( IsPDep[jj] ) coeff = (double)Pflow[jj];
    else if ( capaci[jj] > Jflow[jj] ) coeff = Jflow[jj] - excess;
    else coeff = MAX( 0.0, capaci[jj] - excess );
    lhs += flow[jj] + coeff*(1.0-y[jj]);
  }
  if ( lhs*VICone > Kdem ) { /* inequality is violated->add to list */
    if ( ! ( VICallocut( (nKSet+1)*nJSet, &cur_smi ) ) ) goto RETURN;
    cur_smi->rhs = (double) Kdem;
    num = 0;
    for ( j=0; j < nJSet; j++ ) {
      jj = JSet[j];
      if ( IsPDep[jj] ) { /* coefficients p-customers */
        coeff = (double)(-Pflow[jj]);
        for ( i=0; i < nKSet; i++ ) {
	  ii = KSet[i];
	  if ( (PDep[ii] < 0) || ( PDep[ii] == jj ) ) {
	    cur_smi->nzval[num] = (double)demand[ii];
	    cur_smi->nzind[num++] = ii*n+jj;
	  }
	}
      }
      else if ( Jflow[jj] < capaci[jj] ) { /* coeff q-customers */
        coeff = (double)(excess - Jflow[jj]);
	for ( i=0; i < nKSet; i++ ) {
	  ii = KSet[i], col=ii*n+jj;
	  if ( x[col] > VICtol ) {
	    cur_smi->nzval[num] = (double)demand[ii];
	    cur_smi->nzind[num++] = col;
	  }
	}
      }
      else { /* coefficients other customers in KSet */
        coeff = (double) ( MIN( 0.0, excess-capaci[jj] ) );
	for ( i=0; i < nKSet; i++ ) if ( PDep[KSet[i]] < 0 ) {
	  cur_smi->nzval[num] = (double)demand[KSet[i]];
	  cur_smi->nzind[num++] = KSet[i]*n+jj;
	}
      }
      /* coefficient of depot jj=JSet[j] */
      if ( coeff < 0.0 ) {
        cur_smi->nzval[num] = coeff;
	cur_smi->nzind[num++] = m*n+jj;
	cur_smi->rhs += coeff;
      }
    }
    if ( num < cur_smi->nzcnt ) {
      if ( ! ( VICshrinkcut( num, cur_smi ) ) ) goto RETURN;
    }
    VICaddtolst( first_smi, cur_smi );
  }
  else fail++;

  if ( ( fail==2 ) || ( nK0Set == m ) ) goto NEWJSET; else goto SELNEXT;

RETURN:
  if ( FSet ) free( FSet );
  if ( J0Set ) free( J0Set );
  if ( JSet ) free( JSet );
  if ( PDep ) free( PDep );
  if ( K0Set ) free( K0Set );
  if ( KSet ) free( KSet );
  if ( Pflow ) free( Pflow );
  if ( J0flow ) free( J0flow );
  if ( Jflow ) free( Jflow );
  if ( flow0 ) free( flow0 );
  if ( flow ) free( flow );
  if ( IsPDep ) free( IsPDep );
  while ( JRoot != NULL ) {
    JEnd = JRoot;
    JRoot = JEnd->succ;
    if ( JEnd->inset ) free( JEnd->inset );
    free( JEnd );
  }
}

/*------------------------------------------------------------------------------*/

void VICuflohi( int m, int n, double* x, double* y, VICcut** first_ohi ) {
/*
  Description:
  -------------
  Generates odd-hole inequalities for the UFLP violated by the fractional
  solution (x,y) where y = ( y(j) ) denotes the location variables
  (0 <= y(j) <= 1) and x = ( x(i,j) ) denotes the allocation variables
  (0 <= x(i,j) <= y(j) ). Define z(j) = 1-y(j). Odd-hole inequalities
  for the UFLP have the following structure:

     z(j1) + x(i1,j1) + x(i1,j2) + z(j2) + x(i2,j2) + x(i2,j3) + z(j3)
   + ... + x(ir,jr) + x(ir,j1) <= 3r/2 - 1

  where r >= 3 and odd. The inequality can only be violated if all
  variables appearing in the inequality have a fractional value.
  For details see
  - Cornuejols G, Thizy J-M (1982). Some facets of the simple plant
    location polytope. Mathematical Programming 23, 50-74
  - Caprara A, Salazar JJG (1999). Separating lifted odd-hole inequalities
    to solve the index selection problem. Disrete Applied Mathematics 92,
    111-134

  Scope: Export
  -------------

  Parameters:
  -------------
   - n         : number of potential depot sites
   - m         : number of customers
   - x         : pointer to an array of doubles of size of at leat m*n containing
                 the allocation part of the fractional solution.
	         Let i=0,...,m-1 and j=0,...n-1 be the indices of customers and
	         depot sites, resp. Then x[i*n + j] is the solution value of
    	         the allocation variable x(i,j), where 0 <= x(i,j) <= 1. The
	         variable x(i,j) denotes the fraction of customer i's demand met
	         from facility j.
   - y         : pointer to an array of doubles of size of at least n containing
                 the location part of the fractional solution, which should be
	         separated by a submodular inequality. 0 <= y[j] <= 1, y[j]>=x[i,j]
   - first_ohi : pointer to the first cut (structure VICcut defined in vico.h)
                 of a linked list of violated odd-hole inequalities found
	         by this procedure. The NULL pointer is returned if no violated
	         inequality was found.
*/


  int     nedges,       /* number of edges of auxiliary bipartite graph        */
          nnodes,       /* number of nodes of auxiliary bipartite graph        */
	  onodes,       /* number of nodes on odd cycle                        */
          nfract;       /* number of depots with fractional LP-value           */

  int    *Custs=NULL,  /* indices of customer nodes in artificial graph        */
         *FSet=NULL,   /* set of depots with fractional LP-value               */
         *Queue=NULL,  /* queue of nodes used in shortest-path computation     */
	 *Pred=NULL,   /* predecessor list to store a path                     */
	 *Cycle=NULL,  /* location variable nodes on an odd cycle              */
	 *Cedge=NULL;  /* customers corresponding to an edge of a cycle        */
  double *Dist=NULL,   /* length of edges of bipartite artificial graph        */
         *Len=NULL,    /* length of shortest paths                             */
         *Weight=NULL; /* weight of odd-cycle starting and ending in depot j   */

  int     cust, edge, i, ii, j, jj, k, last, next, node, root;
  int     fn1, fn2, en1, en2, n1, n2, mn=m*n;
  double  ndist, owc, wc, ws, xs;
  char    ok, IsDble;
  int    *nzind;
  double *nzval;
  double  LenEps = 1.0E-8;
  VICcut* cur_ohi;


  *first_ohi = NULL;
  FSet = (int*) calloc( n, sizeof(int) ); if ( FSet==NULL ) goto RETURN;
  Weight = (double*) calloc(n, sizeof(double) ); if (Weight==NULL) goto RETURN;
  Cycle = (int*) calloc(n+1, sizeof(int) ); if (Cycle==NULL) goto RETURN;
  Cedge = (int*) calloc(n+1, sizeof(int) ); if (Cedge==NULL) goto RETURN;

  /* Determine set of depots with fractional LP-value */
  for ( j=0,nfract=0; j < n; j++ ) {
    Weight[j] = VICone;
    if ( (y[j] > VICtol) && (y[j] < VICone) ) FSet[nfract++] = j;
  }

  if ( nfract < 3 ) goto RETURN;

  /* Allocate memory required for artificial bipartite graph */
  nnodes = 2*nfract;
  nedges = ( nfract*(nfract+1) )/2;
  Dist   = (double*) calloc( nedges, sizeof(double) ); if (Dist==NULL) goto RETURN;
  Custs  = (int*) calloc( nedges, sizeof(int) ); if (Custs==NULL) goto RETURN;
  Queue  = (int*) calloc( nnodes+1, sizeof(int) ); if (Queue==NULL) goto RETURN;
  Pred   = (int*) calloc( nnodes, sizeof(int) ); if (Pred==NULL) goto RETURN;
  Len    = (double*) calloc( nnodes, sizeof(double)); if (Len==NULL) goto RETURN;

  /* Get length/weights of edges of auxialary bipartite graph */
  for ( i=0,edge=0; i < nfract; i++ ) {
    for ( j=0; j < i; j++ ) {
      ii = FSet[i], jj = FSet[j];
      for ( k=0, ws=0.0; k < m; k++ ) {
        xs = x[ii] + x[jj];
	if ( xs > ws ) ws=xs, cust=k;
	ii += n, jj += n;
      }
      Dist[edge] = 1.0 + y[FSet[i]] + y[FSet[j]] - 2.0*ws;
      if ( Dist[edge] < 0.0 ) Dist[edge] = 0.0;
      Custs[edge++] = cust;
    }
    Dist[edge] = VICinf;
    Custs[edge++] = 0;
  }

  /* Compute odd cycles by means of computing shortest paths in the
     auxialary bipartite graph, starting with each fractional
     location variable "node" in turn */
  for ( root=0; root < nfract; root++ ) {

    /* compute minimum weight odd cycle with root "root" */
    owc=Weight[root];
    /* initialize queue of nodes whose forward arcs have to be explored */
    for ( i=0; i < nnodes; i++ ) Queue[i] = -1, Len[i]=VICinf;
    Queue[nnodes] = root, Queue[root] = nnodes, last = root, Len[root]=0.0;
    /* As long as queue is not empty: */
    while ( Queue[nnodes] < nnodes ) {
      /* pick node from head of queue */
      node = Queue[nnodes], Queue[nnodes]=Queue[node], Queue[node]=-1;
      if ( last == node ) last = nnodes;
      /* explore forward star of node "node" */
      i = node % nfract;
      if ( node < nfract ) k=nfract; else k=0;
      for ( j=0; j< nfract; j++ ) {
        next = j+k;
  	ndist = Dist[ GETEDGE(i,j) ] + Len[node];
	if ( ndist+LenEps < Len[next] ) {
	  Len[next] = ndist;
	  Pred[next] = node;
	  /* add node next at tail of the queue */
	  if ( Queue[next] < 0 ) Queue[last] = next, Queue[next]=nnodes, last=next;
	}
      }
    }

    /* check if path from node "root" to its copy node "nfract+root
       gives an odd cylce without node repetitions */
    node = nfract+root;
    wc = Len[node];
    Cycle[0] = root;
    onodes = 0;
    while ( node != root ) {
      next = Pred[node] % nfract;
      Cycle[++onodes] = next;
      Weight[next] = MIN( Weight[next], wc );
      node = Pred[node];
    }
    for ( node=0, ok=1; (node<onodes)&&(ok); node++ ) {
      for ( next=node+1; (next<onodes)&&(ok); next++ ) ok=(Cycle[node]!=Cycle[next]);
      Cedge[node] = Custs[ GETEDGE(Cycle[node], Cycle[node+1] ) ];
    }
    if ( ok ) ok = (wc < owc ); /* weight of cycle less than 1 */

    /* Check cycle for customer repetitions */
    if ( ok ) {
      do {
        for ( fn1=0, IsDble=0; (fn1<onodes)&&(IsDble==0); fn1++ )
          for ( fn2=fn1+1; (fn2<onodes)&&(IsDble==0); fn2++ )
	    IsDble = (Cedge[fn1] == Cedge[fn2]);
        if ( IsDble ) {  /* remove customer repition */
          en1=fn1--;
	  en2=fn2--;
  	  n1 = (fn1+1)+(onodes-en2);
	  n2 = fn2-en1+1;
	  if ( (n1 % 2) > 0 ) {
	    for ( node=en2; node <= onodes; node++ )
	      Cycle[node-en2+en1] = Cycle[node], Cedge[node-en2+en1] = Cedge[node];
	    onodes = n1;
	  }
	  else {
	    for ( node = en1; node <= fn2; node++ )
	      Cycle[node-en1] = Cycle[node], Cedge[node-en1] = Cedge[node];
	    Cycle[n2] = Cycle[0];
	    onodes = n2;
	  }
        }
      } while ( IsDble );
    }

    /* Check if inequality corresponding to computed cycle was already generated */
    cur_ohi = *first_ohi;
    while ( (ok) && (cur_ohi != NULL) ) {
      for ( node=0; node < onodes; node++ ) {
	j = FSet[Cycle[node]]+mn;
	for ( jj=0; jj < cur_ohi->nzcnt; jj++ ) if ( j==cur_ohi->nzind[jj] ) break;
	if ( jj==cur_ohi->nzcnt ) break;
      }
      ok = ( node < onodes );
      cur_ohi = cur_ohi->nextcut;
    }

    /* Add the corresponding inequality to the list of generated inequalities */
    if ( ok ) {
      if ( ! ( VICallocut( 3*onodes, &cur_ohi ) ) ) goto RETURN;
      cur_ohi->rhs = (double) ( (onodes*3 )/2 - onodes );
      nzval = cur_ohi->nzval;
      nzind = cur_ohi->nzind;
      for ( node=0; node < onodes; node++ ) {
        j=FSet[Cycle[node]], jj=FSet[Cycle[node+1]], cust=Cedge[node];
	*nzval++ = -1.0, *nzval++=1.0, *nzval++=1.0;
	*nzind++ = mn+j, *nzind++=cust*n+j, *nzind++=cust*n+jj;
      }
      VICaddtolst( first_ohi, cur_ohi );
    }

  }   /* end for root */


RETURN:
  if ( Custs ) free(Custs);
  if ( FSet ) free(FSet);
  if ( Queue ) free(Queue);
  if ( Pred ) free(Pred);
  if ( Dist ) free(Dist);
  if ( Len ) free(Len);
  if ( Weight ) free(Weight);
  if ( Cycle ) free(Cycle);
  if ( Cedge ) free(Cedge);

}

/*------------------------------------------------------------------------------*/

void VICgapfn( int m, int n, int IsGap, int* capaci, int** weight, double* X,
               int* indx, VICcut** fn_cut ) {
/*
  Description:
  -------------
  Given the GAP

       max   \sum_{i\in I} \sum_{j\in J} c_{ij} x_{ij}

       s.t.:         \sum_{i\in I} x_{ij} = 1   \forall j \in J \\

             \sum_{j\in J} w_{ij} x_{ij} <= s_i \forall i \in I \\

	             x_{ij} \in \{0,1\} \forall i,j

  the procedure uses a heuristic described in Farias IR, Nemhauser GL (2001),
  A family of inequalities for the generalized assignment polytope, Oper. Res.
  Letters 29, for finding a violated inequality of the form

  \sum_{j\in J_i} w_{ij} x_{ij} + \sum_{j\in J_i} a_j \sum_{k\ne i} x_{kj} <= s_i

  where J_i \subseteq J, a_j = s_i - ( W_i - w_{ij} ), W_i = \sum_{j\in J_i} w_{ij}

  Scope: Export
  -------------

  Parameters:
  -------------
  - m      : number of agents in the GAP
  - n      : number of jobs in the GAP
  - IsGap  : equals 1 if the problem under consideration is a GAP and 0 if it is
             a LEGAP
  - capaci : pointer to an array of integers of length of at least m containing
             the agents' capacities
  - weight : weight[i][j] is the amount w_{ij} of resources required by agent i
             for performing job j
  - X      : pointer to an array of doubles of length of at least m*n, such that
             X[i*n+j] gives the LP value of variable x_{ij} for i=0,...,m-1
	     and j=0,...,n-1
  - indx   : if not null, indx[i*n+j] contains the index of assignment variable
             x(i,j) for i=0,...,m-1 and j=0,...,n-1
  - fn_cut : pointer to the first cut (structure VICcut defined in vico.h)
             of a linked list of violated inequalities found
             by this procedure. The NULL pointer is returned if no violated
	     inequality was found.
*/
  int     alpha, coeff, i, j, k, num, numnz, wsum, isviol;
  double  lhs, *xij;
  char    *Iserved = (char*) calloc(n, sizeof(char) );
  VICcut  *curcut;

  *fn_cut = NULL;
  for ( i=0; i < m; i++ ) {
    /* check if agent i's capacity is fully used in solution X */
    for (j=0, numnz=0, xij=X+i*n, lhs=0.0, wsum=0; j < n; j++, xij++ ) {
      Iserved[j] = ( *xij > VICtol );
      if ( Iserved[j] ) lhs += *xij*weight[i][j], wsum += weight[i][j], numnz++;
    }
    num = numnz;
    isviol = ( ( (double)capaci[i]-lhs ) < VICtol );
    if ( isviol ) {
      for ( j=0,isviol=0; j < n; j++) if ( Iserved[j] ) {
        coeff = capaci[i] - wsum + weight[i][j];
	if ( coeff > 0 ) {
	  numnz += m-1;
	  for ( k=0; (k < m)&&(isviol==0); k++ ) if ( k != i) {
	    isviol = ( X[k*n+j] > VICtol );
	  }
	}
      }
    }
    if ( isviol ) {
      /* Add the cut */
      if ( IsGap ) { /* In order to reduce number of nonzeros in the cut,
                        substitute \sum_{k\ne i} x_{kj} by 1 - x_{ij} */
        if ( VICallocut( num, &curcut ) ) {
	  curcut->sense = 'L';
	  curcut->rhs = (double)capaci[i];
	  for ( j=0, numnz=0; j < n; j++ ) if ( Iserved[j] ) {
	    alpha = MAX(0, capaci[i] - wsum + weight[i][j]);
	    coeff = weight[i][j] - alpha;
	    if ( coeff != 0 ) {
	      curcut->nzval[numnz] = (double) coeff;
	      curcut->nzind[numnz++] = (indx) ? indx[i*n+j] : i*n+j;
	    }
	    curcut->rhs -= (double) alpha;
	  }
	  curcut->nzcnt = numnz;
        }
      } else if ( VICallocut( numnz, &curcut ) ) { /* problem is a LEGAP! */
        for ( j=0,numnz=0; j < n; j++ ) if ( Iserved[j] ) {
	  curcut->nzval[numnz] = (double)weight[i][j];
	  curcut->nzind[numnz++] = (indx) ? indx[i*n+j] : i*n+j;
	  coeff = capaci[i] - wsum + weight[i][j];
	  if ( coeff > 0 ) for ( k=0; k < m; k++ ) if ( k != i ) {
	    curcut->nzval[numnz] = (double)coeff;
	    curcut->nzind[numnz++] = (indx) ? indx[k*n+j] : k*n + j;
	  }
	}
	curcut->sense = 'L';
	curcut->rhs = (double) capaci[i];
	curcut->nzcnt = numnz;
      }
      VICaddtolst( fn_cut, curcut );
    }
  }

  free(Iserved);

}

/*------------------------------------------------------------------------------*/

void VICuflcov( CPXENVptr Env, int m, int n, char SolveCov, double* x, double* y,
                VICcut** first_cut ) {
/*
  Description:
  ------------
  Let K be a subset of the set of all customers, and let J denote a subset
  of the set of potential depot sites. Define a binary matrix a(i,j) for
  each i \in K and j \in J. Let b denote (a lower bound on) the minimum number
  of depots j \in J required to cover each customer i \in K, that is

   b =   min \sum_{j\in J} y_j
       s.t.: \sum_{j\in J} a_{ij} y_j \ge 1 \forall i \in K
             y_j = 0,1 \forall j \in J

  Then the inequality

  \sum_{i \in K} \sum_{j\in J} a_{ij} x_{ij} - \sum_{j\in J} \le |K| - b

  are valid for the UFLP. These inequalites generalize odd holes for the UFLP
  and have been proposed by D.C. Cho et al. (1983). On the uncapacitated
  facility location problem I: Valid inequalities and facets. Mathematics
  of Operations Research 8, 579-589. See also G. Cornuejols, J.-M. Thizy
  (1982). Some facets of the simple plant location polytope. Mathematical
  Programming 23, 50-74.

  Let (x*, y*) denote a fractional solution. In order to find a violated
  inequality of the above type, the following simple heuristic is tried:
  (i)    Set J = \{ j : 0 < y*_j < 1 \}
  (ii)   Set K = \{ i : x*_{ij} > 0 for more than 1 j \in J \}
  (iii)  Solve the covering problem in order to find b (alternatively, find
         a lower bound on b by solving the LP relaxation of the covering
	 problem and rounding up the objective function value)
  (iv)   Check if point (x*,y*) violates the resulting inequality


  Scope: Export
  -------------

  Parameters:
  -------------
  - m, n     : number of customers and number of potential depot sites
  - SolveCov : if equal to 1, the covering problem is solved exactly;
               otherwise a lower bound on b is computed
  - x        : pointer to an array of doubles of size of at leat m*n containing
               the allocation part of the fractional solution.
	       Let i=0,...,m-1 and j=0,...n-1 be the indices of customers and
	       depot sites, resp. Then x[i*n + j] is the solution value of
    	       the allocation variable x(i,j), where 0 <= x(i,j) <= 1. The
	       variable x(i,j) denotes the fraction of customer i's demand met
	       from facility j.
  - y        : pointer to an array of doubles of size of at least n containing
               the location part of the fractional solution, which should be
	       separated by a submodular inequality. 0 <= y[j] <= 1, y[j]>=x[i,j]
  - first_ohi: pointer to the first cut (structure VICcut defined in vico.h)
               of a linked list of violated odd-hole inequalities found
	       by this procedure. The NULL pointer is returned if no violated
	       inequality was found.
*/

  double   *rhs=NULL, *matval=NULL, *obj=NULL, *lb=NULL, *ub=NULL, lhs, objval;
  int      *matbeg=NULL, *matcnt=NULL, *matind=NULL;
  char     *sense=NULL, isviol;
  char     *probname = (char*) calloc(16,sizeof(char) );
  CPXLPptr covlp=NULL;

  int *JSet = ( int* ) calloc( n, sizeof( int ) );
  int *KSet = ( int* ) calloc( m, sizeof( int ) );
  int col, covnum, i, j, mm, nn, num, numnz, nz, Status;

  /* Determine the customer and depot subsets defined above */
  *first_cut = NULL;
  for ( j=0, nn=0; j < n; j++ ) if ( (y[j] > VICtol) && (y[j] < VICone) ) JSet[nn++] = j;
  for ( i=0, mm=0, numnz=0; i < m; i++ ) {
    for ( j=0, num=0; j < nn; j++ ) num += ( x[i*n+JSet[j]] > VICtol );
    if ( num > 1 ) KSet[mm++] = i, numnz += num;
  }
  if (mm<2) goto RETURN;

  /* Solve the covering problem and its LP relaxation, resp. */
  strcpy( probname, "tmpcover");
  covlp = CPXcreateprob( Env, &Status, probname );
  if ( Status == 0 ) Status = ( (rhs = (double*) calloc(mm,sizeof(double)))==NULL );
  if ( Status == 0 ) Status = ( (sense = (char*) calloc(mm,sizeof(char)))==NULL );
  if ( Status == 0 ) Status = ( (obj = (double*) calloc(nn,sizeof(double)))==NULL );
  if ( Status == 0 ) Status = ( (lb = (double*) calloc(nn,sizeof(double)))==NULL );
  if ( Status == 0 ) Status = ( (ub = (double*) calloc(nn,sizeof(double)))==NULL );
  if ( Status == 0 ) Status = ( (matbeg = (int*) calloc(nn,sizeof(int)))==NULL );
  if ( Status == 0 ) Status = ( (matcnt = (int*) calloc(nn,sizeof(int)))==NULL );
  if ( Status == 0 ) Status = ( (matind = (int*) calloc(numnz,sizeof(int)))==NULL );
  if ( Status == 0 ) Status = ( (matval = (double*) calloc(numnz,sizeof(double)))==NULL );
  if ( Status > 0 ) goto RETURN;

  for ( i=0; i < mm; i++ ) rhs[i]=1.0, sense[i]='G';
  for ( j=0, nz=0; j < nn; j++ ) {
    lb[j] = 0.0, ub[j] = 1.0, obj[j] = 1.0;
    matbeg[j] = nz, matcnt[j] = 0;
    for ( i=0; i < mm; i++ ) if ( x[KSet[i]*n+JSet[j]] > VICtol ) {
      matcnt[j]++;
      matval[nz]=1.0;
      matind[nz++] = i;
    }
  }

  Status = CPXcopylp( Env, covlp, nn, mm, CPX_MIN, obj, rhs, sense, matbeg, matcnt,
                      matind, matval, lb, ub, NULL );
  if ( SolveCov ) {
    Status = CPXmipopt( Env, covlp );
    Status = CPXgetmipobjval( Env, covlp, &objval );
    covnum = (int)(objval + VICtol);
  } else {
    Status = CPXoptimize( Env, covlp );
    Status = CPXgetobjval( Env, covlp, &objval );
    if ( ( objval - floor(objval) ) < VICtol )
      covnum = (int) (objval + VICtol );
    else covnum = (int) (objval + 1.0);
  }

  /* Check if resulting inequality is violated */
  for ( i=0,lhs=0.0; i < mm; i++ )
    for ( j=0; j < nn; j++ ) {
      col = KSet[i]*n + JSet[j];
      if (x[col] > VICtol) lhs += x[col];
    }
  for ( j=0; j < nn; j++ ) lhs -= y[JSet[j]];
  isviol =  (lhs*VICone > mm-covnum);
  if ( ! (isviol) ) goto RETURN;

  /* Create the inequality */
  numnz += nn;
  if ( VICallocut( numnz, first_cut ) ) {
    (*first_cut)->sense = 'L';
    (*first_cut)->rhs = (double) ( mm-covnum);
    (*first_cut)->nzcnt = numnz;
    for (i=0,nz=0; i < mm; i++ )
      for ( j=0; j < nn; j++ ) {
        col = KSet[i]*n+JSet[j];
	if ( x[col] > VICtol ) {
	  (*first_cut)->nzval[nz]=1.0;
	  (*first_cut)->nzind[nz++]=col;
	}
      }
    for ( j=0; j < nn; j++ ) {
      col = n*m + JSet[j];
      (*first_cut)->nzval[nz] = -1.0;
      (*first_cut)->nzind[nz++] = col;
    }
  }


RETURN:
  if ( rhs ) free( rhs );
  if ( sense ) free( sense );
  if ( obj ) free( obj );
  if ( lb ) free( lb );
  if ( ub ) free( ub );
  if ( matbeg ) free( matbeg );
  if ( matcnt ) free( matcnt );
  if ( matind ) free( matind );
  if ( matval ) free( matval );
  if ( probname ) free( probname );
  if ( covlp ) Status = CPXfreeprob( Env, &covlp );
  if ( JSet ) free( JSet );
  if ( KSet ) free( KSet );
}

/*------------------------------------------------------------------------------*/

int VICisintvec( int n, double* pi ) {
/* Checks if vector pi is integer */

  double api;
  int   j;
  for ( j=0; j < n; j++ ) {
    api = fabs( pi[j] );
    if ( (api - trunc( api + VICtol ) ) > VICtol ) return( 0 );
  }
  return ( 1 );
}

/*------------------------------------------------------------------------------*/


static
void VICintcoeff( int n, int maxmult, double* pi, double* rhs ) {
/*
  Decscription: Transforms coefficients pi of the inequality pi*x <= rhs to
  integers following a method suggested by Kaparis and Letchford.
*/

 int    j, mult=1, ok;
 double api, base[n+1], frac, minpi=CPX_INFBOUND;

 for ( j=0; j < n ; j++ ) {
   api = fabs( pi[j] );
   if ( (0 < api) && (api < minpi) ) minpi = api;
 }
 *rhs /= minpi;
 base[n]=*rhs;
 for ( j=0; j < n; j++ ) pi[j] /= minpi, base[j] = pi[j];
 ok = VICisintvec( n, pi );

 while ( (mult < maxmult) && (!(ok) ) ) {
   *rhs = (++mult)*base[n];
   for ( j=0; j < n; j++ ) pi[j] = base[j]*mult;
   ok = VICisintvec( n, pi );
 }

 /* Above heuristic failed. Return a weakened inequality obtained
    by simply rounding down the left-hand and right-hand sides.*/
 for ( j=0; j < n; j++ ) {
   if ( pi[j] > 0.0 )
     pi[j] = trunc( pi[j] + VICtol );
   else if ( pi[j] < 0.0 ) {
     api  = trunc( -pi[j] + VICtol );
     frac = abs(-pi[j] - api);
     pi[j]= -api;
     if ( frac > VICtol ) pi[j] -= 1.0;
   }
 }
 if ( *rhs > 0 )
   *rhs = trunc( *rhs + VICtol );
 else if ( *rhs < 0.0 ) {
   api = trunc( -(*rhs) + VICtol );
   frac= abs( -(*rhs) - api );
   *rhs = -api;
   if ( frac > VICtol ) *rhs -= 1.0;
 }

}

/*------------------------------------------------------------------------------*/

//void VICsnffenchel( CPXENVptr Env, int n, int demand, int* capaci, double* x,
//                    double* y, int* indx, int* indy, VICcut** cut ) {
///* Description:
//   --------------
//   This procedure tries to find a fenchel cutting plane that cuts of a fractional
//   solution (x*,z*) to the single node flow problem
//
//        \sum_{j=1}^n x(j) = demand,
//        0 <= x(j) <= z(j) for j=1,...,n,
//        z(j)=0 or z(j)=capaci(j) for j=1,...,n.
//
//   Let X be the set of all feasible solutions (x,z) to the single-node flow
//   problem above. A fenchel cut is then given by
//
//         \sum_{j=1}^n \pi(j) x(j) - \sum_{j=1} \lambda(j) z(j) <= \pi_0,
//
//   where (\pi, \lambda, \pi_0) is a solution to the separation problem
//          \max  \pi x* - \lambda z*
//          s.t.: \pi x  - \lambda y <= \pi_0 \forall (x,y)\in X
//                0 <= \pi(j) <= 1.0 forall j,
//                0 <= \lambda(j) <= 1.0 forall j.
//
//   The separation problem is solved by means of subgradient optimization.
//
//   Scope: Export
//   -------------
//
//   Parameters:
//   -------------
//   - n      : number of arcs from the n sinks to the single source
//   - demand : the sink's demand
//   - capaci : pointer to an array of integers of size of at least n containing
//              the arcs' capacities
//   - x      : pointer to an array of doubles of size of at leat n containing
//              the values of flow variables in the fractional solution
//   - y      : pointer to an array of doubles of size of at least n containing
//              the values to the binary variables in the fractional solution,
//	      0 <= y[j] <= 1 (the variable z(j) is defined as z(j)=capaci(j)y(j)),
//   - indx   : NULL or a pointer to an integer array of size of at least n containing
//              the indices of the x-variables. If indx is NULL, the indices 0,..,n-1
//              are used.
//   - indy   : NULL or a pointer to an integer array of size of at least n containing
//              the indices of the y-variables. If indy is NULL, the indices n,..., 2n-1
//              are used.
//   - cut    : pointer to a pointer to a cut structure VICcut defined in
//              vico.h. If no violated fenchel cut is found, the null
//              pointer is returned; otherwise **cut contains the cut. Coefficient
//              of the cut are returned in terms of variables (x,y)
//*/
//  int     fail=0, iter, j, jj, n_arcs=n, nn, numnz, zj;
//  double  lowbnd=-VICinf, cur_low, cur_up, cur_up2, upbnd=1.0E30;
//  double  fvalue, dnorm, sigma=FCHELopt.alpha, theta, *tmp, zqj;
//  double  *dwork, *pi, *lambda, *sg0, *sg1, *sg, *sgx, *sgz, *tcost,
//          *fcost, *direc, *zq, *xq;
//  int     *order = (int*) calloc( n, sizeof( int ) );
//  int     *iwork, *flow, *isfree, *freex, *freey, *cap;
//  VICcut  *fenchel=NULL;
//  FILE    *fchel_log=NULL;
//
//  /* Arcs are reindexed such that non-zero arcs are listed first */
//  if ( FCHELopt.Reduce ) {
//    for ( j=n=0; j < n_arcs; j++ ) if ( y[j] >= VICtol ) order[n++] = j;
//  }
//  else for ( j=n=0; j < n_arcs; j++ ) if ( capaci[j] > 0) order[n++]=j;
//  nn = 2*n;
//
//  /* Allocate working space memory */
//  iwork = (int*) calloc( 4*n, sizeof( int ) );
//  dwork = (double*) calloc( 14*n, sizeof(double) );
//  pi = dwork, lambda=pi+n, sg0=lambda+n , sg1=sg0+nn,    sg=sg1+nn,   sgx=sg,
//  sgz=sgx+n , tcost=sg+nn, fcost=tcost+n, direc=fcost+n, zq=direc+nn, xq=zq+n;
//  flow = iwork, isfree=flow+n, freex=isfree, freey=freex+n, cap=freey+n;
//
//  /* Open log file */
//  if ( FCHELopt.CHK ){
//    fchel_log = fopen( "fenchel.log", "a" );
//    if (fchel_log ) {
//      fprintf(fchel_log,"\n*GENERATING FENCHEL CUTTING PLANE*\n");
//      fprintf(fchel_log,
//        "SG_ITER  LOWBOUND  UPBOUND  CURR_LOW  CURR_UP  STEPPARA  STEPLEN\n");
//    }
//  }
//
//  /* Check CPLEX environment if CPLEX has to be used for subproblem */
//  if ( ( Env==NULL ) && (FCHELopt.Algo==ALG_CPLEX ) ) FCHELopt.Algo=ALG_MT1;
//
//  /* initial solution and upper bound. Note that calloc intializes with zeros */
//  for ( j=0, upbnd=0.0; j < n; j++ ) {
//    cap[j] = capaci[jj = order[j]];
//    freex[j] = (int) ( x[jj] > VICtol );
//    if ( freex[j] ) {
//      upbnd   += cap[j]*(1.0-y[jj]);
//      freex[j] = ( x[jj] + VICtol < capaci[jj] );
//      pi[j]    = 1-freex[j];
//    }
//    freey[j] = (int) ( y[jj] < VICone );
//    if ( y[jj] < VICtol ) freey[j]=0, lambda[j]=1.0;
//  }
//
//  /* subgradient steps */
//  for ( iter=0; iter < FCHELopt.maxit; iter++ ) {
//
//    /* Solve the subproblem max{\pi x - \lambda z : (x,z)\in X } for the
//       current (\pi,\lambda). Note that ssfctp_mt1 solves a minization
//       problem in variables (x,y). The subproblem is thus transformed to
//       demand - min (1-\pi(j))x + (\lambda(j)*capaci(j))y
//                        s.t.: \sum_j x(j) = demand,
//                              0 <= x(j) <= capaci(j)*y(j) forall j,
//                              y(j)=0,1 forall j.
//    */
//    for ( j=0; j < n; j++ ) tcost[j]=1.0-pi[j], fcost[j]=cap[j]*lambda[j];
//    switch (FCHELopt.Algo ) {
//      case ALG_MT1:
//        ssfctp_mt1( n, demand, cap, fcost, tcost, 1, flow, &fvalue );
//	break;
//      case ALG_DP:
//        ssfctp_dp( n, demand, cap, fcost, tcost, 2, flow, &fvalue );
//	break;
//      case ALG_CPLEX:
//        ssfctp_cpx( Env, n, demand, cap, fcost, tcost, flow, &fvalue );
//	break;
//    }
//    printf(".");
//
//    fvalue = demand - fvalue;
//
//    /* Get the current lower bound, upper bound and the subgradient */
//    for ( j=0, cur_low=cur_up=cur_up2=0.0; j < n; j++ ) {
//      zj       = flow[j] > 0 ? cap[j] : 0.0;
//      zqj      = cap[j]*y[jj=order[j]];
//      sgx[j]   = x[jj]-flow[j];
//      sgz[j]   = zj-zqj;
//      cur_up  += MAX(0, sgx[j] ) + MAX(0,sgz[j]);
//      cur_low += pi[j]*sgx[j] + lambda[j]*sgz[j];
//      zq[j]    = (iter*zq[j] + (double)zj )/(double)(iter+1);
//      xq[j]    = (iter*xq[j] + (double)flow[j])/(double)(iter+1);
//      cur_up2 += MAX(0,x[jj]-xq[j]) + MAX(0,zq[j]-zqj);
//    }
//    if ( cur_up2 < cur_up ) cur_up = cur_up2;
//
//    /* Check for improvement in lower and upper bound */
//    if ( cur_up < upbnd ) upbnd = cur_up;
//    if ( cur_low > lowbnd + VICzero ) {
//      lowbnd = cur_low;
//      if ( lowbnd > VICtol ) { /* inequality is violated */
//        if ( fenchel==NULL ) if ( !( VICallocut( nn, &fenchel ) ) ) break;
//        fenchel->sense = 'L';
//        fenchel->rhs = fvalue;
//        for ( j=numnz=0; j < n; j++ ) {
//          jj = order[j];
//	  if ( pi[j] > VICzero ) {
//            fenchel->nzval[numnz] = pi[j];
//            fenchel->nzind[numnz++] = (indx) ? indx[jj] : jj;
//	  }
//	  if ( lambda[j] > VICzero ) {
//            fenchel->nzval[numnz] = -lambda[j]*cap[j];
//            fenchel->nzind[numnz++] = (indy) ? indy[jj] : n_arcs+jj;
//	  }
//        }
//	fenchel->nzcnt = numnz;
//      }
//    } else if ( (FCHELopt.H > 0) && (++fail > FCHELopt.H ) ) {
//      fail = 0;
//      sigma /= 2.0;
//      if ( sigma < VICzero ) break;
//    }
//
//    if ( lowbnd - 1.0 > upbnd ) {
//      printf("\n ERROR in FENCHEL CUT GENERATOR.\n");
//      break;
//    }
//
//    if ( ( upbnd - lowbnd )/MAX(1.0, lowbnd ) < VICtol ) break;
//    if ( upbnd < 0.1 ) break; /* continuing seems no more meaningful */
//
//    if ( iter < 1 ) for ( j=0; j < nn; j++ ) direc[j]=sg1[j]=sg[j];
//    if ( iter < 2 ) for ( j=0; j < nn; j++ ) sg0[j] = sg1[j];
//
//    /* Compute direction according to the specified subgradient
//       strategy */
//    switch ( FCHELopt.sg_strat ) {
//      case SG_DEFAULT:
//        for ( j=0, dnorm=0.0; j < nn; j++ ) {
//          direc[j] = sg[j];
//          dnorm   += direc[j]*direc[j];
//        }
//        break;
//      case SG_DEFLECT:
//        for ( j=0, dnorm=0.0; j < nn; j++ ) {
//          direc[j] = (sg[j] + 0.3*sg1[j] + 0.1*sg0[j])/1.4;
//          dnorm   += direc[j]*direc[j];
//        }
//        break;
//      case SG_SMOOTH:
//        for ( j=0, dnorm=0.0; j < nn; j++ ) {
//          direc[j]+= 0.3*(sg[j]-direc[j]);
//          dnorm   += sg[j]*sg[j];
//	}
//    }
//
//    theta = sigma*(upbnd-cur_low)/dnorm;
//    for ( j=0; j < nn; j++ ) if ( isfree[j] ) {
//      pi[j] += theta*direc[j];
//      if ( pi[j] < 0.0 ) pi[j]=0.0; else if ( pi[j] > 1.0 ) pi[j]=1.0;
//    }
//    tmp=sg0, sg0=sg1, sg1=sg, sg=tmp, sgx=sg, sgz=sgx+n;
//
//    if ( fchel_log ) {
//      fprintf(fchel_log,"%7d  %8.4f  %7.4f  %8.4f  %7.4f  %8.6f  %7.6f\n",
//        iter,lowbnd,upbnd,cur_low,cur_up,sigma,theta);
//    }
//
//    if ( FCHELopt.H==0 ) sigma /= 1.05;
//
//  } /*End of subgradient steps */
//
//  free( dwork );
//  free( iwork );
//  free( order );
//  *cut = fenchel;
//  if ( fchel_log ) fclose( fchel_log );
//  printf("\n");
//
//}
/*------------------------------------------------------------------------------*/

static
void VICinitsepmaster( CPXENVptr Env, CPXLPptr* master, int n, int* order,
                       double* xlp ) {
/*
   Description: Initialises the master problem used for exact knapsack separation
                with columns x^j = e_j (j-th unit vector )
*/

  int      j, status;
  double   *dwork = (double*) calloc( 5*n, sizeof( double ) );
  int      *iwork = ( int* ) calloc( 3*n, sizeof( int ) );
  char     *sense = (char*) calloc( n, sizeof( char ) );
  char     *probname = (char*) calloc(16,sizeof(char) );

  double   *rhs=dwork, *matval=rhs+n, *obj=matval+n, *lb=obj+n, *ub=lb+n;
  int      *matbeg=iwork, *matcnt=matbeg+n, *matind=matcnt+n;


  for ( j=0; j < n; j++ ) {
    rhs[j]=xlp[order[j]];
    sense[j] = 'G';
    lb[j] = 0.0;
    ub[j] = CPX_INFBOUND;
    matval[j] = 1.0;
    matcnt[j] = 1;
    matind[j] = j;
    matbeg[j] = j;
    obj[j] = 1.0;
  }

  strcpy( probname, "SEPMASTER");
  *master = CPXcreateprob( Env, &status, probname );
  if ( *master ) {
    status = CPXcopylp( Env, *master, n, n, CPX_MIN, obj, rhs, sense, matbeg, matcnt,
                        matind, matval, lb, ub, NULL );
    if ( status > 0 ) CPXfreeprob( Env, master );
  }
  free( dwork );
  free( iwork );
  free( sense );
  free( probname );
}

/*------------------------------------------------------------------------------*/

static
char VICkpsepsub( int n, int cap, int* order, int* weight, double* pi,
                  item* items, double* rhs ) {
/*
   Description: Checks if the inequality pi*x <= rhs is a valid inequality for the
                knapsack problem {weight x <= cap, x binary} and returns 0 if this
                is the case. The check is performed by maximizing pi*x over the
                knapsack polytope using the combo algorithm for the knapsack
                problem. On output *rhs is replaced with the maximum value of the
                inequalities left-hand side
*/
  int    icoeff, iobj, isint, isviol, j, last=n, nn=0, ub=12000, wsum, ww;
  item   *cur;
  double lhs, scale=10000;

  isint = VICisintvec( n, pi );
  if ( isint ) scale=1.0, ub = (int) trunc( *rhs + VICtol )+1;

  for ( j=wsum=0; j < n; j++ ) {
    ww = weight[order[j]];
    icoeff = (int) round( pi[j]*scale );
    if ( icoeff > 0 )
      wsum += ww, cur=items+nn, nn++, cur->x=1;
    else
      last--, cur = items + last, cur->x=0;
    cur->num = j;
    cur->w = ww;
    cur->p = icoeff;
  }
  if ( wsum > cap ) {
    if ( nn > 1 )
      iobj = combo( items, items+nn-1, cap, 0, ub, 1, 0 );
    else
      items->x=0;
  }
  for ( j=0,lhs=0.0; j < nn; j++ ) if ( items[j].x ) lhs += pi[items[j].num];
  isviol = ( lhs > *rhs + VICtol );
  *rhs   = lhs;
  return ( isviol );

}

/*------------------------------------------------------------------------------*/

static
int VICkpsepaddcol( CPXENVptr Env, CPXLPptr master, int n, item *items ) {
/*
   Description: Add a new column specified by the solution to the knapsack problem
   to the master problem that is used in the exact knapsack separation procedure
*/

  int    j, nzcnt, cmatbeg=0, status=0;
  int    cmatind[n];
  double cmatval[n], obj=1.0, lb=0.0, ub=CPX_INFBOUND;

  for ( j=nzcnt=0; j < n; j++ ) if ( items[j].x ) {
    cmatval[nzcnt]   = 1.0;
    cmatind[nzcnt++] = items[j].num;
  }
  status = CPXaddcols ( Env, master, 1, nzcnt, &obj, &cmatbeg, cmatind, cmatval,
                        &lb, &ub, NULL );
  return ( status );
}

/*------------------------------------------------------------------------------*/

static
int VICkplift ( int justUp, int n, int nn, int rcap, int* weight, int* order,
                int* coeff, double* rco, double* xlp ) {
/*
  Description:
  ------------
  Given the inequality

  coeff[order[0]]*x[order[0]] + ... + coeff[order[nn-1]]*x[order[nn-1]] <= coeff[n]

  that is valid for the knapsack problem if x[order[i]] is fixed to xlp[order[i]=0,1
  for i=nn, ..., n-1, the routine first down-lifts (if justUp=0) variables x[order[i]] that
  are fixed to one in order to obtain a valid inequality and then up-lifts variables
  x[order[i]] fixed to zero to further sharpen the inequality. Within these two
  groups of variables, variables are lifted according to increasing magnitudes of
  reduced costs. In case of a tie, the variable showing smaller weight in the
  knapsack inequality is lifted first.

  A dynamic programming method is used to perform the lifting. The approach exploits
  the recursive function

  omega_i(t) = min\{ omega_{i-1}(t), weight[order[i-1]] + omega_{i-1}(t-coeff[order[i-1]])\}

  for computing the lifted coefficient coeff[order[i]]. The function must be evaluated for
  t=0, ..., \sum_{j=0}^{i-1} coeff[order[j]]. The method may thus get into severe
  troubles if the coefficients in the inequality are too large!

  The return value is zero if no error could be detected; otherwise 1
  is returned. rcap must equal the residual capacity of the knapsack.
*/

  int    cnt, i, j=0, k, nn1, t, T, Tmax, Tprev, osize = coeff[n]+1, pi, w;
  int    *omega=NULL, *tmparr = NULL;
  double rcmin;

  if ( justUp == 0 ) {
    for ( j=T=0, osize=5; j < nn; j++ ) {
      T += (w=coeff[order[j]]);
      if ( w > osize ) osize =w ;
    }
    osize *= 2*n;
  }
  omega  = (int*) calloc( osize, sizeof( int ) );
  if ( omega==NULL ) return( 1 );

  /* Initialize the function omega on stage i=1 */
  for ( t=1, i=order[0], T=coeff[i], w=weight[i]; t <= T; t++ ) omega[t] = w;

  /* Recursively compute the function omega for stages i=1, ..., nn */
  Tmax = osize-1;
  for ( i=2; i <= nn; i++ ) {
    j=order[i-1], pi = coeff[j], w=weight[j], Tprev = T, T += pi;
    if ( T > Tmax ) T = Tmax;
    for ( t=Tprev+1; t <= pi; t++ ) omega[t] = w;
    for ( t=MAX(Tprev,pi)+1; t <= T; t++ ) omega[t] = w + omega[t-pi];
    for ( t=Tprev; t > pi; t-- ) omega[t]= MIN( omega[t], omega[t-pi] + w );
    for ( t=pi; t > 0; t-- ) omega[t] = MIN( omega[t], w );
  }

  /* Rearrange variables in order[nn] ... order[n-1] such that variables to be
     down-lifted are listed first. It is assumed that just variables showing an
     LP value of almost 1 were fixed and have thus to be down-lifted  */
  nn1 = nn;
  if ( justUp == 0 ) {
    j=nn, nn1=n;
    while ( j < nn1 ) {
      if ( xlp[order[j]] < VICone ) VICswap( order+j, order+(--nn1) ); else j++;
    }
  }

  /* Down-lift variables that were fixed to one */
  for ( cnt=nn; cnt <  nn1; cnt++ ) {
    /* Determine variable x[order[i]] to be lifted next */
    for ( i=cnt, j=cnt+1, w=weight[k=order[i]], rcmin=rco[k]; j < nn1; j++ )
      if ( (rco[k=order[j]] < rcmin) || ( (rco[k] < rcmin + VICtol) && (weight[k] < w) ) )
        w = weight[k], rcmin = rco[k], i=j;
    VICswap ( order + cnt, order + i );
    /* Determine the lifted coefficient */
    for ( t=0, Tprev=T, rcap+=w; t <= Tprev; t++ ) if ( omega[t] > rcap ) break;
    coeff[i=order[cnt]] = pi = (--t) - coeff[n];
    coeff[n] += coeff[i];
    T += coeff[i];
    if ( cnt < n-1 ) {
      /* Check if more space is required to store function values omega */
      if ( T >= osize ) {
        osize = MAX( 2*T+1, osize + 100 );
        omega = (int*) realloc ( omega, osize*sizeof(int) );
        if ( omega==NULL ) return( 1 );
      }
      /* Compute function omega for the next stage */
      for ( t=Tprev+1; t <= pi; t++ ) omega[t] = w;
      for ( t=MAX(Tprev,pi)+1; t <= T; t++ ) omega[t] = w + omega[t-pi];
      for ( t=Tprev; t > pi; t-- ) omega[t] = MIN( omega[t], omega[t-pi] + w );
      for ( t=pi; t > 0; t-- ) omega[t] = MIN( omega[t], w );
    }
  }

  // The following is old code that was not really working and is outcommented
  // Up-lift variables that were fixed to zero
  //for ( cnt=nn1; cnt <  n; cnt++ ) {
    // Determine variable x[order[i]] to be lifted next */
  //  for ( i=cnt, j=cnt+1, w=weight[k=order[i]], rcmin=rco[k]; j < n; j++ )
  //    if ( (rco[k=order[j]] < rcmin) || ( (rco[k] < rcmin + VICtol) && (weight[k] < w) ) )
  //     w = weight[k], rcmin = rco[k], i=j;
  //  VICswap ( order + cnt, order + i );
    /* Determine the lifted coefficient */
  //  for ( t=0, Tprev=T, rcap-=w; t <= Tprev; t++ ) if ( omega[t] > rcap ) break;
  //  coeff[i=order[cnt]] = pi = coeff[n] - (--t);
  //  T += coeff[i];
  //  rcap += w;
    /* Check if more space is required to store function values omega */
  //  if ( cnt < n-1 ) {
  //    if ( T+1 > osize ) {
  //      osize += MAX(T+1, osize + 100 );
  //      omega = (int*) realloc ( omega, osize );
  //      if ( omega==NULL ) return( 1 );
  //    }
      /* Compute function omega for the next stage */
  //    for ( t=Tprev+1; t <= pi; t++ ) omega[t] = w;
  //    for ( t=MAX(Tprev,pi)+1; t <= T; t++ ) omega[t] = w + omega[t-pi];
  //    for ( t=Tprev; t > pi; t-- ) omega[t] = MIN( omega[t], omega[t-pi] + w );
  //    for ( t=pi; t > 0; t-- ) omega[t] = MIN( omega[t], w );
  //  }
  //}


  //  Now uplift the valid inequality

  // -> Eliminate first variables that need to be zero because weight[j] > capaci
  j=nn1, nn=n;
  while ( j < nn ) {
    if ( weight[order[j]] > rcap ) VICswap( order+j, order + (--nn) ); else j++;
  }

  // Check if enough memory is allocated to omega array
  if ( osize <= coeff[n] ) {
    tmparr = (int*) calloc( coeff[n], sizeof( int ) );
    memcpy(tmparr, omega, T*sizeof(int) );
    free ( omega );
    omega = tmparr;
  }

  // Uplift the variables in order[nn1 ... nn-1 ]. Note that for doing the uplifting
  // we only need to know at most the omega values in the range from 0 ... coeff[n].
  T = MIN( T, coeff[n] );
  for ( cnt=nn1; cnt <  nn; cnt++ ) {
    // Determine variable x[order[i]] to be lifted next */
    for ( i=cnt, j=cnt+1, w=weight[k=order[i]], rcmin=rco[k]; j < nn; j++ )
      if ( (rco[k=order[j]] < rcmin) || ( (rco[k] < rcmin + VICtol) && (weight[k] < w) ) )
        w = weight[k], rcmin = rco[k], i=j;
    VICswap ( order + cnt, order + i );
    /* Determine the lifted coefficient */
    for ( t=0, pi = 0, Tprev=T, rcap-=w; t <= Tprev; t++ ) {
      if ( t > coeff[n] ) break; // lifted coefficient is zero
      if ( omega[t] > rcap ) {
        pi = coeff[n] - (--t); // t = max_{t'} \{ t' : \omega(t') <= rcap \}
	break;
      }
    }
    coeff[i=order[cnt]] = pi;
    rcap += w;
    T = MIN( T+pi, coeff[n] );
    if ( cnt < nn-1 ) {
    /* Compute function omega in range 0 ... coeff[n] (at most) for the next stage. */
      for ( t=Tprev+1; t <= pi; t++ ) omega[t] = w;
      for ( t=MAX(Tprev,pi)+1; t <= T; t++ ) omega[t] = w + omega[t-pi];
      for ( t=Tprev; t > pi; t-- ) omega[t] = MIN( omega[t], omega[t-pi] + w );
      for ( t=pi; t > 0; t-- ) omega[t] = MIN( omega[t], w );
    }
  }
  // Variables in order[nn]...order[n-1] need to be zero and obtain lifting
  // coefficient rcap+1
  for ( j=nn; j < n; j++ ) coeff[order[j]] = coeff[n]+1;

  free( omega );
  return( 0 );
}

/*------------------------------------------------------------------------------*/

void VICkpsep( CPXENVptr Env, int justUp, int n, int cap, char sense, int* weight,
               int* indx, double* xlp, double* rco, VICcut** cut ) {
/*
  Description:
  --------------
  Exact knapsack separation. Given a fractional solution x*, the most violated valid
  inequality \pi x \le 1 is returned by optimizing over the 1-polar, that is by solving

     max\{ \pi x* : \pi x^t \le 1 for each feasible solution x^t \}

  The dual of the above separation problem is solved by means of column generation.
  To reduce the effort, the separation is performed for a reduced polytope obtained
  by resp. fixing variables x_j to 0 and 1 if x*_j=0 or x^*_j=1. Afterwards coefficients
  of fixed variables are lifted in order to obtain a valid inequality.

  Scope: Export
  --------------

  Parameters:
  -------------
   - Env     : pointer to the CPLEX environment. Env must be open!
   - justUp  : if equal to one, only variables showing a value of zero in the
               LP solution are first fixed to zero. The resulting inequality
	       is then just uplifted to obtain a (strengthened) inequality
	       for the full polytope. Otherwise (justUp=0) both variables
	       showing LP-value of zero and one are fixed. The resulting
	       inequality is then first to be downlifted to obtain a valid
	       inequality, and thereafter uplifting is applied.
   - n       : number of variables in the knapsack inequality
   - cap     : right hand side of the knapsack inequality
   - sense   : sense (direction) of the knapsack inequality. In case of
               sense='L' the inequality is given by
               weight[0]*x[indx[0]] + .... + weight[n-1]*x[indx[n-1]] <= cap
               Otherwise it is assumed that the inequality is
               weight[0]*x[indx[0]] + .... + weight[n-1]*x[indx[n-1]] >= cap
   - weight  : pointer to an array of integers of size of at least n
               containing the coefficients of variables/items in the knapsack
               inequality. The coefficients are not restricted to be
               nonnegative!
   - indx    : pointer to an array of integers of size of at least n
               containing the indices of the variables arising in the
               knapsack inequality.
               If indx is NULL, it is assumed that these variables are
               numbered from 0 to n-1
   - xlp     : pointer to an array of doubles of size of at least n
               containing the LP solution.
   - rco     : pointer to an array of doubles of size of at least n
               containing the absolute values of reduced costs. Reduced
               cost information is used for deciding on the lifting sequence.
   - cut     : pointer to a pointer to a cut structure VICcut defined in
               vico.h. If no violated inequality is found, the null
               pointer is returned; otherwise **cut contains the cut.
 */

  char   *invers = calloc(n, sizeof(char) );
  int    *order = (int*) calloc( n, sizeof( int ) );
  int    gcd, i, j, last, nn, rcap, status, w, wsum, *coeff=NULL;
  double mobjval=0.0, minpi, rhs=1.0, *pi=NULL, myOne=VICone;
  item   *items=NULL;

  CPXLPptr master=NULL;
  if ( justUp ) myOne = 2.0;

  /* Transform to standard binary knapsack structure */
  *cut=NULL;
  rcap=cap;
  VICkidefault( n, sense, invers, &rcap, weight, xlp );

  /* Fix variables and remove variables of weight larger than residual
     knapsack capacity */
  for ( j=nn=0, last=n; j < n; j++ )
    if ( xlp[j] < VICtol )
      order[--last] = j;
    else if ( xlp[j] > myOne ) {
      order[--last] = j;
      rcap -= weight[j];
    }
    else
      order[nn++] = j;
  for ( j=0; j < nn; j++ ) if ( weight[i=order[j]] > rcap ) {
    order[j] = order[--nn];
    order[nn]= i;
  }
  if ( nn==1 ) goto RETURN;

  /* Check if the fractional solution must trivially lie in the
     reduced knapsack polytope */
  for ( j=w=0; (j < nn) && (w==0); j++ ) w = weight[order[j]];
  if ( j == nn ) goto RETURN;
  for ( i=j; i < nn; i++ ) if ( weight[order[i]] != w ) break;
  if ( i==nn ) goto RETURN;
  for ( j=wsum=0; j < nn; j++ ) if ( xlp[i=order[j]] > VICtol ) wsum += weight[i];
  if ( wsum <= rcap ) goto RETURN;

  pi = ( double* ) calloc( nn, sizeof( double ) );
  items = ( item* ) calloc( nn, sizeof( item ) ); /* items of knapsack problem */

  /* Initialise master separation problem */
  VICinitsepmaster( Env, &master, nn, order, xlp );
  if ( master == NULL ) goto RETURN; /* shit happens */

  /* Solve initial master problem */
  status = CPXoptimize( Env, master );
  if ( status == 0 ) status = CPXgetobjval( Env, master, &mobjval );
  if ( status > 0 ) mobjval = 0.0; /* some error occurred */

  /* Apply column generation to solve the complete master problem */
  while ( mobjval > 1.0+VICtol ) {
    /* Retrieve right-hand side coefficients of current inequality, that is
       the dual variables of the current master program */
    if ( ( status = CPXgetpi( Env, master, pi, 0, nn-1 ) ) > 0 ) break;

    /* If the inequality pi*x <= 1 is valid, we are done. Check this
       by solving the knapsack problem */
    rhs=1.0;
    if ( !( VICkpsepsub( nn, rcap, order, weight, pi, items, &rhs ) ) ) break;

    /* Current inequality is invalid. The solution to the last knapsack problem
       is thus added as a new column to the master problem. */
    if ( ( status = VICkpsepaddcol( Env, master, nn, items ) ) > 0 ) break;

    /* Reoptimize the master problem */
    if ( ( status = CPXprimopt( Env, master ) ) > 0 ) break;
    if ( ( status = CPXgetobjval( Env, master, &mobjval ) ) > 0 ) break;

  };
  if ( master ) status = CPXfreeprob( Env, &master );
  if ( (status > 0) || ( mobjval <= 1.0+VICtol) ) goto RETURN;

  /* A violated facet-defining inequality has been found if Status==0 and
     mobjval > 1.0 + VICtol. Make the inequality's coefficients integer */
  rhs = 1.0;
  if ( !( VICisintvec( nn, pi ) ) ) {
    for ( j=0, minpi=CPX_INFBOUND; j < nn; j++ )
      if ( (0 < pi[j]) && ( pi[j] < minpi ) ) minpi = pi[j];
    VICintcoeff( nn, minpi*rcap+1, pi, &rhs );
    /* Check validity of resulting integer inequality */
    VICkpsepsub( nn, rcap, order, weight, pi, items, &rhs );
  }

  coeff = ( int* ) calloc(n+1, sizeof(int) );
  coeff[n] = (int) round( rhs+VICtol ); /* rhs of current inequality */
  for ( j=0, mobjval=0.0; j < nn; j++ ) {
    coeff[i=order[j]] = (int) round( pi[j] + VICtol );
    mobjval += coeff[i]*xlp[i];
  }
  if ( mobjval < (double)coeff[n] + VICtol ) goto RETURN;

  /* Divide the inequality obtained so far by greatest common divisor */
  for ( j=0, gcd=coeff[n]; (j < nn) && (gcd > 1); j++ ) {
    if ( (w =coeff[order[j]]) > 0 ) {
      if ( gcd > w ) VICswap( &gcd, &w );
      do {
        w = w % gcd;
        if ( w > 0 ) VICswap( &gcd, &w );
      } while ( w > 0 );
    }
  }
  coeff[n] /= gcd;
  for ( j=0; j < nn; j++ ) coeff[order[j]] /= gcd;

  /* Next the inequality obtained so far has to be lifted in order
     to obtain a valid inequality */
  if ( ( status = VICkplift(justUp, n, nn, rcap, weight, order, coeff, rco, xlp) ) > 0 ) goto RETURN;

  /* Setting a number of variables to 1 might have implied that some variables had to
     be zero although their LP value is larger than zero. In the lifting this variables
     were uplifted. It might nevertheless be that the lifted inequality is not violated.
     Hence, check if the lifted inequality is still violated. */
  for ( j=0, mobjval=0.0; j < n; j++ ) if ( coeff[j] > 0 ) mobjval += xlp[j]*coeff[j];
  if ( mobjval < (double)coeff[n] + VICtol ) goto RETURN;

  /* Return the cut in appropriate data structure */
  for ( j=nn=0; j < n; j++ ) nn += ( coeff[j] > 0 );
  if ( VICallocut( nn, cut ) ) {
    (*cut)->sense = 'L';
    (*cut)->rhs   = (double)coeff[n];
    (*cut)->nzcnt = nn;
    for ( j=nn=0; j < n; j++ ) if ( coeff[j] > 0 ) {
      (*cut)->nzind[nn] = (indx) ? indx[j] : j;
      if ( invers[j] ) {
        (*cut)->rhs -= (double)coeff[j];
        (*cut)->nzval[nn++] = -(double)coeff[j];
      }
      else
        (*cut)->nzval[nn++] = (double)coeff[j];
    }
  }

RETURN:
 /* Free allocated memory */
  VICkirestore( n, sense, invers, &rcap, weight, xlp );
  free ( invers );
  free ( order );
  if ( pi ) free( pi );
  if ( items ) free ( items );
  if ( coeff ) free ( coeff );
}

