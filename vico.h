/*!
 * \brief Implements seperations routines several different problems.
 *
 * This is the header file for module "vico.c" (Valid Inequalities for selected Combinatorial Optimization problems)
 *
 * \author Andreas Klose
 * \version 2.9.0
 * \note    Programming language: C
 * \note    In order to use the functions listed below from within a SUN Pascal program compiled with option -L using the SUN Pascal compiler,
 *          compile file vico.c with option -DSUNPAS and link with libF77
 *
 * \section License
 *
 * Copyright 2016, Andreas Klose.
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
 *
 * If you use the software in any academic work, please make a reference to "An LP-Based Heuristic for Two-Stage Capacitated Facility Location Problems",
 * Journal of the Operational Research Society, Vol. 50, No. 2 (Feb. 1999), pp. 157-166.
 *
 *
 * \latexonly
 * \section{Change log for vico.h}
 * \begin{center}
 *     \begin{tabularx}{\textwidth}{llr X}\toprule
 *        FILE:          &   \multicolumn{3}{l}{vico.h and vico.c}\\
 *        Version:       &   \multicolumn{3}{l}{2.9.0}\\
 *        \multicolumn{4}{l}{- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -}\\
 *        CHANGE LOG:    &   DATE    &   VER.-NO.    &   CHANGES MADE\\ \midrule
 *                       &   2003--05--12   &   1.0.0       & first implementation\\
 *                       &   2003-07-07     &   1.1.0       & bug in VICuflohi removed\\
 *	                     &   2003-08-06     &   2.0.0       & routine VIClci and VICkconf added\\
 *	                     &   2003-08-07     &   2.1.0       & routine VICgapfn added\\
 *	                     &   2003-08-15     &   2.2.0       & routine VICuflcov added\\
 *	                     &   2003-08-22     &   2.3.0       & VIClci, VICkconf modified to ease handling\\
 *	                     &   2003-08-23     &   2.3.1       & bugs in VICuflohi removed (shortest path computation now based on FIFO, error in cut existence check removed)\\
 *	                     &   2003-08-23     &   2.3.2       & minor changes, some "statics" inserted\\
 *	                     &   2003-08-27     &   2.3.3       & minor changes: routines VICsearchcut, VICclearlst added\\
 *	                     &   2004-02-20     &   2.4.0	    & routine VICeci added\\
 *	                     &   2004-02-24     &   2.4.1       & routine VICintsort, VICdblsort added\\
 *	                     &   2004-02-25     &   2.5.0       & routine VICkpreduce added\\
 *                       &   2004-03-11     &   2.5.1       & small bug in VIClci removed\\
 *                       &   2005-05-10     &   2.5.2       & bug in VICkconf "(if card==n ) ... return" removed\\
 *                       &   2005-08-01     &   2.5.3       & bug in kirestore() removed (all coefficients of a generated cut were modified if at least one variable was inverted!)\\
 *                       &   2007-05-30     &   2.6.0       & Started with implementing Fenchel cuts based on single-node flow structures\\
 *                       &   2007-06-05     &   2.6.1       & Different subgradient strategies for Fenchel cut generation implemented.\\
 *                       &   2007-06-14     &   2.6.2       & Different algorithms for solving the subproblem within Fenchel cut generation may optionally be chosen\\
 *                       &   2007-06-15     &   2.6.3       & Additional parameter "Freq" for Fenchel cut generation introduced\\
 *                       &   2007-06-22     &   2.7.0       & Additional parameter "Reduce" for Fenchel cut generation introduced. If Reduce=0 all "zero arcs" are removed and the inequality is generated for the reduced polytope\\
 *                       &   2007-08-20     &   2.7.1       & small bug in vicsnffenchel removed: capacities are now allowed to be zero (variable can then be ignored)\\
 *                       &   2008-01-09     &   2.8.0       & start to implement routine for exact knapsack separation\\
 *	                     &   2008-06-18     &   2.8.1       & Some bugs removed in VICkplift and VICkpsep\\
 *	                     &   2008-06-19     &   2.8.2       & Bug remove in VICkplift ( (t-pi) could be negative)\\
 *	                     &   2012-08-20     &   2.8.3       & Adjusted the uplifting in VICkplift and included possibility to exclusively fix variables of zero LP value when defining the reduced knapsack polytope in the exact knapsack separation procedure\\
 *                       &   2012-12-18     &   2.9.0       & Inclusion of a procedure suggested by Kaparis and Letchford (2010) to separate extended cover inequalities\\\bottomrule
 *     \end{tabularx}
 * \end{center}
 * \endlatexonly
 *
 * \section{Example of usage}
 * The following is an example of usage for the capacitated facility location problem
 * \code{.c}
 * VICcut* All_Generated_Cuts;    // pointer to first cut in a linked list containing all cuts generated so far //
 * VICcut* Cuts_of_Current_Round; // pointer to first cut in a linked list containing cuts generated in current round of generating cutting planes //
 * VICcut* Curr_Cut;              // pointer to first cut in a linked list of cuts generated by one of the routines described above //
 * int m, n;                      // number of customers and pot. depot sites //
 * double **cost;                 // cost[i][j] = cost of supplying all of customer i's demand from facility j //
 * double *fixcost;               // fixed depot costs //
 * int *demand, *capaci;         // customer demands and depot capacities //
 * double *x, *y;                // pointer to fractional solution (x,y) //
 * All_Generated_Cuts = NULL;
 * Cuts_of_Current_Round = NULL;
 * Read_CFLP_Data_and_Allocate_Mem(&m, &n, demand, capaci, cost );
 *
 * x = (double*) calloc( m*n, sizeof(double) );
 * y = (double*) calloc( n, sizeof(double) );
 *
 * // (x,y) denotes LP-solution of a CFLP
 * // x[i*n+j] is percentage of customer i's demand met from facility j in this solution, where  m = #customers, n =#depots //
 * Solve_First_LP( Data, x, y );
 *
 * if ( xy_Is_Integer ) return(); // LP solution is integer //
 *
 * do {
 *   Cuts_of_CurrentRound = NULL;
 *
 *   VICcflfc( m, n, demand, capaci, x, y, &Curr_Cut );
 *   if ( Curr_Cut ) VICaddtolst ( &Cuts_of_Current_Round, Curr_Cut );
 *
 *   VICcflsmi( m, n, demand, capaci, double* x, double* y, &Curr_Cut );
 *   if ( Curr_Cut ) VICaddtolst ( &Cuts_of_Current_Round, Curr_Cut );
 *
 *   VICuflohi( m, n, x, y, &Curr_Cut );
 *   if ( Curr_Cut ) VICaddtolst ( &Cuts_of_Current_Round, Curr_Cut );
 *
 *   Select_Cuts_You_Think_Are_Helpful ( &Cuts_of_Current_Round );
 *   if ( Cuts_of_Current_Round ) {
 *     Add_Cuts_to_Current_LP_Relaxation( Cuts_of_Current_Round, ... );
 *     VICaddtolst( &All_Generated_Cuts, Cuts_of_Current_Round );
 *   }
 *
 * } while ( Cuts_of_Current_Round != NULL );
 *
 * // Everythink done, free memory //
 * VICfreelst( &All_Generated_Cuts );
 * \endcode
 */

#ifndef __VICO_H
#define __VICO_H

#ifndef __CPXDEFS_H
#include <ilcplex/cplex.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct VICcut {
  char    sense;           //!< sense of inequality: 'L' means <=, 'G' is >=
  double  rhs;             //!< right-hand side of inequality
  int     nzcnt;           //!< number of non-zeros in the inequality
  int*    nzind;           //!< column/variable indices of the non-zeros
  double* nzval;           //!< values of the non-zeros in the inequality
  double  UsrRVal;         //!< may be used to store e.g. a dual variable
  int     UsrIVal;         //!< may be used to store e.g. an integer flag
  void*   UsrDatPtr;       //!< pointer to any additional user data
  struct  VICcut* nextcut; //!< pointer to next cut
} VICcut;

/*-------------------------------------------------------------------------------
The following is used for defining parameters for Fenchel cut generation
-------------------------------------------------------------------------------*/
#define SG_DEFAULT 0
#define SG_DEFLECT 1
#define SG_SMOOTH  2
#define ALG_MT1    0
#define ALG_DP     1
#define ALG_CPLEX  2

typedef struct TFENCHELopt {
  int    maxit;    /* maximum number of subgradient iterations. Default: 100 */
  int    sg_strat; /* determines the subgradient strateqy (Default: 2)
                      SG_DEFAULT: standard subgradient
		      SG_DEFLECT: subgradient deflection
		      SG_SMOOTH : exponential smoothing */
  double alpha;    /* Intial value of step size parameter (Default: 2) */
  int    H;        /* step size parameter is halfed if no improve occured in
                      lower bound after H iterations. If H=0 then the parameter
		      is continuously declined by a factor of 1.05. */
  int    CHK;      /* if > 0 then subgradient iterations are written to a file
                      named "fenchel.log". Default=0 */
  int    Algo;     /* Algorithm used for solving the single-node flow problems:
                      ALG_MT1:    ssfctp_mt1 (implicit enumeration)
		      ALG_DP:    ssfctp_dp  (dynamic programming)
		      ALG_CPLEX: cplex MIP solver
		      (Default = ALG_MT1) */
  int    Freq;     /* Frequency with which these cuts may be tried, e.g.,
                      only if no other cut found (Freq=0) or in every iteration
		      (Freq=1). Default Freg=1 */
  int    Reduce;   /* If Reduce=1, the Fenchel cut is generated for the reduced polytope
                      with all variables (x_j,y_j) with y_j=0 removed. This gives also
                      a valid inequality for the full polytope*/
} TFENCHELopt; //!< The following is used for defining parameters for Fenchel cut generation

/*!
 * Checks if a vector is almost integer
 * \param n integer. Size of the array
 * \param pi pointer to double array. Contains the values of the array
 */
int VICisintvec( int n, double* pi );


/*!
 * Sets the above mentioned parameters to their above mentioned default values
 */
void VICsetdefaults(  );


/*!
 * PURPOSE: Adds a linked list of cuts to an existing linked list of cuts. The new cuts
 * are inserted at the top of the existing list of cuts. Let
 *
 * First -> Second -> ... -> Last -> NULL
 *
 * be the list of already existing cuts and let
 *
 * Firstnew -> Secondnew -> ... -> Lastnew -> NULL
 *
 * denote the linked list of new cuts. After calling the procedure, the old list looks like
 *
 * Firstnew -> Secondnew -> ... > Lastnew -> First -> Second -> Last -> NULL
 *
 * However, the memory allocated for the cuts in the new list is not copied!
 * Therefore, do not delete it after addition to the old list.
 *
 * \param first    : pointer to the pointer to first cut in the existing list of cuts (if the list is empty *first must be NULL)
 * \param firstnew : pointer to first cut in the new list of cuts (which may consists of just one cut)
 *
 * \code{.c}
 * VICcut *MyCutsSoFar = NULL;
 * VICcut *MyNewCuts = NULL;
 *
 * ProcedureForGeneratingNewCuts( &MyNewCuts );
 * if ( MyNewCuts != NULL ) VICaddtolst( &MyCutsSoFar, MyNewCuts );
 * \endcode
 */
void VICaddtolst( VICcut** first, VICcut* firstnew );


/*!
 * Frees the memory allocated for a cut to which the pointer *cut points
 * \param cut  pointer to the pointer to the
 */
void VICfreecut ( VICcut** cut );


/*!
 * Frees the memory allocated for a linked list of cuts and empties the list
 * \param first : pointer to the pointer to the first cut in the list
 *
 * \code{.c}
 * VICcut* MyCutList;
 * many very strong and helpful cuts generated and optimum proven
 * VICfreelst( &MyCutList );
 * \endcode
 */

void VICfreelst ( VICcut** first );


/*!
 * Allocates memory required to store data of a cut
 *
 * \param numnz  number of nonzero coefficients in the cut
 * \param cut    on completion *cut is the pointer to the cut
 * \return 1 on success and 0 otherwise
 * \code{.c}
 * VICcut* MyCutPtr;
 * int     numnz=1000;
 * VICallocut( numnz, &MyCutPtr );
 * \endcode
 */

char VICallocut( int numnz, VICcut** cut );

/*!
 * Searches for the cut to which the pointer "cut" is pointing to in a given list of cuts
 * \param cutlst  pointer to the first cut in the linked list of cuts
 * \param cut pointer to the cut which is search in the list
 * \return NULL if the cut is not contained in the list, and the pointer to the cut in the list of cuts otherwise
 */
VICcut* VICsearchcut ( VICcut* cutlst, VICcut* cut );

/*!
 * Eliminates duplicate cuts from a linked list of cuts such that every cut
 * is only contained once in that list
 * \param first pointer to the pointer to the first cut in the list to be cleared
 */
void VICclearlst( VICcut** first );

/*!
 * Sorting of arrays of integers/doubles
 * \param n  dimension of the array that has to be sorted
 * \param ascending sort in ascending (descending) order if ascending = 1 (0)
 * \param doinit if doinit=1 it is assumed that the array "order" is not initialised
 * \param size if size=sizeof(int) it is assumed that numbers points to integer array otherwise numbers must point to array of doubles
 * \param numbers pointer to array of integers/doubles of dimension of at least n
 * \param order on output the sorted array is given by numbers[ order[0],...,order[n-1] ]
 */
void VICsort( int n, int ascending, int doinit, int size, void* numbers, int* order );

/*!
 * Tries to generate a (single) extended flow cover inequality for the
 * CFLP, which is violated by the current solution (x,y). The procedure
 * may also be used to generate extended flow cover inequalities
 * for the single-node flow problem given by
 * \f{align}{
 *	    &\sum_j z_j = d\\
 *	    &z_j \leq capaci_j*y_j\\
 *	    &0 \leq z_j \leq \min\{d,capaci_j\}\\
 *	    &y_j \in \{0,1\},\quad \forall j\\
 * \f}
 * In this case, call the procedure with m=1, demand = d, x=z/d
 * \param n number of potential depot sites
 * \param m number of customers
 * \param demand pointer to an array of integers of size of at least m containing the customers' demands
 * \param capaci pointer to an array of integers of size of at least n containing the depot capacities
 * \param x pointer to an array of doubles of size of at leat m*n containing the allocation part of the fractional solution, which should be
 *          separated by a flow cover inequality. Let i=0,...,m-1 and j=0,...n-1 be the indices of customers and depot sites, resp. Then x[i*n + j] is the solution value of
 *          the allocation variable x(i,j), where 0 <= x(i,j) <= 1. The variable x(i,j) denotes the fraction of customer i's demand met from facility j.
 * \param y pointer to an array of doubles of size of at least n containing the location part of the fractional solution, which should be separated by a flow cover inequality.
 *          0 <= y[j] <= 1 and x(i,j) <= y[j]
 * \param fc_cut pointer to a pointer to a cut. If no violated flow cover is found, the null pointer is returned; otherwise **fc_cut contains the cut.fc_cut->nzval contains
 *          the nonzeros of the cut, and fc_cut->nzind contains the column indices of the nonzeros, where the index of value i*n+j corresponds to the allocation variable
 *  	    x(i,j) and the index m*n+j corresponds to the location variable y(j).
 *
 * \code{.c}
 * VICcut* MyCutList = NULL;
 * VICcut* MyNewCut = NULL;
 * int m, n, *demand=NULL, *capaci=NULL;
 *
 * Read_Problem_Data_and_Allocate_Space( m, n, demand, capaci, ... );
 * Solve_Something_like_the_LP_Relaxation_and_obtain_x_y;
 *
 * VICcflfc( m, n, demand, capaci, x, y, &MyNewCut );
 * if ( MyNewCut != NULL ) VICaddtolst( &MyCutList, MyNewCut );
 * \endcode
 */
void VICcflfc( int m, int n, int* demand, int* capaci, double* x, double* y, VICcut** fc_cut );


/*!
 * Generates special types of "submodular inequalities" for the CFLP using a separation heuristic of Aardal. Several such inequalities, which cut off
 * the solution (x,y), may be returned.
 * \param n number of potential depot sites
 * \param m number of customers
 * \param demand pointer to an array of integers of size of at least m containing the customers' demands
 * \param capaci pointer to an array of integers of size of at least n containing the depot capacities
 * \param x pointer to an array of doubles of size of at leat m*n containing the allocation part of the fractional solution, which should be
 *          separated by a submodular inequality. Let i=0,...,m-1 and j=0,...n-1 be the indices of customers and
 *	        depot sites, resp. Then x[i*n + j] is the solution value of the allocation variable x(i,j), where 0 <= x(i,j) <= 1. The
 *	        variable x(i,j) denotes the fraction of customer i's demand met from facility j.
 *
 *
 * \code{.c}
 * VICcut* MyCutList = NULL;
 * VICcut* MyNewCut = NULL;
 * int m, n, *demand=NULL, *capaci=NULL;
 *
 * Read_Problem_Data_and_Allocate_Space( m, n, demand, capaci, ... );
 * Solve_Something_like_the_LP_Relaxation_and_obtain_x_y;
 * VICcflsmi( m, n, demand, capaci, x, y, &MyNewCut );
 * if ( MyNewCut != NULL ) VICaddtolst( &MyCutList, MyNewCut );
 * \endcode
 */
void VICcflsmi( int m, int n, int* demand, int* capaci, double* x, double* y, VICcut** first_smi );


/*!
 * Generates odd-hole inequalities violated by the solution (x,y) for the UFLP. Several such cuts may be returned.
 * \param n number of potential depot sites
 * \param m number of customers
 * \param x pointer to an array of doubles of size of at leat m*n containing the allocation part of the fractional solution.
 *	      Let i=0,...,m-1 and j=0,...n-1 be the indices of customers and depot sites, resp. Then x[i*n + j] is the solution value of
 *        the allocation variable x(i,j), where 0 <= x(i,j) <= 1. The variable x(i,j) denotes the fraction of customer i's demand met from facility j.
 * \param y pointer to an array of doubles of size of at least n containing the location part of the fractional solution, which should be
 *	      separated by a submodular inequality. 0 <= y[j] <= 1, y[j]>=x[i,j]
 * \param first_ohi : pointer to the first cut of a linked list of violated odd-hole inequalities found by this procedure. The NULL pointer is returned
 *	      if no violated inequality was found. See routine VICcflfc regarding definition of column indices.
 *
 *
 * \code{.c}
 * VICcut* MyCutList = NULL;
 * VICcut* MyNewCut = NULL;
 *
 * Solve_Something_like_the_LP_Relaxation_and_obtain_x_y;
 *
 * VICuflohi( m, n, x, y, &MyNewCut );
 * if ( MyNewCut != NULL ) VICaddtolst( &MyCutList, MyNewCut );
 *\endcode
 */
void VICuflohi( int m, int n, double* x, double* y, VICcut** first_ohi );

/*!
 * Tries to find violated combinatorial inequalities for the uncapacitated facility location problem.
 * Let K be a subset of the set of all customers, and let J denote a subset of the set of potential depot sites. Define a binary matrix \f$a(i,j)\f$ for
 * each \f$i \in K\f$ and \f$j \in J\f$. Let b denote (a lower bound on) the minimum number of depots \f$j \in J\f$ required to cover each customer \f$i \in K\f$, that is

 * \f{align}{
 *  b =   \min & \sum_{j\in J} y_j\\
 *      s.t.: & \sum_{j\in J} a_{ij} y_j \ge 1,\quad  \forall i \in K\\
 *            & y_j = 0,1 \forall j \in J
 * \f}
 * Then the inequality  \f$ \sum_{i \in K} \sum_{j\in J} a_{ij} x_{ij} - \sum_{j\in J} y_j \leq |K| - b\f$
 * is valid for the UFLP. These inequalites generalize odd holes for the UFLP and have been proposed by D.C. Cho et al. (1983): "On the uncapacitated
 * facility location problem I: Valid inequalities and facets", Mathematics of Operations Research 8, 579-589. See also G. Cornuejols, J.-M. Thizy
 * (1982): "Some facets of the simple plant location polytope", Mathematical Programming 23, 50-74.
 * Let \f$(x*, y*)\f$ denote a fractional solution. In order to find a violated inequality of the above type, the following simple heuristic is tried:
 * \f{enumerate}{
 *   \item Set $ J = \{ j : 0 < y^*_j < 1 \} $
 *   \item Set $ K = \{ i : x^*_{ij} > 0 for more than one j \in J \} $
 *   \item Solve the covering problem in order to find b (alternatively, find a lower bound on b by solving the LP relaxation of the covering
 *	        problem and rounding up the objective function value)
 *   \item Check if point $ (x*,y*)$ violates the resulting inequality
 * \f}
 *
 * \param n number of potential depot sites
 * \param m number of customers
 * \param x pointer to an array of doubles of size of at leat m*n containing the allocation part of the fractional solution.
 *	      Let \f$i=0,...,m-1\f$ and \f$j=0,...n-1\f$ be the indices of customers and depot sites, resp. Then x[i*n + j] is the solution value of
 *   	  the allocation variable x(i,j), where 0 <= x(i,j) <= 1. The variable x(i,j) denotes the fraction of customer i's demand met
 *	      from facility j.
 * \param y pointer to an array of doubles of size of at least n containing the location part of the fractional solution, which should be
 *	      separated by a submodular inequality. 0 <= y[j] <= 1, y[j]>=x[i,j]
 * \param first_cut pointer to cut found by this procedure. The NULL pointer is returned if no violated inequality was found. See routine VICcflfc regarding
 *	      definition of column indices.
 */
void VICuflcov( CPXENVptr Env, int m, int n, char SolveCov, double* x, double* y,
                VICcut** first_cut );

/*!
 * Clique generation and coefficient improvement (reduction) for a single knapsack inequality of the form
 * \f$ \sum_j a_j x_j <= b  or  \sum_j a_j x_j >= b\f$ where \f$x_j=0,1\f$ for all \f$j\f$
 *	The procedure tries to generate a single clique inequality from the knapsack inequality and to reduce the right-hand side together with
 *	the coefficients of variables in the clique.
 *
 * \param WHATRED Determines which coeffiecient improvement scheme is applied. WHATRED = 0 -> no coefficient improvement at all. WHATRED = 1 -> just simple coefficient improvement
 *        (increase w_j to c if x_j=1 implies all other x_k = 0). WHATRED = 2 -> only apply reduction of coefficients of variables in the clique.
 *         WHAT_RED = 3 -> do both types of improvement if possible.
 * \param TRYCLI If TRYCLI=0 no cliques are generated, otherwise cliques are derived if possible. ( In order to use the coefficient reduction, that is WHATRED >= 2;
 *        a clique is required and TRYCLI must equal 1).
 * \param n number of (free) variables in the knapsack.
 * \param cap pointer to an integer containing the capacity of the knapsack (right-hand side of inequality). *cap is possible changed (reduced).
 * \param sense  sense of inequality ('L' for <= and 'G' for >=).
 * \param weight pointer to an array of integers of length of at least n containig the coefficients ("weights") of (free) variables in the inequality. Some of the weights[j]
 *        may be changed (reduced).
 * \param order NULL or a pointer to an array of integers of length of at least n containing an ordering of the items 0,...,n-1 according to increasing weights.
 * \param indx NULL or a pointer to an array of integers of length of at least n containing indices of the variables in the inequality. If indx=NULL, it is assumed that the variables are numbered from 0, ..., n-1
 * \param clique pointer to the the generated clique cut
 */
void VICkpreduce( int WHATRED, int TRYCLI, int n, int* cap, char sense, int* weight, int* order, int* indx, VICcut** clique );


/*!
 * Heuristic procecedure of Gu, Nemhauser and Savelsbergh for generating lifted cover inequalties from a knapsack structure like
 * \f{align}{
 *   & \sum_{j \in N} a_j x_j \leq b\\
 *	 & 0 \leq x_j \leq 1   \forall j \in N\\
 *   &   x_j \in \{0,1\} \forall j \in N
 * \f}
 * or
 * \f{align}
 *   & \sum_{j \in N} a_j x_j \geq b\\
 *	 & 0 \leq x_j \leq 1   \forall j \in N\\
 *	 & x_j \in \{0,1\} \forall j \in N
 * \f}
 * See: Gu Z, Nemhauser GL, Savelsbergh MWP (1998). ``Lifted cover inequalities for 0-1 linear programs: Computation''. INFORMS J. on Computing 10:427-437
 * EXAMPLE:
 * Consider the following MIP
 * \f{align}{
 *  \max   & 7x_0 + 3x_1 + 6x_2 + 9x_3 + 10x_4 + 6x_5 + 8x_6 + 9x_7\\
 *  s.t.:  &   x_0 + 2 x_1 +   x_2                   \leq 2\\
 *         & 3 x_3 + 5 x_4 + 2 x_5 + 6 x_6 + 7 x_7 \leq 11\\
 *	 x_j\in\{0,1\}, \forall j
 * \f}
 * Furthermore, assume that variables x(3) is fixed to one and variable \f$x_5\f$ is fixed to zero. In order to derive a lifted cover inequality from
 * the second inequality, the procedure above has to be called in the following way:
 *
 * n = 3, cap = 11-3 = 8, sense = 'L', weight = (5, 6, 7 ), indx = (4, 6, 7), xlp = \f$( x_4, x_6, x_7 )\f$, rco = ( |redcost(4)|, |redcost(6)|, |redcost(7)| )
 *
 * \param n number of (free) variables appearing in the knapsack inequality.
 * \param cap right-hand side of the knapsack inequality.
 * \param sense sense (that is 'L' for <= or 'G' for >=) of the knapsack inequality.
 * \param weight coefficient a_j of (free) variables in the knapsack inequality (coefficients a_j are not restricted to be nonnegative!).
 * \param idx indices of the (free) variables in the knapsack inequality. If null it is assumed that variables are numbered from 0 to n-1.
 * \param xlp LP solution of (free) variables appearing in the knapsack inequality.
 * \param rco bsolute values of reduced costs of variables appearing in the knapsack inequality
 * \param lci pointer to the generated cut pointer.
 *
 */
void VIClci( int n, int cap, char sense, int* weight, int* indx, double* xlp, double* rco, VICcut** lci );


/*!
 * Tries to get a (1,k)-configuration inequality using the separation heuristic of Crowder, Johnson, Padberg in Oper. Res. 31 (1983).
 * For an alternative separation heuristic for (1,k)-configurations see Carlos E. Ferreira (1997). On Combinatorial Optimization Problems arising in
 * Computer Systems Design. Phd Thesis, Technische Universit\E4t Berlin.
 * Given the Knapsack-Polytop \f$ X = { x : \sum_{j\in N} w[j]*x[j] <= c, x_j = 0,1 }\f$ a (1,k)-configuration is a set \f$ NP \cup {z}\f$, where \f$ NP \subset N\f$, such that
 * \f{enumerate}{
 * \item  $\sum_{j \in NP} w[j] <= c$
 * \item The set $K \cup {z}$ is a cover with respect to N for all subsets K of NP with cardinality k
 * \f}
 * The corresponding (1,k)-configuration inequality is given by
 * \f{equation}{
 *       (r - k + 1)x[z] + \sum_{j\in NP} x[j] <= r, where r=|NP|
 * \f}
 * Crowder, Johnson, Padberg propose the following separation heuristic:
 * \f{enumerate}{
 * \item 1. Let $S \subset N$ be the cover, and $z\in S$ the item with maximum weight. Set $NP = S-{z}$ and $k=\vert NP\vert$.
 * \item For all $j \in N-S $ with $\sum_{l \in N-S} w[l] <= c$ do:
 *    \begin{enumerate}
 *       \item Check if $K \cup {z}$ is a cover for any $K \subseteq NP$ with $\vert K\vert=k$
 *       \item If this is the case build the corresponding (1,k)-configuration inequality and lift it.
 *       \item If the lifted inequality is violated, add the inequality
 *     \end{enumerate}
 *   \f}
 * \param n number of (free) variables appearing in the knapsack inequality.
 * \param cap right-hand side of the knapsack inequality.
 * \param sense sense (that is 'L' for <= or 'G' for >=) of the knapsack inequality
 * \param weight coefficient a_j of (free) variables in the knapsack inequality (coefficients a_j are not restricted to be nonnegative!).
 * \param indx indices of the (free) variables in the knapsack inequality. If null it is assumed that variables are numbered from 0 to n-1.
 * \param xlp LP solution of (free) variables appearing in the knapsack inequality
 * \param rco absolute values of reduced costs of variables appearing in the knapsack inequality
 * \param kconf	pointer to the first cut in a linked list of generated k-configuration cuts
 */
void VICkconf( int n, int cap, char sense, int* weight, int* indx, double* xlp, double* rco, VICcut** kconf );

/*!
 * Separation procedure of Gabrel & Minoux for finding most violated extended cover inequality. See: V. Gabrel, M. Minoux (2002).
 * A scheme for exact separation of extended cover inequalities and application to multidimensional knapsack problems. Oper. Res. Lett. 30:252-264.
 * For an example see VIClci
 *
 * \param n number of (free) variables appearing in the knapsack inequality.
 * \param cap right-hand side of the knapsack inequality.
 * \param sense sense (that is 'L' for <= or 'G' for >=) of the knapsack inequality.
 * \param weight coefficient a_j of (free) variables in the knapsack inequality (coefficients a_j are not restricted to be nonnegative!).
 * \param order NULL or an ordering of the items in the knapsack according to increasing weights.
 * \param indx indices of the (free) variables in the knapsack inequality. If null it is assumed that variables are numbered from 0 to n-1.
 * \param xlp LP solution of (free) variables appearing in the knapsack inequality.
 * \param eci pointer to the generated cut pointer.
 */

void VICeci( int n, int cap, char sense, int* weight, int* order, int* indx, double* xlp, VICcut** eci );

/*!
 * Procedure for generating extended cover inequalities as suggested in K. Kaparis, A.N. Letchford (2010). Separation algorithms for 0-1 knapsack
 * polytopes, Math. Prog. 124:69-91.
 *
 * \param n number of variables appearing in the knapsack inequality.
 * \param cap right-hand side of the knapsack inequality.
 * \param do_exact if equal to 1, exact separation is tried. This requires to repeatedly solve a 0-1 knapsack problem. If equal to 0 these knapsack problems
 *        are solved heuristically using a greedy method.
 * \param sense sense (that is 'L' for <= or 'G' for >=) of the knapsack inequality
 * \param weight coefficient a_j of variables in the knapsack inequality
 * \param indx indices of the variables in the knapsack inequality. If NULL, it is assumed that variables are numbered from 0 to n-1.
 * \param xlp LP solution of (free) variables appearing in the knapsack inequality
 * \param eci pointer to the generated cut pointer (NULL if none found)
 */
void VICecikl( int n, int cap, int do_exact, char sense, int* weight, int* indx,  double* xlp, VICcut** eci );

/*!
 * Given the GAP ( or alternatively LEGAP )
 * \f{align}{
 *       \max   & \sum_{i\in I} \sum_{j\in J} c_{ij} x_{ij}\\
 *        s.t.: & \sum_{i\in I} x_{ij} = 1,\quad \forall j \in J\\
 *              & \sum_{j\in J} w_{ij} x_{ij} \leq s_i,\quad \forall i \in I\\
 *	            & x_{ij} \in \{ 0,1\} \forall i\in I,j\in J
 * \f}
 * the procedure uses a heuristic described in Farias IR, Nemhauser GL (2001), a family of inequalities for the generalized assignment polytope, Oper. Res.
 * Letters 29, for finding a violated inequality of the form
 * \f{equation}{ \sum_{j\in J_i} w_{ij} x_{ij} + \sum_{j\in J_i} a_j \sum_{k\ne i} x_{kj} \leq  s_i\f}
 * where \f$ J_i \subseteq J, a_j = s_i - ( W_i - w_{ij} )\f$ and \f$W_i = \sum_{j\in J_i} w_{ij}\f$
 * \param m number of agents.
 * \param n number of jobs.
 * \param IsGap set equal to 1 if the problem is a GAP, that is if every job has to be assigned to exactly one agent.
 *        Otherwise a LEGAP is assumed, that is some jobs may be not assigned to any agent.
 * \param capaci pointer to an array of integers of length of at least m containing the agents' capacities.
 * \param weight: pointer to an array of pointers of length of at least m such that weight[i][j] gives the amount of resources required by agent i to perform job j.
 * \param X pointer to an array of doubles of length of at least m*n containing the fractional LP solution, X[i*n+j] must contain the LP value of assignment variable x_{ij}
 * \param indx if not NULL this array gives the indices of the assignment variables in the user's problem formulation, that is indx[i*n+j] is the index of variable x_{ij}.
 *        If indx==NULL it is assumed that variables are numbered from 0 to m*n-1, where i*n+j is the index of variable x_{ij}.
 * \param fn_cut pointer to the first cut in a linked list of generated cuts of this type
 */
void VICgapfn( int m, int n, int IsGap, int* capaci, int** weight, double* X, int* indx, VICcut** fn_cut );

/*!
 * Subroutine for exact knapsack separation. Given a fractional solution x*, the most violated valid inequality
 * \f$ \pi x \le 1\f$ is returned by optimizing over the 1-polar of the knapsack polytope, that is by solving the separation problem
 * \f{equation}{
 *    \max\{ \pi x^* : \pi x^t \le 1 for each feasible solution x^t \}.
 * \f}
 * The dual of the this separation problem is solved by means of column generation. In order to reduce the effort, the separation is performed for a reduced knapsack
 * polytope obtained by resp. fixing variables \f$x_j\f$ to 0 and 1 if \f$x^*_j=0\f$ or \f$x^*_j=1\f$. Afterwards, a sequential lifting of fixed variables is applied in order to get a valid cut.
 *
 * \param Env pointer to the CPLEX environment as returned by CPXopenCPLEX. Env must be a valid pointer to the open CPLEX environment
 * \param justUp if equal to one, only variables showing a value of zero in the LP solution are first fixed to zero. The resulting inequality is then just uplifted to obtain
 *        a (strengthened) inequality for the full polytope. Otherwise (justUp=0) both variables showing LP-value of zero and one are fixed. The resulting inequality is then
 *        first to be downlifted to obtain a valid inequality, and thereafter uplifting is applied.
 * \param n number of variables appearing in the knapsack inequality.
 * \param cap right-hand side of the knapsack inequality.
 * \param sense sense (that is 'L' for <= or 'G' for >=) of the knapsack inequality.
 * \param weight coefficient a_j of variables in the knapsack inequality (coefficients a_j are not restricted to be nonnegative!).
 * \param indx indices of the variables in the knapsack inequality. If NULL it is assumed that variables are numbered from 0 to n-1
 * \param xlp LP solution of variables appearing in the knapsack inequality.
 * \param rco absolute values of reduced costs of variables appearing in the knapsack inequality.
 * \param cut pointer to the generated cut pointer.
 */
void VICkpsep( CPXENVptr Env, int justUp, int n, int cap, char sense, int* weight, int* indx, double* xlp, double* rco, VICcut** cut );

#ifdef __cplusplus
}
#endif

#endif // ifdef VICO_H
