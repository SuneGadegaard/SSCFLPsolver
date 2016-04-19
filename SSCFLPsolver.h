#ifndef SSCFLPSOLVER_H_INCLUDED
#define SSCFLPSOLVER_H_INCLUDED

/** \mainpage An introduction to the SSCFLPsolver class
* \author Sune Lauth Gadegaard
* \version 1.2.0
*
* \section License
*
* Copyright 2015, Sune Lauth Gadegaard.
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
* If you use the software in any academic work, please make a reference to
*
* S.L. Gadegaard and A. Klose and L.R. Nielsen, (2015), ``An exact algorithm based on cutting planes
* and local branching for the single source capacitated facility location problem'', Technical report, CORAL, Aarhus University, sgadegaard\@econ.au.dk.
*
* \section Description
*
* Class implementing a solver for the single sourced capacitated facility location problem (SSCFLP)
* The algorithm implememnted here is presented in Gadegaard et al. (2015) ``An exact algorithm based on cutting planes
* and local branching for the single source capacitated facility location problem''. The idea of the algorithm
* is based on the cut--and--solve framweork presented in S. Climer and W. Zhang, (2016), ``Cut-and-solve: An iterative search strategy for combinatorial optimization problems'',
* Artificial Intelligence, 170:714-738. We iterate between solving a constrianed relaxation (the Dense problem) of SSCFLP and a restricted version of SSCFLP (the Sparse problem).
* This leads to a framework where a constrained CFLP is used as the Dense problem and a smaller SSCFLP is used for the Sparse problem.
*
* \section formulation ILP formulation
*
* \latexonly
    Given a set of potentil facility sites, $I$, and a set of customers, $J$, the SSCFLP can be formulated as follows: Let $f_i>0$ be the cost of opening facility $i$, and $c_{ij}$ the cost of
    assigning customer $j$ to facility $i$. Furthermore, let $d_j>0$ the demand for a certain good at customer $j$ and $s_i>0$ the capacity for the good at facility $i$. Introducing binary
    variable $y_i$ which is equal to one if and only if a facility is open at facility $i$ and binary variable $x_{ij}$ indicating if customer $j$ is assigned to facility $i$ ($x_{ij}=1$) or
    not ($x_{ij}=0$), the SSCFLP can be formulated as the linear integer programming problem
    \begin{align}
        \min        \ & \sum_{i \in I} \sum_{j \in J} c_{ij}x_{ij} + \sum_{i \in I} f_iy_i\tag{O}\label{Obj}\\
        \text{s.t.:}\ & \sum_{i \in I} x_{ij} = 1, \quad \forall j \in J,\tag{A}\label{Ass}\\
                    \ & \sum_{j \in J} d_jx_{ij} \leq s_iy_i, \quad \forall i \in I,\tag{C}\label{Cap}\\
                    \ & x_{ij}-y_i \leq 0, \quad \forall i \in I,\ j \in J.\tag{GUB}\label{Gub}\\
                    \ & \sum_{i \in I} s_iy_i \geq D = \sum_{j \in J},\tag{TD}\label{Tot}\\
                    \ & x_{ij},\ y_i \in \{0,1\}, \quad \forall i \in I,\ j \in J.\tag{B}\label{Bin}
    \end{align}
    Here \eqref{Obj} minimizes the total cost, composed of assignment costs and fixed opening costs. Constraints \eqref{Ass} ensure that all customers are assigned to exactly one facility
    while the constraints \eqref{Cap} make sure that the no customer is assined to a closed facility and that each facility's capacity is respected. The generalized upper bounds \eqref{Gub}
    are infact redundant (implied by \eqref{Cap}), but the improve the LP relaxation quite a lot. The total demand constraint \eqref{Tot} is likewise redundant, but its structure is used
    in the cutting plane part of the algorithm implemented here. Lastly, \eqref{Bin} require all variables to be binary.
* \endlatexonly
*
* \section Compiling
* The codes were compiled using the GNU GCC compiler on a Linux Ubuntu 14.04 machine.
* The following flags were used: -Wall -O3 -std=c++11 -DIL_STD.
* The Code::blocks IDE was used as well. Information on how to configure
* Code::Blocks IDE with CPLEX on a Linux machin can be found here: <http://www-01.ibm.com/support/docview.wss?uid=swg21449771>
*
* \latexonly
*   \section{Example of usage}
*   This section contains two examples of how the SSCFLPsolver can be used. The first example loads data from a file provided as the first argument to the main function while the other
*   loads data into the SSCFLPsolver object using pointers to arrays.
*   \subsection{Loading data from a file}
*   This is a simple example using a datafile provided as an argument to the main function
*   \begin{lstlisting}
*       // main.cpp
*       #include"SSCFLPsolver.h"
*       int main(int argc, char** argv){
*           try{
*               SSCFLPsolver solver = SSCFLPsolver();
*               // the second argument specifies the format of the data file.
*               solver.Load(argv[1],1);
*               if(solver.Run()){
*                   std::cout << "Hurra!" << std::endl;
*               }else{
*                   std::cout << "Bummer!" << std::endl;
*               }
*               return 0;
*           }catch(std::exception &e){
*               std::cerr << "Exception: " << e.what() << std::endl;
*               exit(1);
*           }catch(...){
*          	    std::cerr << "An unexpected exception was thrown. Caught in main.\n";
*          	    exit(1);
*           }
*       }
*   \end{lstlisting}
*
*   \subsection{Loading data from pointers}
*   In this example we load the data using a user-provided function called \texttt{GetMyData}.
*   \begin{lstlisting}
*       //main.cpp
*       #include"SSCFLPsolver.h"
*       int main(){
*           try{
*               int n, m;
*               int* f, s, d;
*               int** c;
*               GetMyData((n, m, c, f, d, s);
*               SSCFLPsolver solver = SSCFLPsolver();
*               solver.Load(n, m, c, f, d, s);
*               if(solver.Run()){
*               std::cout << "Hurra!" << std::endl;
*               }else{
*                   std::cout << "Bummer!" << std::endl;
*               }
*           }catch(std::exception &e){
*               std::cerr << "Exception: " << e.what() << std::endl;
*               exit(1);
*           }catch(...){
*               std::cerr << "An unexpected exception was thrown. Caught in main.\n";
*               exit(1);
*           }
*       }
*   \end{lstlisting}
*
*   \subsection{Problematic behaviour}
*   One kind of not so brilliant thing about this implementation is that you cannot call the run function multiple times, as an IloConversion is used after the cutting plane phase.
*   The IloConversion can only be called once per model! A possible workaround would be to simply have a pure cutting plane IloModel internally in the class which could be dedicated
*   solely to the generation of cutting planes.
* \endlatexonly
*
* \latexonly
* \section{Change log for SSCFLPsolver.h and SSCFLPsolver.cpp}
* \begin{center}
*     \begin{tabularx}{\textwidth}{llr X}\toprule
*        FILE:          &   \multicolumn{3}{l}{SSCLPsolver.h and SSCFLPsolver.cpp}\\
*        Version:       &   \multicolumn{3}{l}{1.2.0}\\
*        \multicolumn{4}{l}{- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -}\\
*        CHANGE LOG:    &   DATE    &   VER.-NO.    &   CHANGES MADE\\ \midrule
*                       &   2015--03--01    &   1.0.0       & First implementation\\
*                       &   2015--04--09    &   1.1.0       & Cutting planes added\\
*                       &   2015--04--13    &   1.1.1       & Variable fixation in the sparse problem added.\\
*                       &   2015--04--14    &   1.1.1       & A number of get/set functions and a "setSolution" function which sets a solution internally in cplex was added.\\
*                       &   2015--04--01    &   1.1.2       & Added exact separation from the effective capacity polytope with generalized upper bounds. Does not seem to increase the bound that much\\
*                       &   2016--04--19    &   1.2.0       & Added functionality to solve the sparse problems using dual ascent. Works best when ratio between total capacity and total demand ins small (\leq 3)\\\bottomrule
*     \end{tabularx}
* \end{center}
* \endlatexonly
*/

#include<ilcplex/ilocplex.h>	// Header for IloCplex
#include<exception>		        // Header for exceptions
#include<stdexcept>	    	    // Header for exceptions (runtime_error, logical_error a.s.o)
#include<vector>                // Header for the vector class
#include"vico.h"	        	// Header for the Valid Inequalities for Combinatorial Optimization. Written by Andreas Klose, Department of Mathematics, Aarhus University.
#include"combo.h"
#include"solution.h"            // Header for the solution class.
#include<time.h>                // Header providing clock
#include<chrono>                // Header used for timing stuff
#include"rdDat.h"

typedef std::chrono::high_resolution_clock CPUclock;

typedef IloArray<IloNumVarArray> IloVarMatrix;  //!< An IloArray of IloNumVarArrays
typedef IloArray<IloNumArray> IloNumMatrix;     //!< An IloArray of IloNums

struct testStats{
    int n;                      //!< Number of facilities
    int m;                      //!< Number of customers
    long numberOfIterations;    //!< Number of iterations in the dual ascent algorithm
    double time;                //!< TotalTime
    long itWhereOptWasFound;    //!< The dual ascent iteration where the optimal solution was found. If equal to zero, optimal solution was found by initial heuristic
    double initialUpperBound;   //!< The initial upper bound before dual ascent starts
    double bestUpperBound;      //!< Best upper bound. Equal to optimal solution if no optimality gap is left
    double bestLowerBound;      //!< Best lower bound on the instance
    double percentageGap;       //!< Percentage gap between lower and upper bound. Calculated as ( UB - LB ) / LB * 100.
    double CuttingTime;         //!< Time used in the cutting plane algorithm
    double CutAndSolveTime;     //!< Time used in the cut and solve algorithm
    double avgNumOfDualItPerSparse; //!< The average number of dual ascent iterations per sparse problem
    double numOfSparseSolved;       //!< Total number of sparse problems solved
    double WeakLowerBound;      //!< Lower bound before adding cutting planes
    double LowerBoundAfterCusts;//!< Lower bound after applying the cutting plane algorithm
    double avgPercentLeft;      //!< Average number of assingment variables left after solving Sparse using dual ascent.
};

class SSCFLPsolver{
    private:
        /**
         * @name Cplex
         * This section contains all the cplex gear needed for the algorithm to run.
         */
        ///@{
        IloEnv env;             //!< The environment used throughout the lifetime of the object
        IloModel SparseModel;   //!< The IloModel used for the sparse problem
        IloModel DenseModel;    //!< The IloModel used for the dense problem
        IloCplex SparseCplex;   //!< The IloCplex environment used for the sparse problem
        IloCplex DenseCplex;    //!< The IloCplex environment used for the dense problem
        IloObjective DenseObj;  //!< Extractable holding the objective function of the dense problem
        IloObjective SparseObj; //!< Extractable holding the objective function of the sparse problem

        IloRangeArray AssCst;   //!< Holds the assignemnt constraints of the denseproblem only!
        IloNumArray duals;      //!< Holds the duals of the assingment constraints of the dense problem

        std::vector<std::pair<int,int>> fixAfterCuttingPhase;
        std::vector<std::pair<int,int>> fixBeforeCuttingPhase;
        ///@}

        /**
         * @name Data
         * This section contains all the data for describing the SSCFLP.
         */
        ///@{
        IloInt n;       //!< Number of facility sites.
        IloInt m;       //!< Number of customers.
        IloNumMatrix c; //!< Assignment costs. c[i][j] = cost of assigning customer j to facility i.
        IloNumArray f;  //!< Fixed opening costs. f[i] = cost of opening facility i.
        IloNumArray d;  //!< Demands. d[j] = demand at customer j.
        IloNumArray s;  //!< Capacities; s[i] = capacity at customer i.
        IloInt TD;      //!< Total demand. TD = sum_j d[j].
        ///@}

        /**
         * @name Parameters and flags
         * This section contains a list of parameters and flags used internally in the SSCFLPsolver class.
         */
         ///@{
        bool hasCleanedUp;  //!< Used to flag if the CleanUp function has been called.
        bool hasRun;        //!< Used to flag if the Run() function has been called.
        bool displayStats;  //!< Used to flag if statistics should be printed to screen.
        bool CplexOutOff;   //!< Used to flag that cplex' output should be redirected to the null stream.
        bool debugMe;       //!< Used to flag if debug print outs should be enabled.
        bool solveSparseUsingDualAscentAlg; //!< True if one should use dual ascent for solving the sparse problems. Defaults to false

        int ObjVal;         //!< After the Run() has been called, ObjVal contains the optimal solution value.
        int maxF;           //!< Contains the maximum value for the fixed opening costs after calling Load().
        int minF;           //!< Contains the minimum value for the fixed opening costs after calling Load().
        int maxC;           //!< Contains the maximum value for the assignment costs after calling Load().
        int minC;           //!< Contains the minimum value for the assignment costs after calling Load().
        int minD;           //!< Contains the minimum value for the customer demands after calling Load().
        int maxD;           //!< Contains the maximum value for the customer demands after calling Load().
        int minS;           //!< Contains the minimum value for the facility capacities after calling Load().
        int maxS;           //!< Contains the maximum value for the facility capacities after calling Load().

        double avgF;        //!< Contains the average value of the fixed costs after calling Load().
        double avgC;        //!< Contains the average value of the assignment costs after calling Load().
        double avgD;        //!< Contains the average value of the customer demands after calling Load().
        double avgS;        //!< Contains the average value of the facility capacities after calling Load().
        double StartTime;   //!< Clock time where the SSCFLPsolver object is initialized.
        double StartRunTime;//!< Clock time where the Run() function is called.
        double CutTime;     //!< CPU seconds used in the cutting phase
        double CnSTime;     //!< CPU seconds used in the cut and solve phase.
        double Runtime;      //!< CPU seconds used for the whole procedure.
        double CutLowerBound;//!< The lower bound produced by the cutting plane algorithm

        int* ySol;          //!< Pointer to an array of integers. If ySol[i]=1 facility i is open.
        int* xSol;          //!< Pointer to an array of integers. If xSol[j]=i customer j is assigned to facility i.
        int* NumCuts;      //!< Counts the number of cuts generated.
        ///@}

        /**
         * @name Tolerances
         * Tolerances used to control floating point arithmetic. All values are initialized in the constructor of the class.
         */
         ///@{
         double myZero;     //!< My interpretation of ``zero''. Constructor initialized
         double myOne;      //!< My interpretation of ``one''. Constructor initialized
         double myEpsilon;  //!< My interpretation of ``equal''. Constructor initialized
         ///@}

        testStats* stats; //!< Pointer to an instance of the testStats struct. Contains on termination the statistics gathered during optimization.

        // Private functions

        /*! \brief Builds the IloModels
         *
         * This function builds the model on the IloModels defined above and exports them to approproate cplex environments.
         */
        void BuildModels();

        /*! \brief Prints info about the program and the author
         *
         * Prints information about the program and the author and the program.
         * List of important parameter values is printed
         */
        void PrintProgramInfo();

        /*! \brief Separates an LP solution (x,y) from the capacity constraints
         *
         * Separates an LP solution (x,y) from the capacity constraints. Tries first with
         * a lifted cover inequality, then with an extended cover inequality and finally with a fenchel cut.
         */
        bool SepCuts();

        /*! \brief Runs the cutting plane phase.
         *
         * Runs the cutting plane phase corresponding to Phase 1, in S.L. Gadegaard and A. Klose and L.R. Nielsen, (2015), ``An exact
         * algorithm based on cutting planes and local branching for the single source capacitated facility location problem'',
         * Technical report, CORAL, Aarhus University, sgadegaard\@econ.au.dk.
         */
        void CuttingPhase();

        /*! \brief Checks if all y-variables can be fixed in the sparse problem
         *
         * Checks if all y-variables can be fixed in the sparse problem. Uses the combinatorial argument that if $sum_{i: 0<=y_i<=1}s_i -s_{\tilde{i}}<=D$ for some D,
         * then $y_\tilde{i}$ can be fixed to one.
         * \param ones std::vector of integers. Contains the free variables on input and the variables which can be fixed to one on output.
         */
        void CheckIfFix(std::vector<int> &ones);

        /*!
         * Change the dense problem to a semi lagrangean of the SSCFLP.
         * As Lagrangean dual multiplier, we use the largest dual cost of the assignemnt constraints of the cut--enhanced SSCFLP
         */
        void ChangeDenseToSemiLagrangean ( );

        /*! \brief Runs a simple local branching heuristic
         *
         * Runs a simple local branching heuristic that starts by refining the locational decision followed by a second phase where
         * the allocation of customers is gradually improved
         */
        void RunLBHeur ();

        /*! \brief Performs local branching on the locational variables
         *
         * Method intended to be used by the RunLBHeur() function. It performs local branching on the locational variables.
         * The first local branching constraints is based on the solution to the LP relaxation after adding cuts.
         * \param LBmod reference to an IloModel object. Used to store a copy of the SSCFLP we are solving
         * \param LBcpx reference to an IloCplex object. USed to solve the IloModel LBmod
         */
        void RefineLocation (  );

        /*! \brief Performs local branching on the allocation varibales
         *
         * Method used in the RunLBHeur() function after the RefineLocation() function has been called.
         * \param LBmod reference to an IloModel object. Used to store a copy of the SSCFLP we are solving
         * \param LBcpx reference to an IloCplex object. USed to solve the IloModel LBmod
         */
        void RefineAllocation (  );

        /*! \brief A greedy algorithm for SSCFLP
         *
         * Algorithm used if no feasible solution could be found. The method tries in a greedy way to assign the customer with largest demand
         * to facility with largest residual capacity. The method usually finds a (low quality) solution!
         * The method starts by solving the binary knapsack problem max{ sum(i) f(i)*z(i) : sum(i) s(i)*z(i)<= sum(i)s(i)-sum(j)d(j)} in order to
         * find an initial set of open facilities. An facility is initially open if z(i)=0 in an optimal solution to the knapsack problem.
         * Then the set of customers is sorted in non--increasing order of demand and the customers are added to the cheapest open facility
         * with enough capacity. If no open facility has enough capacity, the
         */
        void greedy ( );

        /*! \brief A dual ascent algorithm for the sparse problems
         *
         * Algortihm that implements a semi-Lagrangean based dual ascent algorithm for solving the sparse problems.
         * It start by writing the constraints $ \sum_{i\in I} x_{ij} =1$ as the two sets of constriants $ \sum_{i\in I} x_{ij} \leq 1$ and
         * $ \sum_{i\in I} x_{ij} \qeq 1$. The latter set is then first surrogate relaxed by multipliers $e=(1,\dots,1)$ and then the resulting constriaint
         * is relaxed in a Lagrangean manner by lagrangean multiplyer $\mu$. This results in a Lagrangean dual having only one variable and
         * a dual which closes the duality gap.
         * \param NewObjective reference to a double. Contains on output the new objective function value \emph{if} it improves the current best solution value.
         * \return true if the dual ascent algorithm solved the problem to optimality. Note, we return false if we could not find an improving solution!
         */
        bool solveSparseUsingDualAscent ( double & NewObjective , solution &incumbent );

        /*!
         * Two exchange local search. Swaps two customers if a gain is obtained
         */
        void twoSwap( );

        /*! \brief Cleans up the memory
         *
         * This function is used to clean up the class SSCFLPsolver in case of exceptions or other unpleasanties.
         * Also used as the base function in the destructor.
         */
        void CleanUp();



    public:
    // Public variables
         /**
         * @name Public variables
         * This section contains all publicly accesible variables
         */
        ///@{
            IloNumVarArray y;       //!< IloNumVarArray used for the location variables
            IloVarMatrix x;         //!< IloVarMatrix used for the assignment variables
        ///@}
	// Public functions
        /*! \brief The constructor of the class
         *
         * The constructor of the class SSCFLPsolver. Initializes the data containers, IloModels, and IloCplex environments.
         */
        SSCFLPsolver();

        /*! \brief The destructor of the class
         *
         * The destructor of the class SSCFLPsolver. Frees all memory allocated during the lifetime of the class object.
         */
        ~SSCFLPsolver();

        /*! \brief Public function for loading data
         *
         * Public function for loading data describing an instance of the SSCFLP. Data is given as a data file.
         * \param FileName pointer to char array. Contains the path (relative or absolute) to the data file.
         * \param format integer. Specifies the format of the data file. 0 = DiazFernandez. 1,2,3 = Holmberg, Yang and Gadegaard. Default is 1.
         */
        void Load(char* FileName, int format);

        /*! \brief Public function for loading data.
         *
         * Public function for loading data describing an instance of the SSCFLP. Data is given directly.
         * \param nn integer. The number of facility sites.
         * \param mm integer. The number of customers.
         * \param cc pointer to a pointer to array of integers. cc[i][j] = cost of assigning customer j to facility i. Must have dimmensions at least n times m
         * \param ff pointer to an array of integers. ff[i] = cost of opening facility at site i. Must have dimmensions at least n.
         * \param dd pointer to an array of integers. dd[j] = demand at customer j. Must have dimmensions at least m.
         * \param ss pointer to an array of integers. ss[i] = capacity at facility site i. Must have dimmensions at least n.
         */
        void Load(int nn, int mm, int **cc, int* ff, int* dd, int* ss);

        /*! \brief Public function for loading the data of a SSCFLP
         *
         * Public function for loading data describing an instance of the SSCFLP. Data is given as a const reference to a rdDat object.
         * \param data pointer to a rdDat object containing the data of a SSCFLP
         */
        void Load( rdDat* data );

        /*! \brief Hands a solution to the solver.
         *
         * This function can be used to set a (heuristic) solution internally in cplex.
         * \param nn integer. Length of the array yval. Must be less than or equal to the number of facility sites.
         * \param mm integer. Length of the array xass. Must be less than or equal to the number of customers.
         * \param yval pointer to an array of integers. yval[i]=1 iff facility i is open in the solution you are providing.
         * \param xass pointer to an array of integers. xass[j]=i iff customer j is assigned to facility i in the solution you provide.
         * \param UBval integer. Contains the objective function value of the solution you provide. Will be used as cutoff value by cplex.
         * \return true if the solution provided was infact feasible.
         * \note Should be called prior to calling run. Cplex is instructed to solve the problem as an LP with the solution given by (yval, xass) fixed.
         */
        bool setSolution(int nn, int mm, int* yval, int* xass, int UBval);

        /*! \brief Returns the optimal objective function value.
         *
         * Returns the optimal objective function value. Can only be called after the Run() function has been called.
         * An exception is thrown if getObjVal() is called prior to Run()
         */
        int getObjVal(){ if(!hasRun){throw std::runtime_error("No solution yet. Call Run() befor getObjVal()!\n");} return ObjVal; }

        /*! \brief Displays statistics for an instance of SSCFLP loaded into the SSCFLPsolver object
         *
         * Displays statistics for the instance of SSCFLP loaded into the SSCFLPsolver object. If the Load function has not been called,
         * the method will just print rubbish
         */
        void printStats();

        /*! \brief This is the main function in the SSCFLPsolver class.
         *
         * This is the main function in the SSCFLPsolver class. The function executes the cutting plane algorithm as well as
         * the cut and solve algorithm.
         * \return true if everything goes as planned. false otherwise.
         */
        bool Run();

        /*! \brief Runs a local branching heuristic
         *
         * Runs a local branching heuristic for the SSCFLP. First the a cutting plane algorithm is run. The y-variables which are positive in the cut-enhanced LP-relaxation
         * is used to indicate which facilities are likely to be in an optimal solution. A local branching constraint is added and cplex is run for a limited amount of time.
         * The solution found after this time-limited run is used as the first solution in a local branching constraint refining the locational deccision. That is, only local
         * branching is performed on the locational variables. When a satisfactory solution is found, a second phase is entered where the allocation-decions is refined. The final
         * solution can be accessed by getSolution( int*y, int* x).
         * \return True if a feasible solution is found, and false otherwise.
         */
        bool RunAsHeuristic ( );

        /*! \brief Method for retrieving an optimal solution to the instance of SSCFLP.
         *
         * Method for retrieving an optimal solution to the instance of SSCFLP. Note that one needs to call Run() before calling getSolution as there will be no solution otherwise.
         * \param getY Pointer to integer array of size at least n. Contains on output an optimal solution to the y-variables. getY[i]=1 if facility i is open in an optimal solution. If for some reason no solution is stored in the class, getY = 0 on output, that is the null-pointer is returned.
         * \param getX Pointer to integer array of size at least m. Contains on output an optimal solution to the x-variables. getX[j]=i if customer j is assigned to facility i in an optimal solution. If for some reason no solution is stored in the class, getX = 0 on output, that is the null-pointer is returned.
         */
        void getSolution(int* getY, int* getX);

        /*! \brief Returns the assignment cost c[i][j]
         *
         * Method returning the cost of assigning customer j to facility i
         * \param i integer. Index of the facility you want.
         * \param j integer. Index of the customer you want.
         * \note The index i has to be less than getNumFac() and j has to be less than getNumCust(). Both i and j need to be non-negative. No exceptions are thrown if you fuck it up!
         * \return The cost of assigning customer j to facility i
         */
        inline
        int getC(int i, int j){return c[i][j];}

        /*! \brief Returns the fixed opening cost f[i]
         *
         * \param i integer. The index of the facility who's fixed cost you want.
         * \note i>=0 and i<getNumFac()! No exceptions are thrown if you fuck it up!
         * \return The fixed incurred while opening facility i
         */
        inline
        int getF(int i){return f[i];}

        /*! \brief Sets the assignment cost c[i][j]
         *
         * Sets the cost of assigning customer j to facility i
         * \param i integer. Index of the facility.
         * \param j integer. Index of the customer.
         * \param cc integer. Value of the assignment cost c[i][j] after execution.
         * \note 0<= i,j,  i<getNumFac(), and j<getNumCust(). Note also that setC does not! change the objective of the IloCplex objects. It only changes the internal data.
         */
        void setC(int i, int j, int cc);

        /*! \brief Sets the fixed opening cost f[i]
         *
         * Sets the cost of opening facility i
         * \param i integer. Index of the facility.
         * \param ff integer. Value of the fixed cost f[i] after execution.
         * \note 0<= i < getNumFac(). Note also that setF does not! change the objective function coefficient of the IloCplex objects. It only changes the internal data.
             */
        inline
        void setF(int i, int ff){f[i] = ff;}

        /*! \brief Returns the number of facilities
         * \return The number of facility sites.
         */
        int getNumFac(){return n;}

        /*! \brief Returns the number of customers
         * \return The number of customers.
         */
        int getNumCust(){return m;}

        /*!
         * Returns the optimal value of the dual  variable of the specified assignment constraint.
         * The duals are of the cut-enhanced problem
         */
        inline
        IloNum getDual ( int j ) { return duals[j]; }

        /*!
         * Returns the demand of cutomer specified by index
         * \param indx integer. Index of the customer for which you want the demand
         */
        inline
        int getDemand ( int indx ){ return d[indx]; }

        /*!
         * Returns the capacity of the facility specified by the index
         * \param indx integer. Index of the facility for which you want the capacity
         */
        inline
        int getCapacity ( int indx ) { return s[indx]; }

        /*!
         * Returns the total demand of the instances. Thas is sum(j) d(j).
         */
        inline
        int getTotalDemand ( ) { return TD; }

        /*!
         * Returns the lower bound obtained by the cutting plane algorithm
         */
        inline
        double getCutLower ( ) { return CutLowerBound; };

        /*! \brief Fixes the x[i][j] variable to the zero before adding cuts
         *
         * Fixes the x[i][j] variable to the zero before adding cuts. This means that if a x[i][j] should be fixed to one, then all other assignments
         * must be set to zero for the corresponding customer.
         * \param i integer. Index of the facility.
         * \param j integer. Index of the customer.
         */
        void fixBeforeCuts ( int i, int j){ fixBeforeCuttingPhase.push_back ( std::pair<int,int>(i,j) ); }

        /*! \brief Fixes the x[i][j] variable to the zero after adding cuts
         *
         * Fixes the x[i][j] variable to the zero after adding cuts. This means that if a x[i][j] should be fixed to one, then all other assignments
         * must be set to zero for the corresponding customer.
         * \param i integer. Index of the facility.
         * \param j integer. Index of the customer.
         */
        void fixAfterCuts ( int i, int j){ fixAfterCuttingPhase.push_back ( std::pair<int,int>(i,j) ); }

        /*! \brief Returns the test statistics gathered during the optimization
         *
         * Returns the test statistics as a pointer to the struct-type testStats.
         */
        testStats* getTestStatistics ( ) { return stats; };

        /*! \brief Enables the dual ascent algorithm
         *
         * Enables the dual ascent algorithm for solving the sparse problems. It relaxes the assignment constraints in a semi-Lagrangean manner,
         * and increases the dual multiplier until an optimal solution has been found.
         */
        inline
        void setDualAcentOn ( ) { solveSparseUsingDualAscentAlg = true; }
};


#endif // SSCFLPSOLVER_H_INCLUDED
