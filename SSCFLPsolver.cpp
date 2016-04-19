#include"SSCFLPsolver.h"


bool initializeWithDual1 = false;
bool initializeWithPercent1 = true;
bool initializeWithSolution1 = false;


/*!
 * \brief Rounds double to two digits.
 * \param d Double. The double which should be rounded.
 */
double roundToTwo(double d){
    double rd = std::floor(d*100)/100;
    if(rd<=0.001) return 0.00;
    else return rd;
}

/*!
 * \brief Calculates the percentage gab between ub and lb.
 * \param ub Double. A number larger than lb, for example an upper bound on SSCFLP.
 * \param lb Double. A number smaller than ub, for example a lower bound on SSCFLP.
 */
double calcPercent(double ub, double lb){
    double percent = ((ub -lb)/lb)*100;
    return roundToTwo(percent);
}

using namespace std::chrono;
long cutsAddedInCutCallback = 0;
long doubleKPcutsAdded = 0;


ILOUSERCUTCALLBACK1( KnapsackSep , SSCFLPsolver& , solver )
{
    try
    {
        const double myZero = 1e-4;
        bool FoundViolatedCuts = false;
        const int NumNodes = getNnodes();
        if ( NumNodes > 1 || cutsAddedInCutCallback >= 100 )
        {
            return;
        }
        int errStat,
            n = solver.getNumFac(),
            m = solver.getNumCust();

        CPXENVptr CpxEnv; // Cplex environment, C-callable library style
        CpxEnv = CPXopenCPLEX(&errStat); // Open the environment


        int* weight     = new int[ m ];
        int* TDweight   = new int[ n ];
            for(int j=0; j<m; ++j) weight[j] = int( solver.getDemand( j ) );

        double* xlp = new double[m];
        double* rco = new double[m];
        double* ylp = new double[n];
        double* yrc = new double[n];

        VICcut* cutlist = NULL;
        VICcut* thiscut = NULL;
        VICcut* TDcut   = NULL;

        for(int i=0; i<n; ++i){ // look at each individual capacity constraint
            thiscut = NULL; // reset the thiscut to NULL
            // Retrive information on the y-variables
            ylp[i] =  getValue( solver.y[i]);
            TDweight[i] = solver.getCapacity( i );

            /* Use the ratio between the fixed opening cost and the capacity of the facility as proxy for reduced costs
             * Used for determining the lifting sequence in VIClci and VICkpsep
             */
            for(int i=0; i<n; ++i) yrc[i] = double( solver.getF( i ) ) / double ( solver.getCapacity( i ) );

            if ( ylp[i]>=myZero ){ // Check if we have a positive y-variable
                for(int j=0; j<m; ++j)
                { // Copy data into appropriate data structure
                    xlp[j] = getValue ( solver.x[i][j] );

                    /* Use the ratio between the cost of assigning customer j to facility i and the demand of customer j as proxy for reduced costs
                     * Used for determining the lifting sequence in VIClci and VICkpsep
                     */
                    rco[j] = double ( solver.getC(i , j ) ) / double ( solver.getDemand( j ) );
                }

                VIClci(m, solver.getCapacity( i ), 'L', weight, NULL, xlp, rco,&thiscut);
                if(NULL == thiscut)
                {
                    VICecikl(m, solver.getCapacity( i ), 0, 'L', weight, NULL, xlp, &thiscut);

                    if(NULL == thiscut)
                    {
                        VICkpsep ( CpxEnv , 0 , m , solver.getCapacity( i ) , 'L' , weight , NULL , xlp , rco , &thiscut );
                    }
                }


                /*
                 * If we have found no cuts try to separate from conv { (x,y) : sum_{j\in J} d_j ( x[i][j] + x[ii][j] ) <= s[i]*y[i] + s[ii]*y[ii],  x[i][j] + x[ii][j] <= 1, all binary }
                 */

                if ( !thiscut ) // We can only form pairs when i < n - 1
                {
                    if ( i < n -1 )
                    {

                    int cap;
                    int* w = new int[ m + 2 ];
                    double* z = new double[ m + 2 ];
                    double* rc = new double[ m + 2 ];

                    for ( int j = 0; j < m; ++j )
                    {
                        w [ j ] = weight [ j ];
                    }
                    //std::cout << "Trying to generate a double knapsack cut\n";
                    /*
                     * Fill in the information of the x variables
                     */
                    for ( int ii = i+1; ii < n; ++ii )
                    {
                        if ( getValue ( solver.y[ii] ) >= myZero )
                        thiscut = NULL;
                        cap = solver.getCapacity( i )  + solver.getCapacity( ii );
                        for ( int j = 0; j < m; ++j )
                        {
                            z [ j ] = getValue ( solver.x [ i ][ j ] ) + getValue( solver.x [ ii ][ j ] );
                            rc[ j ] = double ( solver.getC( i , j ) + solver.getC( ii , j ) ) / double( w[ j ] );
                        }
                        /*
                         * Fill in the information for the y variables
                         */
                        z[ m ] = 1.0 - getValue( solver.y[ i ] );
                        z[ m + 1 ] = 1.0 - getValue( solver.y[ ii ] );
                        rc [ m ] = double ( solver.getF( i ) ) / double ( solver.getCapacity( i ) );
                        rc [ m+1 ] = double ( solver.getF ( ii ) ) / double ( solver.getCapacity( ii ) );
                        w[ m ] = solver.getCapacity( i );
                        w [ m + 1 ] = solver.getCapacity( ii );

                        //VICkpsep( CpxEnv , 0 , m + 2 , cap ,'L' , w , NULL , z , rc , &thiscut );
                        VIClci( m +2 , cap , 'L' , w, NULL, z , rc, &thiscut );
                        if ( thiscut )
                        {
                            thiscut->UsrIVal = i;
                            thiscut->UsrRVal = ii;
                            //std::cout << "Adding cuts for i1 = " << i << " and i2 = " << ii << std::endl;
                            VICaddtolst( &cutlist, thiscut );
                            break;
                        }
                    }


                    delete [] z;
                    delete [] w;
                    delete [] rc;
                    }
                }
                else
                { // If we generated a cut, add this cut to the cutlist
                    thiscut->UsrIVal = i;
                    thiscut->UsrRVal = -1.0;
                    VICaddtolst(&cutlist, thiscut);
                }
            }

        }

        // Consider the total demand constraint
        int TD = solver.getTotalDemand ( );
        VIClci(n, TD, 'G', TDweight, NULL, ylp, yrc, &TDcut);
        if( !TDcut )
            VICecikl(n,TD, 0, 'G', TDweight, NULL, ylp, &TDcut);
        if( !TDcut )
            VICkpsep(CpxEnv, 0, n, TD, 'G', TDweight, NULL, ylp, yrc, &TDcut);
        if( TDcut ){
            TDcut->UsrIVal = -1;
            VICaddtolst(&cutlist, TDcut);
        }

        // Now add all the cuts we have generated
        int NumOfCuts = 0;
        // Now add all the cuts we have generated
        for( VICcut* cut = cutlist; cut!=NULL; cut=cut->nextcut ){
            NumOfCuts ++;
            cutsAddedInCutCallback++;
            IloExpr aCut = IloExpr( getEnv ( ) );
            int idx = cut->UsrIVal;
            int dIndx = cut->UsrRVal;

            if ( dIndx <= 0.5 )
            {   // We have a cut for a single knapsack structure
                if( idx>=0 ){ // We have a cut for the capacity of facility idx
                    for(int nz = 0; nz < cut->nzcnt; ++nz)
                    {
                        aCut+=(cut->nzval[nz])*solver.x[idx][cut->nzind[nz]];
                    }
                    add( aCut <= (cut->rhs)*solver.y[idx] );
                }else{
                    for(int nz = 0; nz < cut->nzcnt; ++nz)
                    {
                        aCut+= (cut->nzval[nz])*solver.y[cut->nzind[nz]];
                    }
                    add( aCut <= (cut->rhs) );
                }
            }
            else //if ( false )
            { // We have a cut for a pair of knapsacks + conflict graph
                doubleKPcutsAdded++;
                for ( int nz = 0; nz <cut->nzcnt; ++nz )
                {
                    if ( cut->nzind [ nz ] < m )
                    { // We are considering x variables
                        aCut += ( cut->nzval [ nz ] ) * ( solver.x[ cut->UsrIVal ][ cut->nzind [ nz ] ] + solver.x[ int( cut->UsrRVal ) ][ cut->nzind [ nz ] ] );
                    }
                    else if ( cut->nzind [ nz ] == m )
                    { //We are constidering y[idx]
                        aCut += ( cut->nzval [ nz ] ) *(1 -  solver.y[ idx ]);

                    }
                    else
                    { // We are considering the variable y[ dIndx ]
                        aCut += ( cut->nzval [ nz ] ) * (1-  solver.y[ dIndx ]);
                    }
                }
                add ( aCut <= cut->rhs );
            }


            aCut.end();
        }


        VICfreelst(&cutlist);


        errStat = CPXcloseCPLEX(&CpxEnv); // Close the cplex environment, and free the memory

        // Delete memory allocated to arrays
        delete[] weight;
        delete[] TDweight;
        delete[] xlp;
        delete[] rco;
        delete[] ylp;
        delete[] yrc;
        if ( NumOfCuts < 2 ) abortCutLoop() ;
        return;
    }
    catch ( IloException ie )
    {
        std::cerr << "IloException in user cut callback KnapsackSep : " << ie.getMessage ( ) << std::endl;
    }
    catch ( std::exception &e )
    {
        std::cerr << "Exception in user cut callback KnapsackSep : " << e.what ( ) << std::endl;
    }
    catch ( ... )
    {
        std::cerr << "Unknown exception in user cut callback KnapsackSep\n";
    }
}


SSCFLPsolver::SSCFLPsolver() :
    n(0),
    m(0),
    TD(0),
    hasCleanedUp(0),
    hasRun(0),
    displayStats(1),
    CplexOutOff(1),
    debugMe(0),
    solveSparseUsingDualAscentAlg(0),
    ObjVal(INT_MAX),
    myZero(1.0E-9),
    myOne(0.9999),
    myEpsilon(1.0E-4)
{
    try{
        // Initializing the sparse model and its cplex environment
        SparseModel = IloModel(env);
        SparseCplex = IloCplex(SparseModel);
        // Initializing the dense model and its cplex environment
        DenseModel  = IloModel(env);
        DenseCplex  = IloCplex(DenseModel);
        // Initializing the variables
        y = IloNumVarArray(env);
        x = IloVarMatrix(env);
        // Initializing the data containers
        c = IloNumMatrix(env);
        f = IloNumArray(env);
        d = IloNumArray(env);
        s = IloNumArray(env);
        duals = IloNumArray ( env );
        AssCst = IloRangeArray ( env );
        // Initializing pointers
        ySol = xSol = nullptr;
        NumCuts = new int[5];
        NumCuts[0] = NumCuts[1] = NumCuts[2] = NumCuts[3] = NumCuts[4] = 0;
        StartTime = double(clock())/double(CLOCKS_PER_SEC);

        // Initialize the test statistics
        stats = new testStats;
        stats->numberOfIterations =
        stats->time =
        stats->itWhereOptWasFound =
        stats->initialUpperBound =
        stats->bestUpperBound =
        stats->bestLowerBound =
        stats->percentageGap =
        stats->CuttingTime =
        stats->CutAndSolveTime =
        stats->avgNumOfDualItPerSparse =
        stats->numOfSparseSolved,
        stats->avgPercentLeft = 0;

    }catch(std::exception &e){
        std::cout << "Exception in the constructor of the SSCFLPsolver class: " << e.what() << std::endl;
        return;
    }catch(IloException &ie){
        std::cerr << "IloException in the constructor of the SSCFLPsolver class: " << ie.getMessage() << std::endl;
        return;
    }catch(...){std::cerr << "Exception originates from the constructor of the SSCFLPsolver class!\n"; throw;}// If you don't know the exception rethrow and hope for better times
}

/*********************************************************************************************************/
SSCFLPsolver::~SSCFLPsolver(){
    try{
        if(!hasCleanedUp) CleanUp();
        return;
    }catch(std::exception &e){
        std::cerr << "Exception in in the destructor of the SSCFLPsolver class: " << e.what() << std::endl;
        exit(1);
    }catch(IloException &ie){
        std::cerr << "IloException in the destructor of the SSCFLPsolver class: " << ie.getMessage() << std::endl;
        exit(1);
    }catch(...){std::cerr << "Exception originates from the destructor of the SSCFLPsolver class!\n"; throw;}// If you don't know the exception rethrow and hope for better times
}

/*********************************************************************************************************/
void SSCFLPsolver::BuildModels(){
    try{
        IloExpr expr(env);
        IloExpr TDexpr(env);
        // Build the objective function
        for(IloInt i=0; i<n; ++i) expr+=IloScalProd(x[i],c[i]);
        expr+=IloScalProd(y,f);
        // Add the objective function to the two models
        DenseObj = IloObjective ( env );
        SparseObj = IloObjective ( env );
        DenseObj = IloMinimize ( env, expr );
        SparseObj = IloMinimize ( env, expr );


        SparseModel.add ( SparseObj );
        DenseModel.add ( DenseObj );

        // Build and add the assignment constraints
        expr.clear();
        for(IloInt j=0; j<m; ++j)
        {
            for(IloInt i=0; i<n; ++i) expr+=x[i][j];

            AssCst.add ( expr == 1 );
            SparseModel.add( AssCst[j] );

            duals.add ( 0.0 );
            DenseModel.add( AssCst[j] );
            expr.clear();
            TD+=d[j];
        }


        // Build and add the capacity constraints
        for(IloInt i=0; i<n; ++i){
            for(IloInt j=0; j<m; ++j) expr+=d[j]*x[i][j];
            SparseModel.add(expr<=s[i]*y[i]);
            DenseModel.add(expr<=s[i]*y[i]);
            expr.clear();
            TDexpr+=s[i]*y[i];
        }
        for(IloInt i=0; i<n; ++i) for(IloInt j=0; j<m; ++j) SparseModel.add(x[i][j]-y[i]<=0);
        // Adding the total demand constraint
        SparseModel.add(TDexpr>=TD);
        DenseModel.add(TDexpr>=TD);

        if ( !SparseCplex.isMIP() ){
            SparseModel.add(IloConversion(env,y,ILOBOOL));
            for(IloInt i=0; i<n; ++i) SparseModel.add(IloConversion(env, x[i],ILOBOOL));
        }

    }catch(std::exception &e){
        std::cerr << "Exception in BuildModels in the SSCFLPsolver class: " << e.what() << std::endl;
        exit(1);
    }catch(IloException &ie){
        std::cerr << "IloException in BuildModels in the SSCFLPsolver class: " << ie.getMessage() << std::endl;
        exit(1);
    }catch(...){std::cerr << "Exception originates from the BuildModels function in the SSCFLPsolver class!\n"; throw;}// If you don't know the exception rethrow and hope for better times
}

/*********************************************************************************************************/
void SSCFLPsolver::PrintProgramInfo(){
    if( !CplexOutOff ){
        std::cout << "\n\n";
        std::cout << "||==============================================================================================||\n";
        std::cout << "|| --> General information about the SSCFLP solver <--                                          ||\n";
        std::cout << "||==============================================================================================||\n";
        std::cout << " | Output from cplex has been redirected to the null stream                                     |\n";
        std::cout << " | Cutting planes of the following type is separated at the root node                           |\n";
        std::cout << " |    Lifted cover inequalities (heuristic separation)                                          |\n";
        std::cout << " |    Extended cover inequalities (heuristic separation)                                        |\n";
        std::cout << " |    Fenchel cutting planes (exact separation)                                                 |\n";
        std::cout << " |    Variable upper bounds (exact separation)                                                  |\n";
        std::cout << " | An initial solution is generated using a local branching heuristic                           |\n";
        std::cout << " | A potential optimality gap is closed by a cut and solve algorithm                            |\n";
        std::cout << "||==============================================================================================||\n";
    }
}

/*********************************************************************************************************/
bool SSCFLPsolver::SepCuts(){
    try{
        bool CutsSeparated = false;
        int errStat;

        CPXENVptr CpxEnv; // Cplex environment, C-callable library style
        CpxEnv = CPXopenCPLEX(&errStat); // Open the environment

        int* weight     = new int[m];
        int* TDweight   = new int[n];
            for(int j=0; j<m; ++j) weight[j] = int(d[j]);

        double* xlp = new double[m];
        double* rco = new double[m];
        double* ylp = new double[n];
        double* yrc = new double[n];

        std::vector<std::pair<int,int>> VariableUBs;

        //VICcut* cutlist = NULL;
        VICcut* thiscut = NULL;
        VICcut* TDcut   = NULL;
        VICcut* cutlist = NULL;


        for(int i=0; i<n; ++i){ // look at each individual capacity constraint
            thiscut = NULL; // reset the thiscut to NULL
            // Retrive information on the y-variables
            ylp[i] = DenseCplex.getValue(y[i]);
            TDweight[i] = s[i];
            /* If DenseCplex has not yet been converted to a mip, get the absolute value of the reduced costs
             * else use the ratio between the fixed opening cost and the capacity of the facility
             * Used for determining the lifting sequence in VIClci and VICkpsep
             */
            if(!DenseCplex.isMIP()) for(int i=0; i<n; ++i) yrc[i] = std::abs(DenseCplex.getReducedCost(y[i]));
            else for(int i=0; i<n; ++i) yrc[i] = double(f[i])/double(s[i]);

            if(ylp[i]>myZero){ // Check if we have a positive y-variable
                for(int j=0; j<m; ++j)
                { // Copy data into appropriate data structure
                    xlp[j] = DenseCplex.getValue(x[i][j]);
                    /* If DenseCplex has not yet been converted to a mip, get the absolute value of the reduced costs
                     * else use the ratio between the cost of assigning customer j to facility i and the demand of customer j
                     * Used for determining the lifting sequence in VIClci and VICkpsep
                     */
                    if(! DenseCplex.isMIP()) rco[j] = std::abs(DenseCplex.getReducedCost(x[i][j]));
                    else rco[j] = double(c[i][j])/double(d[j]);
                    if ( xlp[j] >= ylp[i] + 0.01 ) VariableUBs.push_back ( std::pair<int,int>(i,j) );
                }

                VIClci(m, s[i], 'L', weight, NULL, xlp, rco,&thiscut);
                if(NULL == thiscut)
                {
                    VICecikl(m, s[i], 0, 'L', weight, NULL, xlp, &thiscut);

                    if(NULL == thiscut)
                    {
                        VICkpsep(CpxEnv, 0, m,s[i],'L',weight,NULL,xlp, rco,&thiscut);
                        if ( thiscut ) ++ NumCuts[2];
                    }
                    else ++ NumCuts[1];
                }
                else ++NumCuts[0];


                if ( thiscut )
                { // If we generated a cut, add this cut to the cutlist
                    thiscut->UsrIVal = i;
                    thiscut->UsrRVal = -1.0;
                    VICaddtolst(&cutlist, thiscut);
                    CutsSeparated = true;
                }
            }

        }



        if ( !CutsSeparated )
        {
            int cap;
            int* w = new int[ m + 2 ];
            double* z = new double[ m + 2 ];
            double* rc = new double[ m + 2 ];

            for ( int j = 0; j < m; ++j )
            {
                w [ j ] = weight [ j ];
            }

            for ( int i = 0; !CutsSeparated && ( i < n ); ++i ) //!CutsSeparated &&
            {
                if ( ( ylp[i] >= 0.01 ) && ( i < n -1 ) ) // We can only form pairs when i < n - 1
                {
                    for ( int ii = i+1; !CutsSeparated && ( ii < n ); ++ii ) //!CutsSeparated &&
                    {
                        if ( true ) //ylp[ii] >= 0.01
                        {
                            thiscut = NULL;
                            cap = s [ i ] + s [ ii ];
                            for ( int j = 0; j < m; ++j )
                            {
                                z [ j ] = DenseCplex.getValue ( x [ i ][ j ] ) + DenseCplex.getValue( x [ ii ][ j ] );
                                if ( !DenseCplex.isMIP() ) rc [ j ] = DenseCplex.getReducedCost( x [ i ][ j ]  ) + DenseCplex.getReducedCost( x[ ii ][ j ] );
                                else rc[ j ] = double ( c[ i ][ j ] + c[ ii ][ j ] ) / double( d[ j ] );
                            }
                            /*
                             * Fill in the information for the y variables
                             */
                            z[ m ] = 1.0 - DenseCplex.getValue( y[ i ] );
                            z[ m + 1 ] = 1.0 - DenseCplex.getValue( y[ ii ] );
                            rc [ m ] = double ( f[ i ] ) / double ( s[ i ] );
                            rc [ m+1 ] = double ( f[ ii ] ) / double ( s[ ii ] );
                            w[ m ] = s[ i ];
                            w [ m + 1 ] = s[ ii ];


                            VIClci( m +2 , cap , 'L' , w, NULL, z , rc, &thiscut );
                            // If we did not separate using lifted cover inequality, try with a fenchel cut
                            if ( !thiscut )
                                VICkpsep( CpxEnv , 0 , m + 2 , cap ,'L' , w , NULL , z , rc , &thiscut );

                            if ( thiscut )
                            {
                                ++NumCuts[4];
                                CutsSeparated = true;
                                thiscut->UsrIVal = i;
                                thiscut->UsrRVal = ii;
                                //std::cout << "Adding cuts for i1 = " << i << " and i2 = " << ii << std::endl;
                                VICaddtolst( &cutlist, thiscut );
                                break;
                            }
                        }
                    }
                }
            }

            delete [] z;
            delete [] w;
            delete [] rc;
        }


        // Consider the total demand constraint
        VIClci(n, TD, 'G', TDweight, NULL, ylp, yrc, &TDcut);
        if( !TDcut )
            VICecikl(n,TD, 0, 'G', TDweight, NULL, ylp, &TDcut);
        if( !TDcut )
            VICkpsep(CpxEnv, 0, n, TD, 'G', TDweight, NULL, ylp, yrc, &TDcut);
        if( TDcut ){
            TDcut->UsrIVal = -1;
            VICaddtolst(&cutlist, TDcut);
        }

        // Now add all the cuts we have generated
        for( VICcut* cut = cutlist; cut!=NULL; cut=cut->nextcut ){

            IloExpr aCut = IloExpr(env);
            int idx = cut->UsrIVal;
            int dIndx = cut->UsrRVal;

            if ( dIndx <= 0.5 )
            {   // We have a cut for a single knapsack structure (capacity or total demand)
                if( idx>=0 )
                { // We have a cut for the capacity of facility idx
                    for(int nz = 0; nz < cut->nzcnt; ++nz)
                    {
                        aCut+=(cut->nzval[nz])*x[idx][cut->nzind[nz]];
                    }
                    DenseModel.add( aCut <= (cut->rhs)*y[idx] );
                    SparseModel.add( aCut <= (cut->rhs)*y[idx] );
                }
                else
                { // This is a cut for the total demand constraints
                    for(int nz = 0; nz < cut->nzcnt; ++nz)
                    {
                        aCut+= (cut->nzval[nz])*y[cut->nzind[nz]];
                    }
                    DenseModel.add( aCut <= cut->rhs);
                    SparseModel.add( aCut <= (cut->rhs) );
                }
            }
            else //if ( tr )
            { // We have a cut for a pair of knapsacks + conflict graph
                for ( int nz = 0; nz <cut->nzcnt; ++nz )
                {
                    if ( cut->nzind [ nz ] < m )
                    { // We are considering x variables
                        aCut += ( cut->nzval [ nz ] ) * ( x[ cut->UsrIVal ][ cut->nzind [ nz ] ] + x[ int( cut->UsrRVal ) ][ cut->nzind [ nz ] ] );
                    }
                    else if ( cut->nzind [ nz ] == m )
                    { //We are constidering y[idx]
                        aCut += ( cut->nzval [ nz ] ) *(1 -  y[ idx ]);

                    }
                    else
                    { // We are considering the variable y[ dIndx ]
                        aCut += ( cut->nzval [ nz ] ) * (1-  y[ dIndx ]);
                    }
                }
                std::string CstName = "DKP_" + std::to_string ( doubleKPcutsAdded );

                IloRange constraint ( aCut <= cut->rhs );
                constraint.setName ( CstName.c_str() );
                DenseModel.add ( constraint );
                SparseModel.add ( aCut <= cut->rhs );
                ++doubleKPcutsAdded;
                //constraint.end ( );
            }


            aCut.end();
        }

        //DenseCplex.exportModel( "AfterCuts.LP" );


        for ( auto it = VariableUBs.begin ( ); it!=VariableUBs.end ( ); ++it )
        {
            CutsSeparated = true;
            DenseModel.add ( x[it->first][it->second] - y[it->first] <= 0);
            SparseModel.add ( x[it->first][it->second] - y[it->first] <= 0);
            ++NumCuts[3];
        }

        errStat = CPXcloseCPLEX(&CpxEnv); // Close the cplex environment, and free the memory
        delete[] weight; // Delete memory allocated to arrays
        delete[] xlp;
        delete[] rco;
        delete[] ylp;

        return CutsSeparated;
    }catch(std::exception &e){
        std::cerr << "Exception in SepCuts in the SSCFLPsolver class: " << e.what() << std::endl;
        exit(1);
    }
}

/*********************************************************************************************************/
void SSCFLPsolver::CuttingPhase(){
    try{
        if ( !CplexOutOff )
        {
            std::cout << "========== RUNNING CUTTING PLANE ALGORITHM ==========\n";
        }
        bool CONTINUE = true;
        std::vector<double> LBs;
        std::stringstream err;

        auto CutStart = CPUclock::now ( );
        while( CONTINUE ){
            if( DenseCplex.solve() ){
                for ( int j=0; j<m; ++j ) DenseCplex.getDuals( duals, AssCst );
                LBs.push_back(DenseCplex.getObjValue());
                CONTINUE = SepCuts();
            }
            else{ err << "Could not solve the Dense Problem. Status is: " << DenseCplex.getStatus() << std::endl; throw std::runtime_error(err.str());} // If not solve, throw error

            if(LBs.size() > 5){
                if(LBs[LBs.size()-1] - LBs[LBs.size()-2] <= 0.10 ) CONTINUE = false;
            }
            std::cout << "Lower bound : " << LBs.back ( ) << std::endl;
        }
        CutLowerBound = LBs.back ( );

        stats->WeakLowerBound = *LBs.begin ( );
        stats->LowerBoundAfterCusts = LBs.back ( );

        auto CutEnd = CPUclock::now ( );
        CutTime = duration_cast<duration<double>>( CutEnd - CutStart ).count ( );


        if ( true )
        {
            std::cout << "Initial lower bound      : " << *LBs.begin ( ) << std::endl;
            std::cout << "Lower bound after cuts   : " << LBs.back ( ) << std::endl;
            std::cout << "Number of cuts generated :\n";
            std::cout << "  LCI                    : " << NumCuts[0] << "\n";
            std::cout << "  ECI                    : " << NumCuts[1] << "\n";
            std::cout << "  FCP                    : " << NumCuts[2] << "\n";
            std::cout << "  VUB                    : " << NumCuts[3] << "\n";
            std::cout << "  ECP                    : " << NumCuts[4] << "\n";
            std::cout << "Total time on cutting    : " << CutTime << std::endl;
            std::cout << "=====================================================\n";
        }

    }catch(std::exception &e){
        std::cerr << "Exception in CuttingPhase in the SSCFLPsolver class: " << e.what() << std::endl;
        exit(1);
    }catch(IloException &ie){
        std::cerr << "IloException in CuttingPhase in the SSCFLPsolver class: " << ie.getMessage() << std::endl;
        exit(1);
    }
}

/*********************************************************************************************************/
void SSCFLPsolver::CheckIfFix(std::vector<int> &ones){
    try{
        int CapSum;
        std::vector<int> FixedOnes;
        std::vector<int>::iterator it;

        for(it=ones.begin(), CapSum = 0; it!=ones.end(); ++it) CapSum+=s[*it];

        it = ones.begin();
        while(it!=ones.end()){
            if(CapSum - s[*it] >= TD) ones.erase(it);
            else ++it;
        }

    }catch(std::exception &e){
        std::cerr << "Exception in CheckIfFix in the SSCFLPsolver class: " << e.what() << std::endl;
        exit(1);
    }catch(IloException &ie){
        std::cerr << "IloException in CheckIfFix in the SSCFLPsolver class: " << ie.getMessage() << std::endl;
        exit(1);
    }
}

/*********************************************************************************************************/
void SSCFLPsolver::ChangeDenseToSemiLagrangean ( )
{
    try
    {
        double maxDual = 0;
        int NewCoef = 0;
        for ( int j=0; j<m; ++j )
        {
            if ( duals[j] > maxDual ) maxDual = duals[j];
            AssCst[j].setBounds ( 0.0 , 1.0 );
        }

        for ( int i=0; i<n; ++i )
        {
            for ( int j=0; j<m; ++j )
            {
                NewCoef = int ( c[i][j] - maxDual + 0.5 );
                DenseObj.setLinearCoef( x[i][j] , NewCoef );
            }
        }
    }
    catch ( IloException &ie )
    {
        std::cerr << "IloException in ChangeDenseToSemiLagragnean in the SSCFLPsolver class : " << ie.getMessage ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( std::exception &e )
    {
        std::cerr << "Exception in ChangeDenseToSemiLagragnean in the SSCFLPsolver class : " << e.what ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
}

/*********************************************************************************************************/
void SSCFLPsolver::CleanUp(){
    try{
        // Freeing memory that has been new'ed
        if ( ySol ) delete[] ySol;
        if ( xSol ) delete[] xSol;

        // Freeing memory occupied by variables
        if(y.getImpl()) y.end();

        if(x.getImpl()){
            if(n>0) for(IloInt i=0; i<n; ++i) x[i].end();
            x.end();
        }

        if ( NumCuts ) delete[] NumCuts;
        // Freeing memory occupied by data
        if(c.getImpl()){
            if(n>0) for(IloInt i=0; i<n; ++i) c[i].end();
            c.end();
        }

        if(f.getImpl()) f.end();

        if(d.getImpl()) d.end();

        if(s.getImpl()) s.end();


        // Freeing memory occupied by models and cplex environments
        //if ( SparseCplex.getImpl() )SparseCplex.end();
        //if ( SparseModel.getImpl() )SparseModel.end();


        //if ( DenseCplex.getImpl() )DenseCplex.end();
        //if ( DenseModel.getImpl() )DenseModel.end();


        // Freeing the memory consumed by IloEnv
        //if ( env.getImpl() ) env.end();
        hasCleanedUp = true;
    }catch(...){std::cerr << "Exception originates from the CleanUp function in the SSCFLPsolver class!\n"; throw;} // If you don't know the exception rethrow and hope for better times
}

/*********************************************************************************************************/
void SSCFLPsolver::Load(char* Filename, int format){
    try{
        double anyNum;


        maxF = maxC = maxD = maxS = avgF = avgC = avgD = avgS = 0;
        minF = minC = minD = minS = INT_MAX;

        std::stringstream err;

        std::ifstream data_file(Filename); // Opening the data file
        if(!data_file) throw std::runtime_error(std::string("Could not open the file: ") + Filename + std::string("\n\n I will terminate now!\n\n"));
        else std::cout << " | Solving problem specified by the file: " << Filename << std::endl;

        if(1 == format){
            data_file >> m;
            if(!data_file) throw std::runtime_error("Could not read the number of customers. Termintaing!\n");

            data_file >> n;
            if(!data_file) throw std::runtime_error("Could not read the number of facilities. Termintaing!\n");

            if(n<=0 || m<=0){
                err << "Non positve values for number of facilities and customers. Provided was n=" << n << ", and m=" << m << ". Terminating\n";
                throw std::runtime_error(err.str().c_str());
            }
            for(IloInt i=0; i<n; ++i) c.add(IloNumArray(env,m));
            for(IloInt j=0; j<m; ++j){ // Reading the assingment costs
                for(IloInt i=0; i<n; ++i){
                    data_file >> anyNum;
                    if(!data_file || anyNum < 0){
                        err << "Could not read the assignment cost of index (" << i << "," << j << "j). Terminating!\n";
                        throw std::runtime_error(err.str().c_str());
                    }
                    c[i][j]=int(anyNum);

                    // Gathering stats
                    if(true){
                        if(c[i][j]>maxC) maxC = c[i][j];
                        else if(c[i][j]<minC) minC = c[i][j];
                        avgC+=c[i][j];
                    }
                }
            }

            for(IloInt j=0; j<m; ++j){ // Reading the demands
                data_file >> anyNum;
                if(!data_file || anyNum<=0){
                    err << "Could not read the demand of customer " << j << ". Terminating!\n";
                    throw std::runtime_error(err.str().c_str());
                }
                d.add( int(anyNum) );

                // Gathering stats
                if(true){
                    if(d[j]>maxD)       maxD = d[j];
                    else if(d[j]<minD)  minD = d[j];
                    avgD+=d[j];
                }
            }

            for(IloInt i=0; i<n; ++i){ // Reading the fixed opening costs
                data_file >> anyNum;
                if(!data_file || anyNum<0 ){
                    err << "Could not read the fixed opening cost of facility " << i << ". Terminating!\n";
                    throw std::runtime_error(err.str().c_str());
                }
                f.add( int(anyNum) );

                // Gathering stats
                if(true){
                    if(f[i]>maxF)       maxF = f[i];
                    else if(f[i]<minF)  minF = f[i];
                    avgF+=f[i];
                }
                x.add(IloNumVarArray(env, m, 0, 1, ILOFLOAT));
                for ( int j=0; j<m; ++j )
                {
                    std::string xName = "x_" + std::to_string ( i ) + "_" + std::to_string ( j );
                    x[i][j].setName( xName.c_str ( ) );
                }
                std::string yName = "y_" + std::to_string ( i );
                y.add ( IloNumVar ( env , 0 , 1 , ILOFLOAT , yName.c_str ( ) ) );
            }

            for(IloInt i=0; i<n; ++i){
                data_file >> anyNum;
                if(!data_file || anyNum<=0 ){
                    err << "Could not read the capacity of facility " << i << ". Terminating!\n";
                    throw std::runtime_error(err.str().c_str());
                }
                s.add( int(anyNum) );

                // Gathering stats
                if(true){
                    if(s[i]>maxS)       maxS = s[i];
                    else if(s[i]<minS)  minS = s[i];
                    avgS += s[i];
                }
            }

            if(true){
                avgF = avgF / ( double(n)   );
                avgC = avgC / ( double(m*n) );
                avgD = avgD / ( double(m)   );
                avgS = avgS / ( double(n)   );
            }
        }else if(2 == format || 3 == format || 4 == format){ // Instances of the Holmberg, Yang or Gadegaard type

            data_file >> n;
            if(!data_file) throw std::runtime_error("Could not read the number of customers. Termintaing!\n");

            data_file >> m;
            if(!data_file) throw std::runtime_error("Could not read the number of facilities. Termintaing!\n");

            if(n<=0 || m<=0)
            {
                err << "Non positve values for number of facilities and customers. Provided was n=" << n << ", and m=" << m << ". Terminating\n";
                throw std::runtime_error(err.str().c_str());
            }

            for ( int i=0; i<n; ++i )
            {   // Reading the capacities and the demands
                data_file >> anyNum;
                if(!data_file || anyNum<=0 ){
                    err << "Could not read the capacity of facility " << i << ". Terminating!\n";
                    throw std::runtime_error(err.str().c_str());
                }
                s.add( int(anyNum) );

                data_file >> anyNum;
                if(!data_file || anyNum<0 ){
                    err << "Could not read the fixed opening cost of facility " << i << ". Terminating!\n";
                    throw std::runtime_error(err.str().c_str());
                }
                f.add( int(anyNum) );
            }

            for(IloInt j=0; j<m; ++j)
            { // Reading the demands
                data_file >> anyNum;
                if(!data_file || anyNum<=0){
                    err << "Could not read the demand of customer " << j << ". Terminating!\n";
                    throw std::runtime_error(err.str().c_str());
                }
                d.add( int(anyNum) );
            }
            avgC = 0;
            for ( int i=0; i<n; ++i )
            {
                c.add ( IloNumArray ( env , m ) );
                x.add(IloNumVarArray(env, m, 0, 1, ILOFLOAT));
                y.add(IloNumVar(env,0,1,ILOFLOAT));
                for ( int j=0; j<m; ++j )
                {
                    data_file >> anyNum;
                    if(!data_file || anyNum < 0){
                        err << "Could not read the assignment cost of index (" << i << "," << j << "j). Terminating!\n";
                        throw std::runtime_error(err.str().c_str());
                    }
                    c[i][j]=int(anyNum);
                    avgC += c[i][j];
                }
            }
            avgC = avgC / double ( m*n );
        }else{
            err << "You specified the format " << format << " which does not exists. Terminating!\n";
            throw std::runtime_error(err.str().c_str());
        }
        stats->n = n;
        stats->m = m;
        BuildModels();

    }catch(std::exception &e){
        std::cerr << "Exception in Load(Filename) in the SSCFLPsolver class: " << e.what() << std::endl;
        exit(1);
    }catch(IloException &ie){
        std::cerr << "IloException in Load(Filename) in the SSCFLPsolver class: " << ie.getMessage() << std::endl;
        exit(1);
    }catch(...){std::cerr << "Exception originates from the Load(Filename) function in the SSCFLPsolver class!\n"; throw;} // If you don't know the exception rethrow and hope for better times
}

/*********************************************************************************************************/
void SSCFLPsolver::Load(int nn, int mm, int** cc, int* ff, int* dd, int* ss){
    try{
        std::cout << " | Data is loaded as pointers to integer arrays" << std::setw(66) << "| \n";
        std::stringstream err;

        n = nn;
        m = mm;
        if(n<=0 || m<=0){
            err << "Non positve values for number of facilities and customers. Provided was n=" << n << ", and m=" << m << ". Terminating\n";
            throw std::runtime_error(err.str().c_str());
        }
        if(cc){ // Check if pointing to something useful
            for(IloInt i=0; i<n; ++i){ // if
                if(cc[i]){ // Check if pointing to something useful
                    c.add(IloNumArray(env,m));
                    for(IloInt j=0; j<m; ++j) c[i][j] = cc[i][j]; // copy the data for the cost matrix
                }
                else{
                    err << "The pointer to the " << i << "'th row of the cost matrix seems to point to null. Terminating the process...\n" << "c[" << i<<"]=" << c[i] << std::endl;
                    throw std::runtime_error(err.str());
                }
            }
        }else{
            throw std::runtime_error("The pointer to the cost matrix seems to point to NULL. Terminating the process...\n");
        }

        if(ff){ // Check if pointing to something useful
            for(IloInt i=0; i<n; ++i) f.add(ff[i]); // copy data
        }else{
            throw std::runtime_error("The pointer to the fixed costs is pointing to NULL. Terminating the process...\n");
        }

        if(dd){ // Check if pointing to something useful
            for(IloInt j=0; j<m; ++j) d.add(dd[j]); // copy data
        }else{
            throw std::runtime_error("The pointer to the customer demands seems to point to NULL. Terminating the process...\n");
        }

        if(ss){// Check if pointing to something useful
            for(IloInt i=0; i<n; ++i) s.add(ss[i]); //copy data
        }else{
            throw std::runtime_error("The pointer to the facility capacities seems to point to NULL. Terminating the process...\n");
        }
        stats->n = n;
        stats->m = m;

        BuildModels();
    }catch(std::exception &e){
        std::cerr << "Exception in Load(data-pointers) in the SSCFLPsolver class: " << e.what() << std::endl;
        exit(1);
    }catch(IloException &ie){
        std::cerr << "IloException in Load(data-pointers) in the SSCFLPsolver class: " << ie.getMessage() << std::endl;
        exit(1);
    }
}

/*********************************************************************************************************/
void SSCFLPsolver::Load ( rdDat* data )
{
    try
    {
        n = data->getNumFac();
        m = data->getNumCust();

        for ( int i = 0; i < n; ++i )
        {
            f.add ( data->getF ( i ) ); // Copy fixed costs
            s.add ( data->getS ( i ) ); // Copy capacities
            c.add ( IloNumArray ( env, m) ); // allocate memory for assignment costs
            for ( int j = 0; j < m; ++j )
            {
                c[i][j] = data->getC ( i , j ); // Copy assignment costs
            }
        }
        for ( int j = 0; j < m; ++j)
        {
            d.add ( data->getD ( j ) ); // Copy the demands
        }

        stats->n = n;
        stats->m = m;

        /*
         * Build the model
         */

        for ( IloInt i = 0; i < n; ++i )
        { // Creating variables with nane
                x.add(IloNumVarArray(env, m, 0, 1, ILOFLOAT));
                for ( int j=0; j<m; ++j )
                {
                    std::string xName = "x_" + std::to_string ( i ) + "_" + std::to_string ( j );
                    x[i][j].setName( xName.c_str ( ) );
                }
                std::string yName = "y_" + std::to_string ( i );
                y.add ( IloNumVar ( env , 0 , 1 , ILOFLOAT , yName.c_str ( ) ) );
            }
        BuildModels();

    }
    catch(std::exception &e)
    {
        std::cerr << "Exception in Load(rdDat) in the SSCFLPsolver class: " << e.what() << std::endl;
        exit(1);
    }
    catch(IloException &ie)
    {
        std::cerr << "IloException in Load(rdDat) in the SSCFLPsolver class: " << ie.getMessage() << std::endl;
        exit(1);
    }
}
/*********************************************************************************************************/
bool SSCFLPsolver::setSolution(int nn, int mm, int* yval, int* xass,int UBval){
    IloNumVarArray vars = IloNumVarArray(env);
    IloNumArray vals = IloNumArray(env);
    int* CapSums =  new int[n]();

    try{
        for(int j=0; j<m; ++j)
        {
            vars.add(x[xass[j]][j]);
            vals.add(1.0);
            CapSums[xass[j]]+=d[j];
        }
        for(int i=0; i<n; ++i)
        {
            if( CapSums[i] >= s[i]+1) throw std::runtime_error("The solution provided was not feasible\n");
            vars.add(y[i]);
            vals.add(yval[i]);
        }
        // Add the mip start
        SparseCplex.addMIPStart(vars, vals, IloCplex::MIPStartSolveFixed);
        // Deallocate memory
        delete[] CapSums;
        vars.end();
        vals.end();
        return true;
    }catch(std::exception &e)
    {
        std::cerr << "Exception in setSolution in the SSCFLPsolver class: " << e.what() << std::endl << "No solution set!\n";
        delete[] CapSums;
        vars.end();
        vals.end();
        return false;
    }
    catch(IloException &ie)
    {
        std::cerr << "IloException in setSolution in the SSCFLPsolver class: " << ie.getMessage() << std::endl << "No solution set!\n";
        delete[] CapSums;
        vars.end();
        vals.end();
        return false;
    }
}

/*********************************************************************************************************/
void SSCFLPsolver::printStats(){

    maxC = maxF = maxD = maxS = avgF = avgC = avgD = avgS = 0;
    minC = minF = minD = minS = INT_MAX;
    for ( int i = 0; i < n; ++i )
    {
        avgF += f[ i ];
        if ( minF > f[ i ] ) minF = f[ i ];
        if ( maxF < f[ i ] ) maxF = f[ i ];

        avgS += s[i];
        if ( minS > s[ i ] ) minS = s[ i ];
        if ( maxS < s[ i ] ) maxS = s[ i ];

        for ( int j = 0; j < m; ++j )
        {
            avgC += c[ i ][ j ];
            if ( minC > c[ i ][ j ] ) minC = c[ i ][ j ];
            if ( maxC < c[ i ][ j ] ) maxC = c[ i ][ j ];

            if ( 0 == i )
            {
                avgD += d[ j ];
                if ( minD > d[ j ] ) minD = d[ j ];
                if ( maxD < d[ j ] ) maxD = d[ j ];
            }
        }
    }
    avgC = avgC / double( m * n );
    avgF = avgF / double ( n );
    avgS = avgS / double ( n );
    avgD = avgD / double ( m );

    std::cout << "||==============================================================================================||\n";
    std::cout << "|| --> Statistics for the loaded instance <-- \n";
    std::cout << "||==============================================================================================||\n";
    std::cout << " | No. of facilities:\t"  << n    << std::endl;
    std::cout << " | No. of customers:\t"   << m    << std::endl;
    std::cout << " | Fixed opening cost:"   << std::endl;;
    std::cout << " |\t Maximum:\t"          << maxF  << std::endl;
    std::cout << " |\t Minimum:\t"          << minF  << std::endl;
    std::cout << " |\t Average:\t"          << avgF  << std::endl;
    std::cout << " | Assignment costs:"     << std::endl;;
    std::cout << " |\t Maximum:\t"          << maxC  << std::endl;
    std::cout << " |\t Minimum:\t"          << minC  << std::endl;
    std::cout << " |\t Average:\t"          << avgC  << std::endl;
    std::cout << " | Customer demands:"     << std::endl;;
    std::cout << " |\t Maximum:\t"          << maxD  << std::endl;
    std::cout << " |\t Minimum:\t"          << minD  << std::endl;
    std::cout << " |\t Average:\t"          << avgD  << std::endl;
    std::cout << " | Facility capacities:"  << std::endl;
    std::cout << " |\t Maximum:\t"          << maxS  << std::endl;
    std::cout << " |\t Minimum:\t"          << minS  << std::endl;
    std::cout << " |\t Average:\t"          << avgS  << std::endl;
    std::cout << "||==============================================================================================||\n";
}

/*********************************************************************************************************/
bool SSCFLPsolver::Run(){
    try{
        /*
         * Print info about the solver
         */
        PrintProgramInfo();
        printStats ( );
        hasRun = true;
        std::vector<int> ONES;
        std::vector<int>::iterator iter;

        auto RuntimeStart = CPUclock::now ( );
        bool CONTINUE = true;
        int it=0;
        double DenseObjVal, endTime, denseTime, sparseTime;

        ySol = new int[n];
        ySol[0] = -1; // Used to indicate no solution has been found
        xSol = new int[m];

        SparseCplex.setParam(   IloCplex::EpAGap,       0.9);
        SparseCplex.setParam(   IloCplex::EpGap,        0.0);
        SparseCplex.setParam(   IloCplex::ParallelMode, 1);
        //SparseCplex.setParam(   IloCplex::MIPSearch,    1);

        DenseCplex.setParam(    IloCplex::EpGap,        0.0);
        DenseCplex.setParam(    IloCplex::EpAGap,       0.0);
        DenseCplex.setParam(    IloCplex::ParallelMode, 1);

        DenseCplex.setOut(env.getNullStream());
        DenseCplex.setWarning(env.getNullStream());
        //SparseCplex.setOut(env.getNullStream());
        SparseCplex.setParam ( IloCplex::MIPEmphasis , 0 );
        SparseCplex.setWarning ( env.getNullStream ( ) );

        /*
         * Starting to fix the variables that should be fixed after cutting phase
         */
        IloRangeArray Fixations ( env );
        for ( auto it = fixBeforeCuttingPhase.begin ( ); it != fixBeforeCuttingPhase.end ( ); ++it )
        {
            Fixations.add ( x[ it->first ][ it->second ] == 0 );
            DenseModel.add ( Fixations[ Fixations.getSize ( ) - 1 ] );
            SparseModel.add ( Fixations[ Fixations.getSize ( ) - 1 ] );
        }

        /*
         * Starting the cutting phase
         */
        auto CutStart = CPUclock::now ( );
        CuttingPhase();
        auto CutEnd =CPUclock::now ( );
        stats->CuttingTime = duration_cast<duration<double>>( CutEnd - CutStart ).count ( );

        DenseModel.add(IloConversion(env,y,ILOBOOL));
        /*
         * Starting to fix the variables that should be fixed after cutting phase
         */
        for ( auto it = fixAfterCuttingPhase.begin ( ); it != fixAfterCuttingPhase.end ( ); ++it )
        {
            Fixations.add ( x[ it->first ][ it->second ] == 0 );
            DenseModel.add ( Fixations[ Fixations.getSize ( ) - 1 ] );
            SparseModel.add ( Fixations[ Fixations.getSize ( ) - 1 ] );
        }

        // Tell SparseCplex to use the cut callback
        SparseCplex.use ( KnapsackSep( env, *this ) );

        /*
         * Run the local branching heuristic
         */
        //RunLBHeur ( );

        IloNumVarArray Zeros(env);
        IloNumArray lb(env);
        IloNumArray ub(env);

        solution incumbent = solution ( n , m);
        if ( ObjVal <= INT_MAX -10 ) incumbent.setSolution( ySol, xSol, ObjVal );
        std::cout << "\n\n";
        std::cout << std::setw(10) << std::left<< "It" << std::setw(10) << std::left << "Lb" << std::setw(10) << std::left << "ub" << std::setw(10) << std::left << "gap" << std::setw(10) << std::left << "D-time" << std::setw(10) << std::left << "S-time" << std::setw(10) << std::left << "time" << std::endl;
        std::cout << "============================================================================\n"<< std::setw(0) << std::left;
        auto CnSStart = CPUclock::now ( );
        stats->initialUpperBound = ObjVal;
        while(CONTINUE){
            cutsAddedInCutCallback = 0;
            ++it;
            stats->numberOfIterations = it;
            ONES.clear();
            auto DenseStart = CPUclock::now ( );
            DenseCplex.setParam(IloCplex::CutUp, ObjVal-0.9);
            if(DenseCplex.solve())  DenseObjVal = DenseCplex.getObjValue();
            else{ // We are optimal as we cannot find a solution improving the incumbent in the feasible set represented by the dense problem
                DenseObjVal = roundToTwo(double(ObjVal));
            } //throw std::runtime_error("Cplex could not solve the dense problem. Cplex status is: " + DenseCplex.getStatus() );
            auto DenseEnd = CPUclock::now ( );
            denseTime = duration_cast<duration<double>>( DenseEnd - DenseStart ).count ( );

            if(ObjVal - DenseObjVal <=0.9 ){ std::cout << it << "\t" << DenseObjVal << std::endl; break;} // The ObjVal Solution is optimal
            IloExpr piercingcut(env);
            for(IloInt i=0; i<n; ++i){
                if( DenseCplex.getValue(y[i]) <= 0.0001 ){
                    Zeros.add(y[i]);
                    for(IloInt j=0; j<m; ++j){ Zeros.add(x[i][j]); lb.add(0); ub.add(1);}
                    piercingcut+=y[i];
                    lb.add(0);
                    ub.add(1);
                }else ONES.push_back(i);
            }

            Zeros.setBounds(lb,lb); // Fix variables to zero

            CheckIfFix(ONES); // Check if we can fix some location variables to one
            for( iter=ONES.begin(); iter!=ONES.end(); ++iter) y[*iter].setLB(1.0);

            SparseCplex.setParam(IloCplex::CutUp, ObjVal -1 );

            auto SparseStart = CPUclock::now ( );

            if ( solveSparseUsingDualAscentAlg )
            {
                stats->numOfSparseSolved += 1.0;
                double NewObj;
                solveSparseUsingDualAscent ( NewObj , incumbent );
                auto SparseEnd = CPUclock::now ( );
                sparseTime = duration_cast<duration<double>>( SparseEnd - SparseStart ).count ( );
            }
            else
            {
                stats->numOfSparseSolved += 1.0;
                if ( SparseCplex.solve ( ) ){
                    auto SparseEnd = CPUclock::now ( );
                    sparseTime = duration_cast<duration<double>>( SparseEnd - SparseStart ).count ( );
                    if(SparseCplex.getObjValue() < ObjVal){
                        stats->itWhereOptWasFound = it;
                        ObjVal = SparseCplex.getObjValue();
                        incumbent.setObjVal(ObjVal);
                        for(int i=0; i<n; ++i)
                        {
                            ySol[i] = int( SparseCplex.getValue ( y[i] ) );
                        }
                        for ( int j=0; j<m; ++j )
                        {
                            for ( int i=0; i<n; ++i )
                            {
                                if  ( SparseCplex.getValue( x[i][j] >= 0.5 ) )
                                {
                                    xSol[j] = i;
                                    break;
                                }
                            }
                        }
                    }
                }
                else
                {
                    auto SparseEnd = CPUclock::now ( );
                    sparseTime = duration_cast<duration<double>>( SparseEnd - SparseStart ).count ( );
                }
            }


            // Reset bounds
            Zeros.setBounds(lb,ub);
            for(iter = ONES.begin(); iter!=ONES.end(); ++iter) y[*iter].setBounds(0.0,1.0);
            DenseModel.add(piercingcut>=1);
            piercingcut.end();
            Zeros.clear();
            lb.clear();
            ub.clear();

            auto TimeSoFar = CPUclock::now ( ) - RuntimeStart;
            endTime = duration_cast<duration<double>>( TimeSoFar ).count ( );

            std::cout << std::setw(10) << std::left << it
                      << std::setw(10) << std::left << DenseObjVal
                      << std::setw(10) << std::left << ObjVal
                      << std::setw(10) << std::left << std::min(100.0, calcPercent(ObjVal,DenseObjVal))
                      << std::setw(10) << std::left << roundToTwo(denseTime)
                      << std::setw(10) << std::left << roundToTwo(sparseTime)
                      << std::setw(10) << std::left << roundToTwo(endTime) << std::endl;
        }
        auto CnSEnd = CPUclock::now ( );
        auto RuntimeEnd = CPUclock::now ( );
        std::cout << "============================================================================\n";
        std::cout << "The optimal solution value is " << incumbent.getObjVal() << std::endl;

        CnSTime = duration_cast<duration<double>>( CnSEnd - CnSStart ).count ( );
        Runtime = duration_cast<duration<double>>( RuntimeEnd - RuntimeStart ).count ( );

        std::cout << "Solution time           : " << Runtime << " seconds" << std::endl;
        std::cout << "Cutting plane alg. time : " << CutTime << "\n";
        std::cout << "Cut and solve time      : " << CnSTime << "\n";

        stats->bestLowerBound = DenseObjVal;
        stats->bestUpperBound = ObjVal;
        stats->percentageGap = double (ObjVal - DenseObjVal) / double ( DenseObjVal ) * 100;
        stats->time = Runtime;
        stats->CutAndSolveTime = CnSTime;
        stats->CuttingTime = CutTime;
        stats->avgNumOfDualItPerSparse = stats->avgNumOfDualItPerSparse / double ( stats->numOfSparseSolved );
        stats->avgPercentLeft /=stats->numOfSparseSolved;

        return true;
    }catch(std::exception &e){
        std::cerr << "Exception in the Run function in the SSCFLPsolver class: " << e.what() << std::endl;
        return false;
    }catch(IloException &ie){
        std::cerr << "IloException in the Run function in the SSCFLPsolver class: " << ie.getMessage() << std::endl;
        return false;
    }catch(...){std::cerr << "Exception originated in the Run function of the SSCLPsolver class!\n"; return false;}
}

/*********************************************************************************************************/
void SSCFLPsolver::getSolution(int* getY, int* getX){
    try{
        if ( nullptr == xSol )
        {
            std::cout << "Setting the null pointer\n";
            getY = getX = nullptr;
            throw std::runtime_error("No solution available\n");
        }
        else if ( ySol[0] < 0 )
        {
            std::cout << "ySol[0] < 0\n";
            getY = getX = nullptr;
            throw std::runtime_error("No solution available\n");
        }
        else
        {
            std::cout << "Copying the solution\n";
            if ( ( getY == nullptr ) || ( getX == nullptr ) ) std::cout << "Parsing null pointer\n";
            for(int i=0; i<n; ++i) getY[i] = ySol[i];
            for(int j=0; j<m; ++j) getX[j] = xSol[j];
        }
    }catch(std::runtime_error &e){
        std::cerr << "ERROR :::: in getSolution in the SSCFLPsolver class: " << e.what() << std::endl;
    }
}

/*********************************************************************************************************/
void SSCFLPsolver::setC( int i, int j, int cost )
{
    try
    {
        if ( ( i >= n ) || ( i < 0 ) || ( j >= m ) || ( j < 0 ) )
        {
            throw std::runtime_error ("Index out of bounds. No changes made to the C-matrix!\n");
        }
        else
        {
            SparseObj.setLinearCoef( x[i][j] ,  cost );
            DenseObj.setLinearCoef( x[i][j] ,  cost );
            c[i][j] = cost;
        }
    }
    catch ( IloException &ie )
    {
        std::cerr << "IloException in setC in the SSCFLPsolver class : " << ie.getMessage ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( std::exception &e )
    {
        std::cerr << "Exception in setC in the SSCFLPsolver class : " << e.what ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( ... )
    {
        std::cerr << "Unknown error in setC in the SSCFLPsolver class.\n";
        exit ( EXIT_FAILURE );
    }
}

/*********************************************************************************************************/
void SSCFLPsolver::RunLBHeur ( )
{
    try
    {
        auto start = CPUclock::now ( );
        // Record the lower bound
        std::cout << "========= RUNNING LOCAL BRANCHING HEURISTIC =========\n";
        std::cout << "ObjValue before starting LBheur " << ObjVal << std::endl;
        /*
         * Refine the locational decision
         */
        std::cout << "Run refineLocation\n";
        RefineLocation (  );
        std::cout << "After running refineLocation\n";
        /*
         * Run the allocation refiner
         */
        RefineAllocation ( );
        twoSwap ( );

        auto End = CPUclock::now ( );
        auto TotalTime =    std::chrono::duration_cast<std::chrono::duration<double>>( End - start ).count ( );
        std::cout << "\nSolution info: \n";
        std::cout << "Lower bound    : " << CutLowerBound << std::endl;
        std::cout << "Upper bound    : " << ObjVal << std::endl;
        std::cout << "Optimality gap : " << calcPercent ( ObjVal , CutLowerBound ) << std::endl;
        std::cout << "Total time     : " << TotalTime << std::endl;
        std::cout << "=====================================================\n";

        SparseCplex.setParam( IloCplex::MIPEmphasis , 0 ); // Set the focus back to default*/
        SparseCplex.setParam(IloCplex::Param::TimeLimit, 1e75 );
        // Also set bounds on y-variables back
        for ( int i = 0; i<n; ++i )
        {
            y[i].setBounds ( 0 , 1 );
        }

    }
    catch ( IloException &ie )
    {
        std::cerr << "IloException in RunLBHeur in the SSCFLPsolver class : " << ie.getMessage ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( std::exception &e )
    {
        std::cerr << "Exception in RunLBHeur in the SSCFLPsolver class : " << e.what ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( ... )
    {
        std::cerr << "Unknown error in RunLBHeur in the SSCFLPsolver class.\n";
        exit ( EXIT_FAILURE );
    }
}

/*********************************************************************************************************/
void SSCFLPsolver::RefineLocation(  )
{
    try
    {
        bool solutionFound = true,
             greedyHasRun  = false;
        int numOfFacilitiesAtZero   = 0,
            NumOfFlips              = int ( 0.1 * double ( n )  );
        double xVal;
        IloExpr lbLHS ( env ), slbLHS ( env );
        IloRangeArray localBranches ( env );

        DenseCplex.solve ( );

        localBranches.add ( lbLHS <= ( double( numOfFacilitiesAtZero ) / 10.0 ) );

        //SparseModel.add (  localBranches [localBranches.getSize () - 1 ] );
        SparseCplex.setParam( IloCplex::Param::TimeLimit, 15 );
        SparseCplex.setParam ( IloCplex::MIPEmphasis , 1 ); // Focus on finding a feasible solution!


        // Solve for 25 seconds and check if we have a solution!
        SparseCplex.solve ();

        if ( ( SparseCplex.getStatus() == IloAlgorithm::Feasible ) == false )
        {
            // If not, make the local branching constraint redundant, and set full focus on finding a feasible solution
            localBranches[ localBranches.getSize ( ) - 1 ].setBounds ( 0 , IloInfinity );
            SparseCplex.setParam( IloCplex::MIPEmphasis , 4 );

            if ( !SparseCplex.solve() )
            {
                // Now we really should have a feasible solution, if we don't then we simply need to resort to a creazy greedy algorithm!
                greedy () ;
                greedyHasRun = true;
                if ( ObjVal >= INT_MAX - 10 )
                {
                    solutionFound = false;
                }
                else
                {
                    IloNumVarArray start ( env );
                    IloNumArray startVal ( env );
                    for ( int j=0; j<m; ++j )
                    {
                        start.add ( x[xSol[j]][j] );
                        startVal.add ( 1.0 );
                    }
                    for ( int i=0; i<n; ++i )
                    {
                        start.add ( y[i] );
                        startVal.add ( ySol[i] );
                    }
                    SparseCplex.addMIPStart( start , startVal , IloCplex::MIPStartEffort::MIPStartSolveFixed );
                    start.end ( );
                    startVal.end ( );
                }
            }
        }

        // If we have found no solution yet, just give up and go home!
        if ( ! solutionFound ) return;
        // Record the best known solution
        // If the greedy has run, we already know the solution! Else we exstract from cplex
        if ( !greedyHasRun )
        {
            ObjVal = SparseCplex.getObjValue( );
            for ( int j=0; j<m; ++j )
            {
                for ( int i=0; i<n; ++i )
                {
                    if ( j == 0 ) ySol[i] = SparseCplex.getValue ( y[i] );
                    xVal = SparseCplex.getValue ( x[i][j] );
                    if ( xVal >= 0.5 )
                    {
                        xSol[j] = i;
                    }
                }
            }
        }

        if ( !CplexOutOff ) std::cout << "Initial incumbent    : " << ObjVal << std::endl;
        if ( !greedyHasRun ) twoSwap ( );
        // Reset time limit!
        SparseCplex.setParam ( IloCplex::Param::TimeLimit , 5 );

        // Enter the main loop of the local branching algorithm
        for ( int iteration = 0; iteration < 5; ++ iteration )
        {
            // Create the local branching constraint
            for ( int i=0; i<n; ++i )
            {
                if ( ySol[i] == 0 ) lbLHS += y[i];
                else lbLHS += ( 1-y[i] );
            }

            if ( greedyHasRun )
            {
                // The starting solutions is presumably very bad, and we need more room!
                localBranches.add ( lbLHS <= 2*NumOfFlips );
                // Set greedy to false so we do no enter here again
                greedyHasRun = false;
            }
            else
            {
                localBranches.add ( lbLHS <= NumOfFlips );
            }

            SparseModel.add ( localBranches [ localBranches.getSize () - 1 ] );
            if ( iteration == 0 )  for ( int i = 0; i < n; ++i ) for ( int j = 0; j < m; ++j ) x[i][j].setBounds ( 0 , 1 );
            SparseCplex.setParam ( IloCplex::CutUp , ObjVal - 1 );

            if ( !SparseCplex.solve () ) break;

            // Recover the solution for later use
            ObjVal = SparseCplex.getObjValue( );
            for ( int i=0; i<n; ++i ) ySol[i] = SparseCplex.getValue( y[i] );
            for ( int j=0; j<m; ++j )
            {
                for ( int i=0; i<n; ++i )
                {
                    if ( SparseCplex.getValue( x[i][j] ) >= 0.5 )
                    {
                        xSol[j] = i;
                    }
                }
            }

            // Change the local branching constraint to a Right constraint
            localBranches[ localBranches.getSize ( ) - 1 ].setBounds ( NumOfFlips +1 , IloInfinity );

        }

        if ( !CplexOutOff ) std::cout << "Best after loc. ref. : " << ObjVal << std::endl;

        // Fix the closed facilities in the best known solution to zero
        for ( int i=0; i<n; ++i )
        {
            if ( ySol[i] == 0 ) y[i].setUB( 0.0 );
        }
        // Remove all local branching constraints
        for ( int i=0; i<localBranches.getSize (); ++i)
        {
            localBranches[ i ].removeFromAll ( );
        }
        twoSwap ( );

        lbLHS.end ( );
        slbLHS.end ( );
        localBranches.endElements ( );

    }
    catch ( IloException &ie )
    {
        std::cerr << "IloException in RefineLocation in the SSCFLPsolver class : " << ie.getMessage ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( std::exception &e )
    {
        std::cerr << "Exception in RefineLocation in the SSCFLPsolver class : " << e.what ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( ... )
    {
        std::cerr << "Unknown error in RefineLocation in the SSCFLPsolver class.\n";
        exit ( EXIT_FAILURE );
    }

}

/*********************************************************************************************************/
void SSCFLPsolver::RefineAllocation(  )
{
    try
    {
        int iteration     = 0,
            FlipsAllowed  =  int ( 0.1 * double ( m ) );


        IloRangeArray branchArray ( env );
        IloExpr lbLHS ( env );



        for ( iteration = 0; iteration < 5; ++ iteration )
        {
            // Make the local branching constraint
            for ( int j=0; j<m; ++j )
            {
                lbLHS += x[ xSol[j] ][ j ];
            }
            branchArray.add ( lbLHS >= m - FlipsAllowed );
            // Add the local branching constraint!
            SparseModel.add ( branchArray [ branchArray.getSize ( )-1 ] );
            SparseCplex.setParam ( IloCplex::CutUp , ObjVal - 1.0 );

            if (! SparseCplex.solve () ) break;

            // Record the solution
            ObjVal = SparseCplex.getObjValue ( );
            for ( int i=0; i<n; ++i ) ySol[i] = SparseCplex.getValue ( y[i] );
            for ( int j=0; j<m; ++j )
            {
                for ( int i=0; i<n; ++i )
                {
                    if ( SparseCplex.getValue ( x[i][j] ) >= 0.5 )
                    {
                        xSol[j] = i;
                        break;
                    }
                }
            }
            // Make the "right" (as oppose to left) local branching constraint around the allocation variables
            branchArray[ branchArray.getSize ( ) - 1 ].setBounds ( 0 , m - FlipsAllowed - 1);
        }
        if ( !CplexOutOff ) std::cout << "Best after all. ref. : " << ObjVal << std::endl;

        // Remove all the branches we have added!
        for ( int l=0; l<branchArray.getSize ( ); ++l )
        {
            branchArray[ l ].removeFromAll ( );
        }

        branchArray.endElements();
        lbLHS.end();

    }
    catch ( IloException &ie )
    {
        std::cerr << "IloException in RefineAllocation in the SSCFLPsolver class : " << ie.getMessage ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( std::exception &e )
    {
        std::cerr << "Exception in RefineAllocation in the SSCFLPsolver class : " << e.what ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( ... )
    {
        std::cerr << "Unknown error in RefineAllocation in the SSCFLPsolver class.\n";
        exit ( EXIT_FAILURE );
    }
}

/*********************************************************************************************************/
void SSCFLPsolver::greedy()
{
    try
    {
        // Locally defined struct used to stor a customer. Only for ease of readability!
        struct customer{
            int demand;
            int index;
            bool operator<(const customer& c) const
            {
                return demand < c.demand;
            }
        };

        bool hasFoundAssignment = false;
        int  bestC = 0,
             aC = 0,
             TotalCap = 0,
             nextBestC = 0;
        customer cust;
        std::vector<customer> sortC;
        std::vector<customer> unassignedCustomers;
        std::vector<int> resCap;

        item   *items = ( item* ) calloc( n, sizeof( item ) );
        // Start by sorting the customers in non--decreasing order of regret cost!
        for ( int j=0; j<m; ++j )
        {
            bestC = nextBestC = INT_MAX;
            for ( int i = 0; i < n; ++i )
            {
                aC = c[i][j];
                if ( aC < bestC )
                {
                    nextBestC = bestC;
                    bestC = aC;
                }
                else if ( aC < nextBestC )
                {
                    nextBestC = aC;
                }
            }
            cust.demand = nextBestC - bestC;
            cust.index = j;
            sortC.push_back( cust );
        }
        for ( int j=0; j<m; ++j ) sortC[j].demand = d[j];
        std::sort ( sortC.begin(), sortC.end() );
        // Now fill in the residual capacities
        for ( int i=0; i<n; ++i )
        {
            resCap.push_back( s[i] );
            ySol[i] = 0; // For now, assume all closed

            TotalCap += s[i];

            items[i].num = i;
            items[i].p = f[i];
            items[i].w = s[i];
            items[i].x = 0;
        }
        combo ( items, items + n - 1, (TotalCap - TD), 0 , INT_MAX, 1 , 0 );

        for ( int i=0; i<n; ++i )
        {
            ySol[ items[i].num ] = 1 - items[i].x;
        }
        delete[] items;
        // Loop as long as we have unassigned customers!

        while ( !sortC.empty ( ) )
        {
            // Get the last customer ( the one with the largest demand! )
            cust = sortC.back ( );
            // Delete that guy from the back
            sortC.pop_back();
            // Now, find the best insertion of last guy by searching all facilities
            bestC = INT_MAX;
            hasFoundAssignment = false;

            // We start by checking those facilities that are already open
            for ( int i = 0; i < n; ++i )
            {
                if ( ( ySol[i]==1 ) && ( resCap[i] > cust.demand ) )
                {
                    // Facility i is a potential server of cutomer cust
                    hasFoundAssignment = true;
                    // Now check if this is better than best assignemnt
                    if ( c[i][cust.index] < bestC )
                    {
                        xSol[cust.index] = i;
                        bestC = c[i][cust.index];
                    }
                }
            }
            if ( !hasFoundAssignment )
            {
                // If we have not found a solution yet, we need to consider closed facilities
                bestC = INT_MAX;
                for ( int i = 0; i < n; ++i )
                {
                    if ( ( resCap[i] > cust.demand ) )
                    {
                        // Facility i is a potential server of cutomer cust
                        hasFoundAssignment = true;
                        // Now check if this is better than best assignemnt
                        if ( c[i][cust.index] + f[i] < bestC )
                        {
                            xSol[cust.index] = i;
                            bestC = c[i][cust.index] + f[i];
                        }
                    }
                }
            }


            if ( !hasFoundAssignment ) unassignedCustomers.push_back( cust );
            else
            {
                resCap[ xSol[cust.index] ] -= d[cust.index];
                ySol[ xSol[cust.index] ] = 1;
            }
        }

        if ( unassignedCustomers.empty ( ) )
        {
            // We need to set the solution
            ObjVal = 0;
            for ( int j=0; j<m; ++j )
            {
                ObjVal += c[ xSol[j] ][ j ];
                ySol[ xSol[j] ] = 1; // Open the facility customer j is assigned to!
            }
            // Add the fixed costs of the solution
            for ( int i=0; i<n; ++i )
            {
                ObjVal += f[i] * ySol[i];
            }
            std::cout << "Greedy found solution of value : " << ObjVal << std::endl;
            twoSwap ( );
        }
        else
        {
            std::cout << "Number of unassigned customers : " << unassignedCustomers.size ( ) << std::endl;
            ObjVal = INT_MAX;
        }

    }
    catch ( IloException &ie )
    {
        std::cerr << "IloException in greedy in the SSCFLPsolver class : " << ie.getMessage ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( std::exception &e )
    {
        std::cerr << "Exception in greedy in the SSCFLPsolver class : " << e.what ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( ... )
    {
        std::cerr << "Unknown error in greedy in the SSCFLPsolver class.\n";
        exit ( EXIT_FAILURE );
    }
}

/*********************************************************************************************************/
bool SSCFLPsolver::solveSparseUsingDualAscent ( double & NewObjective, solution &incumbent )
{
    try
    {
        bool ProblemSolvedToOptimality = true,
             Continue = true,
             SolutionFeasibel = true;
        int DualLowerBound = 0,
            cstSum = 0;
        double mu = 0,
               tmpMu = 0,
               tmpAllMu = 0,
               maxMu = 0;
        std::vector< int > UnAssigned;
        IloNumVarArray vars ( env );
        IloNumArray vals ( env );

        for ( int i=0; i<n; ++i)
        {
            vars.add ( y[i] );
            cstSum += s[i];
            if ( cstSum < TD ) vals.add ( 1 );
            else vals.add ( 0 );
            maxMu += f[ i ];
            for ( int j=0; j<m; ++j )
            {
                vars.add( x[i][j] );
                vals.add ( 0 );
            }
        }

        SparseCplex.addMIPStart(vars, vals, IloCplex::MIPStartRepair);

        // We start by relaxing the constraints sum(i) x_(i,j)=1 to sum(i) x_(i,j)<=1
        for ( int j = 0; j < m; ++j )
        {
            AssCst[j].setBounds ( 0.0 , 1.0 );

        }


        /* Now we should find an initial value for the lagrangean multiplier
         * If no solution is found yet, we use the duals of the assignments constrins,
         * otherwise we use the largest assignments cost used in the current best solution
         */

        // Initialize the semi lagrangean dual variable

        if ( initializeWithDual1 )
        {
            std::cout << "****************** Initializing with dual information ******************\n";
            mu = 0;
            for ( int j = 0; j < m; ++j )
            {
                if ( duals[ j ] > mu ) mu = std::ceil ( duals[ j ] );
            }

        }
        else if (  initializeWithPercent1 )
        {
            std::cout << "****************** Initializing with percent ******************\n";
            mu = 0;
            int NumOfOpen = 0,
                NumOfAllowed = 0;
            std::vector<int> tmpVec;
            auto nthIt = tmpVec.begin ( );
            for ( int j = 0; j < m; ++j )
            {
                tmpVec.clear ( );
                for ( int i = 0; i < n; ++i )
                {
                    if ( y[i].getUB( ) >= 0.5 )
                    {
                        if ( 0 == j ) ++NumOfOpen;
                        tmpVec.push_back( c[ i ][ j ] );
                    }
                }
                NumOfAllowed = int (  double ( NumOfOpen ) * 0.75 );
                std::nth_element(tmpVec.begin ( ), tmpVec.begin ( ) + NumOfAllowed , tmpVec.end ( ) );
                int cValue =  *( tmpVec.begin () + NumOfAllowed );

                if ( cValue > mu ) mu = cValue;
            }

        }
        else if ( initializeWithSolution1 )
        {
            std::cout << "****************** Initializing with feasible solution ******************\n";
            mu = 0;
            // If we do not have a solution, make one!
            if ( ObjVal >= INT_MAX - 10 ) greedy();
            // If the greedy did not produce a solution resort to the duals!
            if ( ObjVal >= INT_MAX - 10 )
            {
                for ( int j = 0; j < m; ++j )
                {
                    if ( duals[ j ] > mu ) mu = std::ceil ( duals[ j ] );
                }
            }
            else
            {
                for ( int j = 0; j < m; ++j )
                {
                    if ( c[ xSol[ j ] ][ j ] > mu ) mu =  c[ xSol[ j ] ][ j ];
                }
            }

        }




        /*
         * Start the main loop of the dual ascent algorithm
         */
        int iteration = 0;
        double deviation = 0.002;
        while ( Continue )
        {
            stats->avgNumOfDualItPerSparse += 1;
            std::cout << "Iteration " << ++iteration << " of dual ascent. Multiplier value = " << mu << "\n";
            Continue = false; // Assume we solve to optimality this time
            SolutionFeasibel = true; //Assume we produce a feasible solution
            tmpMu = tmpAllMu = INT_MAX;
            UnAssigned.clear ( );

            // Note that we want the dual ascent to produce a solution less than the current best obj: m*mu + L(mu) < ObjVal -1 <=> L(mu) < ObjVal - m*mu - 1
            SparseCplex.setParam ( IloCplex::CutUp , ObjVal - m*mu -1 );

            // Start out by setting the objective coefficients using the multiplier mu
            for ( int i = 0; i < n; ++i )
            {
                if ( y[i].getUB ( ) >= 0.5 )
                {
                    for ( int j = 0; j < m; ++j )
                    {
                        SparseObj.setLinearCoef( x[ i ][ j ] , c[ i ][ j ] - mu );
                        if ( c[ i ][ j ] > mu ) x[ i ][ j ].setBounds ( 0 , 0 );
                        else x[ i ][ j ].setBounds( 0 , 1 );
                    }
                }
            }

            // Tell cplex to abort when solution within 1% is found
            if ( mu  <= maxC )
            {
                deviation /= 1.1;
                SparseCplex.setParam( IloCplex::EpGap , deviation );
            }
            else
            {
                SparseCplex.setParam( IloCplex::EpGap , 0.0 );
            }

            // Solve the problem heuristically
            if ( SparseCplex.solve ( ) )
            {
                // Check the slack of each of the assignment constraints
                SolutionFeasibel = true;
                for ( int j = 0; j < m; ++j )
                {
                    cstSum = 0;
                    for ( int i = 0; i < n; ++i ) cstSum += SparseCplex.getValue ( x[i][j] );

                    if ( cstSum == 0 )
                    {
                        UnAssigned.push_back( j );
                        SolutionFeasibel = false;
                        break;
                    }
                }

                std::cout << "Solution is feasible : " << SolutionFeasibel << std::endl;

                // Now test if solution is feasible as we then need to solve to optimality
                if ( SolutionFeasibel )
                {
                    // Set the optimality gap to zero
                    SparseCplex.setParam ( IloCplex::EpGap , 0.0 );
                    // Solve to optimality. Note we do not need to test if SparseCplex contains a solution as we already have a solution within 1% of optimality!
                    SparseCplex.solve ( );

                    // Update the dual lower bound
                    DualLowerBound = m*mu + SparseCplex.getObjValue ( );
                    std::cout << "Dual lower bound : " << DualLowerBound << std::endl;
                    if ( DualLowerBound >= ObjVal )
                    {
                        for ( int j = 0; j < m; ++j )
                        {
                            AssCst[j].setBounds ( 1 , 1 );
                            for ( int i = 0; i < n; ++i )
                            {
                                SparseObj.setLinearCoef( x[ i ][ j ] , c[ i ][ j ] );
                                 x[ i ][ j ].setBounds ( 0 , 1 );
                            }
                        }

                        return false;
                    }
                    else
                    {
                        UnAssigned.clear ( );
                        SolutionFeasibel = false;
                        for ( int j = 0; j < m; ++j )
                        {
                            cstSum = 0;
                            for ( int i = 0; i < n; ++i ) cstSum += SparseCplex.getValue ( x[i][j] );

                            if ( cstSum == 0 )
                            {
                                UnAssigned.push_back( j );
                                SolutionFeasibel = false;
                                break;
                            }
                        }
                    }
                }



                // If the flag Continue = false, then the solution is optimal and we can break the while loop
                if ( UnAssigned.empty ( ) )
                {
                    std::cout << "It seems no customers are unassigned\n";
                    break;
                }
                else if ( !UnAssigned.empty ( ) ) // Otherwise we must update the dual multiplier
                {
                    Continue = true;
                    for ( auto it = UnAssigned.begin ( ); it != UnAssigned.end ( ); ++it )
                    {
                        for ( int i = 0; i < n; ++i )
                        {
                            if ( ( c[ i ][ *it ] > mu ) && ( c[ i ][ *it ] < tmpMu ) ) tmpMu = c[ i ][ *it ];
                        }
                    }
                    // Test if we improved the multiplier
                    if ( tmpMu < INT_MAX -100 )
                    {
                        mu = std::max ( mu + 10.0 , tmpMu );
                    }
                    else
                    { // We need to look through all possible assignments!
                        tmpMu = INT_MAX;
                        maxMu = 0;
                        for ( int i = 0; i < n; ++i )
                        {
                            maxMu += f[ i ];
                            for ( int j = 0; j < m; ++j )
                            {
                                if ( ( c[ i ][ j ] > mu ) && ( c[ i ][ j ] < tmpMu ) )
                                {
                                    tmpMu = c[ i ][ j ];
                                }
                            }
                        }

                        // Check if we improved using this rule, else simply max out!
                        if ( tmpMu < INT_MAX - 100 ) mu = std::max ( mu + 10.0 , tmpMu );
                        else
                        {
                            for ( int j = 0; j < m; ++j ) AssCst[ j ].setBounds ( 0.0 , 1.0 );
                            mu = maxMu + maxC + 1;

                        }

                    }

                }


            }
            else
            {
                for ( int j = 0; j < m; ++j )
                {
                    AssCst[ j ].setBounds ( 1 , 1 );
                    for ( int i = 0; i < n; ++i )
                    {
                        SparseObj.setLinearCoef( x[ i ][ j ] , c[ i ][ j ] );
                         x[ i ][ j ].setBounds ( 0 , 1 );
                    }
                }
                return false;
            }
        }

        // Retrieve the solution
        std::vector<int> xx ( m );
        std::vector<int> yy ( n );
        int NumOfAssLeft = 0;
        for ( int j = 0; j < m; ++j )
        {
            if ( 0 == j ) for ( int i = 0; i < n; ++i ) yy[i] = SparseCplex.getValue ( y[i] ); // Only check this the first time
            for ( int i = 0; i < n; ++i )
            {
                if ( ( y[i].getUB() >= 0.5 ) && ( c[ i ][ j ] < mu ) ) ++NumOfAssLeft;
                if ( SparseCplex.getValue( x[i][j] ) >= 0.5 )
                {
                    xx[j] = i;
                    break;
                }
            }
        }

        stats->avgPercentLeft += double ( NumOfAssLeft ) / double ( m * n );

        DualLowerBound = 0;
        for ( int i = 0; i < n; ++i ) DualLowerBound += f[i] * yy[i];
        for ( int j = 0; j < m; ++j ) DualLowerBound += c[ xx[ j ] ][ j ];
        if ( DualLowerBound < ObjVal )
        {
            stats->itWhereOptWasFound = stats->numberOfIterations;
            ObjVal = DualLowerBound;
            for ( int i = 0; i < n; ++i ) ySol[i] = yy[i];
            for ( int j = 0; j < m; ++j ) xSol[j] = xx[j];
            incumbent.setObjVal(ObjVal);
        }

        // Set the constraints back to equality constraints
        for ( int j = 0; j < m; ++j )
        {
            AssCst[j].setBounds ( 1 , 1 );
            for ( int i = 0; i < n; ++i )
            {
                SparseObj.setLinearCoef( x[i][j] , c[i][j] );
                x[ i ][ j ].setBounds ( 0 , 1 );
            }
        }

        return ProblemSolvedToOptimality;
    }
    catch ( IloException &ie )
    {
        std::cerr << "IloException in solveSparseUsingDualAscent in the SSCFLPsolver class : " << ie.getMessage ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( std::exception &e )
    {
        std::cerr << "Exception in solveSparseUsingDualAscent in the SSCFLPsolver class : " << e.what ( ) << "\n";
        exit ( EXIT_FAILURE );
    }
    catch ( ... )
    {
        std::cerr << "Unknown error in solveSparseUsingDualAscent in the SSCFLPsolver class.\n";
        exit ( EXIT_FAILURE );
    }
}


/*********************************************************************************************************/
void SSCFLPsolver::twoSwap ( )
{
    try
    {
        bool foundImprovement = true;
        int CurCost = 0,
            newCost = 0,
            Imp     = 0,
            BestImp = 0,
            tempI    = 0;
        std::vector<int> resCap ( n );

        // Calculate the residual capacity of the facilities before we start the two swap heuristic
        for ( int i = 0; i < n; ++i )
        {
            resCap[i] = s[i];
        }
        for ( int j=0; j<m; ++j )
        {
            resCap[ xSol[j] ] -= d[j];
        }

        /*
         * Start of main loop
         */
        while ( foundImprovement )
        {
            foundImprovement = false;
            BestImp = 0;
            for ( int j = 0; j < m - 1; ++j )
            {
                for ( int k = j; k < m; ++k )
                {
                    // Test if there is enough capacity to make the swap
                    if ( ( resCap[ xSol[j] ] + d[j] >= d[k] ) && ( resCap[ xSol[k] ] + d[k] >= d[j] ) )
                    {
                        CurCost = c[ xSol[j] ][ j ] + c[ xSol[k] ][ k ];
                        newCost = c[ xSol[j] ][ k ] + c[ xSol[k] ][ j ];
                        // Test if that makes an improvement
                        Imp = CurCost - newCost;
                        if ( Imp >  0 )
                        {
                            resCap[ xSol[j] ] += ( d[j] - d[k] );
                            resCap[ xSol[k] ] += ( d[k] - d[j] );
                            tempI = xSol[j];
                            xSol[j] = xSol[k];
                            xSol[k] = tempI;
                            ObjVal -= Imp;
                            foundImprovement = true;
                        }
                    }
                }
            }
        }
        std::cout << "ObjVal after two swap : " << ObjVal << std::endl;
    }
    catch( std::exception &e )
    {
        std::cerr << "Error in twoSwap in the SSCFLPsolver class : " << e.what ( ) << std::endl;
    }
}

/*********************************************************************************************************/
bool SSCFLPsolver::RunAsHeuristic ( )
{
    // Initialize arrays to hold solution
    ySol = new int[n];
    ySol[0] = -1; // Used to indicate no solution has been found
    xSol = new int[m];

    // Tell Cplex to shut up
    DenseCplex.setOut( env.getNullStream ( ) );
    DenseCplex.setWarning( env.getNullStream ( ) );
    SparseCplex.setOut ( env.getNullStream ( ) );
    SparseCplex.setWarning( env.getNullStream ( ) );
    SparseCplex.setParam (   IloCplex::MIPSearch, 1);
    SparseCplex.setParam ( IloCplex::ParallelMode , 1 );

    /*
     * Fix variables before adding cutting planes
     */
    IloRangeArray Fixations ( env );
    for ( auto it = fixBeforeCuttingPhase.begin ( ); it != fixBeforeCuttingPhase.end ( ); ++it )
    {
        Fixations.add ( x[ it->first ][ it->second ] == 0 );
        DenseModel.add ( Fixations[ Fixations.getSize ( ) - 1 ] );
        SparseModel.add ( Fixations[ Fixations.getSize ( ) - 1 ] );
    }

    //Run the cutting phase to strengthen then lower bound
    CuttingPhase();

    /*
     * Starting to fix the variables that should be fixed after cutting phase
     */
    for ( auto it = fixAfterCuttingPhase.begin ( ); it != fixAfterCuttingPhase.end ( ); ++it )
    {
        Fixations.add ( x[ it->first ][ it->second ] == 0 );
        DenseModel.add ( Fixations[ Fixations.getSize ( ) - 1 ] );
        SparseModel.add ( Fixations[ Fixations.getSize ( ) - 1 ] );
    }


    // Run the local branching heuristic!
    RunLBHeur ( );

    Fixations.removeFromAll ( );
    hasRun = true;
    return ( ObjVal <=  INT_MAX -10 );
}

