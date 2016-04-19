#include "rdDat.h"

rdDat::rdDat(std::string Filename, int ProblemType){
    try{
        TheProblemType = ProblemType;

        rdSSCFLP(Filename);


    }catch(std::exception & e){
        std::cerr << "Exception in CFLP constructor: " << e.what() << std::endl;
        exit(1);
    }
}
//======================================================================================================
void rdDat::rdPmed(std::string DataFile){
    try{
        double anyNum;
        std::ifstream InFile;
        std::stringstream err;
        if(!InFile) throw std::runtime_error("Cannot create the ifstream!\n");

        InFile.open(DataFile);
        if(!InFile){
            err << "Cannot open the input file:\t" << DataFile << "\n";
            throw std::runtime_error(err.str());
        }

        /* Format:
         * n = number of nodes in the graph. We set n=m just to be sure.
         * p = The number of open facilites in an optimal solution.
         * c_ij = cost matrix
         * t_ij = travel time matrix
         */
        InFile>>anyNum;
        if(!InFile) throw std::runtime_error("Could not retrieve the number of facilities.\n");
        if( anyNum < 1 ) throw std::runtime_error("A non-positive number of vertices has been specified. That is not meaningfull!\n");
        n = anyNum;
        InFile>>anyNum;
        if(!InFile) throw std::runtime_error("Could not retrieve the number of facilities.\n");
        if( anyNum < 1 ) throw std::runtime_error("A non-positive number of vertices has been specified. That is not meaningfull!\n");
        m = anyNum;

        InFile >> anyNum;
        if(!InFile) throw std::runtime_error("Could not retrieve the number of open facilities, p.\n");
        p = anyNum;

        c = std::vector<std::vector<double>>( n );

        for(int i=0; i<n; ++i){ //  Reading the cost matrix
            c[i] = std::vector<double>(m);
            for(int j=0; j<n; ++j){
                InFile >> anyNum;
                if( !InFile ){
                    err << "Could not read the cost matrix entry " << i << "," << j << std::endl;
                    throw std::runtime_error(err.str());
                }
                c[i][j] = anyNum;
            }
        }



    }catch(std::exception &e){
        std::cerr << "Exception in rdPmed in the rdDat class: " << e.what() << std::endl;
    }
}

//======================================================================================================
void rdDat::rdSSCFLP(std::string DataFile){
    try{
        double anyInt;
        std::ifstream InFile;
        std::stringstream err;
        TotalDemand = 0;
        if(!InFile) throw std::runtime_error("Cannot create the ifstream!\n");

        InFile.open(DataFile);

        if(!InFile){
            err << "Cannot open the input file:\t" << DataFile << "\n";
            throw std::runtime_error(err.str());
        }

        if ( TheProblemType != 1 )
        {
            // Format: Holmberg, yang, gadegaard
            // n (number of facilities)
            // m (number of customers)
            // s_1 f_1
            // s_2 f_2
            // ...
            // s_n f_n
            //d_1 d_2 d_3 ... d_m
            // c_ij
            InFile>>anyInt;
            if(!InFile) throw std::runtime_error("Could not retrieve the number of facilities.\n");
            n = anyInt;
            InFile >> anyInt;
            if(!InFile) throw std::runtime_error("Could not retrieve the number of customers.\n");
            m = anyInt;

            c = std::vector<std::vector<double>>( n );
            d = std::vector<int>( m );
            s = std::vector<int>( n );
            f = std::vector<double> ( n );

            for(int i=0; i<n; ++i){
                InFile >> anyInt;
                if(!InFile){
                    err << "Cannot read capacity " << i << ". Terminating!\n";
                    throw std::runtime_error(err.str());
                }
                s[i] = anyInt;
                InFile >> anyInt;
                if(!InFile){
                    err << "Cannot read fixed cost " << i << ". Terminating!\n";
                    throw std::runtime_error(err.str());
                }
                f[i] = anyInt;
            }

            for(int j=0; j<m; ++j){
                InFile >> anyInt;
                if(!InFile){
                    err << "Cannot read demand " << j << ". Terminating!\n";
                    throw std::runtime_error(err.str());
                }
                d[j] = anyInt;
                TotalDemand += d[j];
            }
            for(int i=0; i<n; ++i){
                c[i] = std::vector<double>( m );
                for(int j=0; j<m; ++j){
                    InFile >> anyInt;
                    if(!InFile){
                        err << "Cannot read assignement cost " << i << " " << j << ". Terminating!\n";
                        throw std::runtime_error(err.str());
                    }
                    c[i][j]=anyInt;
                }
            }
        }
        else
        {
            TotalDemand = 0;
            double anyNum;
            InFile >> m;
            if(!InFile) throw std::runtime_error("Could not read the number of customers. Termintaing!\n");

            InFile >> n;
            if(!InFile) throw std::runtime_error("Could not read the number of facilities. Termintaing!\n");

            if(n<=0 || m<=0){
                err << "Non positve values for number of facilities and customers. Provided was n=" << n << ", and m=" << m << ". Terminating\n";
                throw std::runtime_error(err.str().c_str());
            }
            for(int i=0; i<n; ++i) c.push_back(std::vector<double>( m ) );
            for(int j=0; j<m; ++j){ // Reading the assingment costs
                for(int i=0; i<n; ++i){
                    InFile >> anyNum;
                    if(!InFile || anyNum < 0){
                        err << "Could not read the assignment cost of index (" << i << "," << j << "j). Terminating!\n";
                        throw std::runtime_error(err.str().c_str());
                    }
                    c[i][j]=int(anyNum);
                }
            }

            for(int j=0; j<m; ++j){ // Reading the demands
                InFile >> anyNum;
                if(!InFile || anyNum<=0){
                    err << "Could not read the demand of customer " << j << ". Terminating!\n";
                    throw std::runtime_error(err.str().c_str());
                }
                d.push_back( int(anyNum) );
                TotalDemand+=d[j];
            }

            for(int i=0; i<n; ++i){ // Reading the fixed opening costs
                InFile >> anyNum;
                if(!InFile || anyNum<0 ){
                    err << "Could not read the fixed opening cost of facility " << i << ". Terminating!\n";
                    throw std::runtime_error(err.str().c_str());
                }
                f.push_back( int(anyNum) );
            }

            for(int i=0; i<n; ++i){
                InFile >> anyNum;
                if(!InFile || anyNum<=0 ){
                    err << "Could not read the capacity of facility " << i << ". Terminating!\n";
                    throw std::runtime_error(err.str().c_str());
                }
                s.push_back( int(anyNum) );
            }
        }
    }catch(std::exception &e){
        std::cerr << "Exception in rdSSCFLP in rdDat class: " << e.what() << std::endl;
    }
}

//======================================================================================================
void rdDat::addFacility( int cap, int fixed, std::vector<double> newC )
{
    f.push_back( fixed );
    s.push_back( cap );
    c.push_back( newC );
    ++n;
}

//======================================================================================================
rdDat::~rdDat(){

}
