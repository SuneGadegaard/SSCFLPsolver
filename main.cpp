#include"SSCFLPsolver.h"
#include <chrono>
int main(int argc, char** argv){
    try{
        SSCFLPsolver solver = SSCFLPsolver();

        solver.Load(argv[1],2);
        solver.printStats();
        if(solver.Run()){
            std::cout << "Hurra!\n";
        }
        int n = solver.getNumFac();
        int m = solver.getNumCust();

        int* y = new int[n];
        int* x = new int[m];

        solver.getSolution(y,x);

        delete[] y;
        delete[] x;

        return 0;
    }catch(std::exception &e){
        std::cerr << "Exception: " << e.what() << std::endl;
    }catch(...){ std::cerr << "An unexpected exception was thrown. Caught in main.\n"; return 1;}
}
