#include"solution.h"

solution::solution(int nn, int mm){
    try{
        n = nn;
        m = mm;
        solY = new int[n]();
        solX = new int[m]();
        ObjVal=-1;
    }catch(std::exception &e){
        std::cerr << "Exception in constructor (no 1) of the solution class: " << e.what() << std::endl;
        return;
    }catch(...){std::cerr << "Exception originated in constructor (no 1) of the solution clas!\n"; return;}
}

/***********************************************************************************************************/
solution::solution(int nn, int mm, int* newY, int* newX, int obj){
    try{
        n = nn;
        m = mm;
        solY = new int[n];
        solX = new int[m];
        if(!newY && !newX){
            for(int k=0; k<std::max(n,m); ++k){
                if(k<n) solY[k] = newY[k];
                if(k<m) solX[k] = newX[k];
            }
            ObjVal = obj;
        }else{
            ObjVal = -1;
        }


    }catch(std::exception &e){
        std::cerr << "Exception in the constructor (no 2) of the solution class: " << e.what() << std::endl;
        return;
    }catch(...){std::cerr << "Exception originated in the constructor (no 2) of the solution class!\n"; return;}
}

/***********************************************************************************************************/
void solution::setSolution(int* newY, int* newX, int obj){
    try{
        if(newY && newX){ for(int k=0; k<std::max(n,m); ++k){ if(k<n) solY[k] = newY[k]; if(k<m) solX[k] = newX[k];} ObjVal = obj;}
        else ObjVal =-1;
    }catch(std::exception &e){
        std::cerr << "Exception in setSolution in the solution class: " << e.what() << std::endl;
        return;
    }
}


/***********************************************************************************************************/
void solution::getSolution(int* getY, int* getX){
    try{
        for(int k=0; k<std::max(n,m); ++k){ if(k<n) getY[k] = solY[k]; if(k<m) getX[k] = solX[k]; }
    }catch(std::exception &e){
        std::cerr << "Exception in getSolution in the solution class: " << e.what() << std::endl;
        return;
    }
}

/***********************************************************************************************************/
int solution::getY(int i){
    if(i<0 || i>=n) std::cout << "\n\nAsking for index " << i << " in getY. This is out of bounds\n";
    if(i<0)return solY[0];
    else if(i>=n) return solY[n-1];

    return solY[i];
}

/***********************************************************************************************************/
int solution::getX(int j){
    if(j<0 || j>=m) std::cout << "\n\nAsking for index " << j << " in getX. This is out of bounds\n";
    if(j<0)return solX[0];
    else if(j>=m)return solX[m-1];

    return solX[j];
}

/***********************************************************************************************************/
solution::~solution(){
    delete[] solY;
    delete[] solX;
}
