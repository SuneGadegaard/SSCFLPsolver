#ifndef SOLUTION_H_INCLUDED
#define SOLUTION_H_INCLUDED

/**
* \author Sune Lauth Gadegaard
* \file solution.h
* \date 2015--04-09
* \version 1.0.0
*
*  Class for storing a solution to the single source capacitated facility location problem. Implemented in C++.
*/


#include<exception>
#include<stdexcept>
#include<iostream>
#include<algorithm>
class solution{
    private:
        int ObjVal;     //!< Objective function value of the solutuion.
        int* solY;      //!< Pointer to array of integers. solY[i] = 1 if facility i is open in the solution. Otherwise solY[i]=0.
        int* solX;      //!< Pointer to array of integers. solX[j] = i if and only customer j is assingned to facility i in the solution.
        int n;          //!< Integer. The number of facilities in the instance of the SSCFLP for which this is a solution.
        int m;          //!< Integer. The number of customers in the instance of the SSCFLP for which this is a solution.
    public:
        /*! \brief Overloaded constructor of the solution class.
         *
         * Constructor of the solution class. Initializes the number of facilities and the number of customers.
         * Allocattes memory for the solY and the solX arrays. Initializes arrays to zero, that is solY[i]=solX[j]=0 for all i and j
         * \param nn Integer. The number of facilities.
         * \param mm Integer. The number of customers.
         */
        solution(int nn, int mm);

        /*! \brief Overloaded constructor of the solution class
         *
         * Constructor of the solution class. Initializes the number of facilities and the number of customers.
         * Allocattes memory for the solY and the solX arrays. Sets the value of solY[i]=newY[i] and solX[j]=newX[j].
         * If newY==NULL or newX==NUL ObjVal is set to -1 indicating the solution is rubbish.
         * \param nn Integer. The number of facilities.
         * \param mm Integer. The number of customers.
         * \param newY Pointer to array of integers. Contains the y-values for the solution you want to store. newY[i]=1 iff facility i is open.Ŕ
         * \param newX Pointer to array of integers. Contains the assignment of the solution you want to store. newX[j]=i iff customer j is assingned to facility i.
         * \param obj Integer. The objective function value of the soltuion.
         */
        solution(int nn, int mm, int* newY, int* newX, int obj);

        /*! \brief Destructor of the solution class.
         *
         * Destructor of the solution class. Clears all allocated memory.
         */
        ~solution();

        /*! \brief Set the solution
         *
         * Sets the solution equal to the new solution provided by (newY, newX).
         * If either newY = NULL or newX = NULL ObjVal is set to -1 to indicate no valid solution is present.
         * \param newY pointer to array of integers. Must be of size at least n. newY[i]=1 iff facility i is open in the solution you want to store.
         * \param newX pointer to array of integers. Must be of size at least m. newX[j]=i iff customer j is assigned to faiclity i in the solution you want to store.
         * \param obj Integer. The objective function value of the solution provided.
         */
        void setSolution(int* newY, int* newX, int obj);

        /*! \brief Returns the solution to the y-variables.
         *
         * Returns the solution to the y-variables stored in the object. If setSolution() or solution(int,int,*int,*int) has not been called, getSolution(int*,int*) will return rubbish.
         * \param   getY pointer to integer array of size at least n. Contains on output the solution to the y variables stored in the solution object.
         *          getY[i]=1 iff faiclity i is open in the solution.
         * \param   getX pointer to integer array of size at least m. Contains on output the solution to the assignments stored in the solution object.
         *          getX[j]=i iff customer j is assinged to facility i in the solution.
         */
        void getSolution(int* getY, int* getX);

        /*!
         * \brief Sets the objective function value of the solution
         */
         void setObjVal(int &obj){ObjVal=obj;}

        /*! \brief Gives the solution value of solY[i].
         *
         * Gives the solution value of solY[i]. If i>=n or i<0 the last or first element will be returned, respectively. No error is thrown, but a message is displayed.
         */
        int getY(int i);

        /*! \brief Gives the solution value of solX[j].
         *
         * Gives the solution value of solX[j]. If j>=m or j<0 the last or first element will be returned, respectively. No error is thrown, but a message is displayed.
         */
        int getX(int j);

        /*! \brief Returns the objective function value of the solution.
         *
         * Returns the objective function value of the solution. Remember if ObjVal=-1 it means the solution is rubbish!
         */
        int getObjVal(){return ObjVal;}
};

#endif // SOLUTION_H_INCLUDED
