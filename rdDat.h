#ifndef RDDAT_H_INCLUDED
#define RDDAT_H_INCLUDED

#include<random>
#include<exception>
#include<stdexcept>
#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<string>


/* An introduction to the rdDat class
* \author Sune Lauth Gadegaard
* \version 1.0.0
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
*
\section{Description}
*
* Class implementing a data reader for different kinds of discrete facility location problems.
*
* \section Compiling
* The codes were compiled using the GNU GCC compiler on a Linux Ubuntu 14.04 machine.
* The following flags were used: -Wall -O3 -std=c++11 -DIL_STD.
* The Code::blocks IDE was used as well.
*
* \latexonly
* \section{Change log for rdDat.h and rdDat.cpp}
* \begin{center}
*     \begin{tabularx}{\textwidth}{llr X}\toprule
*        FILE:          &   \multicolumn{3}{l}{rdDat.h and rdDat.cpp}\\
*        Version:       &   \multicolumn{3}{l}{1.0.0}\\
*        \multicolumn{4}{l}{- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -}\\
*        CHANGE LOG:    &   DATE    &   VER.-NO.    &   CHANGES MADE\\ \midrule
*                       &   2016--01--28    &   1.0.0       & First implementation\\ \bottomrule
*     \end{tabularx}
* \end{center}
* \endlatexonly
*/

class rdDat{
    private:
        int n;      //!< Number of facilities
        int m;      //!< Number of customers
        int p;      //!< Number of open faiclities in a solution to the p-median problem
        int TotalDemand; //!< Sum of the demands
        std::vector<int> d;     //!< Demands. d[j] demand of customer
        std::vector<int> s;     //!< Capacities. s[i] capacity of facility i
        std::vector<double> f;     //!< Fixed opening cost. f[i] cost of opening facility i
        std::vector<std::vector<double>> c; //!< Assingment cost. c[i][j] is the cost of supplying all of customer j's demand from facility i


        int TheProblemType; //!< Integer indicating which problem type is in qustion
    public:

        /*!
         * Constructor of the rdDat class.
         * \param Filename  String. Contains the path to a data file of appropriate format
         */
        rdDat(std::string Filename, int ProblemType);

        /*!
         * Destructor of the class. Cleans up after the us.
         */
        ~rdDat();

        /*!
         * Reads the data of a p-median problem
         */
        void rdPmed(std::string DataFile);

        /*!
         * Reads the data of an uncapacitated facility location problem problem
         */
        void rdUFLP(std::string DataFile);

        /*!
         * Reads the data of a capacitated facility location problem problem
         */
        void rdCFLP(std::string DataFile);

        /*!
         * Reads the data of a single source capacitated facility location problem problem
         */
        void rdSSCFLP(std::string DataFile);

        /*!
         * Returns the number of facilities
         */
        inline
        int getNumFac(){return n;}

        /*!
         * Returns the number of customers
         */
        inline
        int getNumCust(){return m;}

        /*!
         * Returns the number facilities which must be open in a p-median problem
         */
        inline
        int getP(){return p;}

        /*!
         * Returns the demand of customer j
         * \param j integer. Index of the customer who's demand you want
         */
        inline
        int getD(int j){return d[j];}

        /*!
         * Returns the capacity of facility i
         * \param i integer. Index of the facility who's capacity you want
         */
        inline
        int getS(int i){return s[i];}

        /*!
         * Returns the fixed opening cost of facility i
         * \param ji integer. Index of the facility who's fixed cost you want
         */
        inline
        int getF(int i){return f[i];}

        /*!
         * Returns the assignment cost of the the facility--customer pair (i,j)
         * \param i integer. Index of the facility
         * \param j integer. Index of the customer
         */
        inline
        int getC(int i, int j){return c[i][j];}



        /*!
         * Returns a pointer to the first element in the integer array containing the demands
         */
        inline
        std::vector<int> getAllD(){return d;}

        /*!
         * Returns a pointer to the first element in the integer array containing the capacities
         */
        inline
        std::vector<int> getAllS(){return s;}

        /*!
         * Returns a pointer to the first element in the integer array containing the fixed opening costs
         */
        inline
        std::vector<double> getAllF(){return f;}

        /*!
         * Returns a pointer to a pointer to the an integer array containing the assignment costs
         */
        inline
        std::vector<std::vector<double>> getAllC(){return c;}



        /**!
         * Adds an extra facility to the data
         * \param cap constant re
         */
        void addFacility ( int cap, int fixed, std::vector<double> newC );

        /**!
         * Set the assignment cost c[i][j]
         */
        void setC ( int i, int j, int cost ) { c[i][j] = cost; }

        /*!
         * Returns the total demand
         */
        int getTotalDemand ( ) { return TotalDemand; }
};

#endif // RDDAT_H_INCLUDED
