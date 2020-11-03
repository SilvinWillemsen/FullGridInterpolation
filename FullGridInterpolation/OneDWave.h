//
//  OneDWave.hpp
//  OneDWave
//
//  Created by Silvin Willemsen on 31/10/2020.
//  Copyright © 2020 Silvin Willemsen. All rights reserved.
//

/*
 Implementation of the 1D wave equation with changing parameters and full grid interpolation. The most important thing to note is that the boundary points are not included in the implementation. So u[1][0] is u_1^n and u[1][N-2] is u_{N-1}^n. Or, in other words there are N-1 moving points ranging from idx = 0 : N-2.
 This also means that in the cubic interpolation we start at idx = -1 so that we can reach a point before the first described point.
 */

#ifndef OneDWave_h
#define OneDWave_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "math.h"

enum InterpolationType
{
    none,
    linear,
    cubic
};

class OneDWave
{
public:
    OneDWave (int startN, int endN,
               double fs, double outLength,
               double excitationLoc, double excitationWidth,
               double outputLocStart, InterpolationType interpolationType);
    
    ~OneDWave();
    
    void recalculateCoeffs (int n);
    void scheme();
    void updateStates();
    
    void retrieveState (int end);
    void retrieveOutput (int indexFromBoundary, bool fromRightBoundary);

    void fullGridInterpolation (bool isNGrowing);
    double cubicInterpolation (double* uVec, int l, double alph);
    double linearInterpolation (double* uVec, int l, double alph);


private:

    int startN, endN, lengthSound;
    double fs, outLength, excitationLoc, excitationWidth, outputLocStart;
    double startC, endC, cDiff;
    
    int N, NPrev, outLoc;
    double k, h, c, lambdaSq;
    double B0, B1, C;
    std::vector<double> emptyVector;
    std::vector<std::vector<double>> uVecs1;
    std::vector<std::vector<double>> uVecs2;

    std::vector<double*> u;
    
    std::vector<double> out;
    
    std::ofstream stateAt, plotIdx, output, cSave, NSave, NChange, lambdaSqSave;
    
    bool pointingAtuVecs1;
    int curPlotIdx = 1;
    InterpolationType interpolationType;
};

#endif /* OneDWave_h */
