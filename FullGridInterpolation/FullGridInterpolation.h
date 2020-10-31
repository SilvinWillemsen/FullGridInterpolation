//
//  FullGridInterpolation.hpp
//  FullGridInterpolation
//
//  Created by Silvin Willemsen on 31/10/2020.
//  Copyright © 2020 Silvin Willemsen. All rights reserved.
//

#ifndef FullGridInterpolation_h
#define FullGridInterpolation_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "math.h"

class FullGridInterpolation
{
public:
    FullGridInterpolation (int startN, int endN,
                           double fs, double outLength,
                           double excitationLoc, double excitationWidth,
                           double outputLocStart);
    
    ~FullGridInterpolation();
    
    void calculate();

    void scheme();
    void updateStates();
    
    void retrieveState();
private:
    int startN, endN, lengthSound;
    double fs, outLength, excitationLoc, excitationWidth, outputLocStart;
    double startC, endC, cDiff;
    
    int N, NPrev, outLoc;
    double k, h, c, lambdaSq;
    double B0, B1, C;
    std::vector<std::vector<double>> uVecs;
    std::vector<double*> u;
    
    std::vector<double> out;
    
    std::ofstream stateAt;
};

#endif /* FullGridInterpolation_h */
