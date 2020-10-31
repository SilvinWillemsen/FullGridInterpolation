//
//  OneDWave.hpp
//  OneDWave
//
//  Created by Silvin Willemsen on 31/10/2020.
//  Copyright © 2020 Silvin Willemsen. All rights reserved.
//

#ifndef OneDWave_h
#define OneDWave_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "math.h"

class OneDWave
{
public:
    OneDWave (int startN, int endN,
               double fs, double outLength,
               double excitationLoc, double excitationWidth,
               double outputLocStart);
    
    ~OneDWave();
    
    void calculate();

    void scheme();
    void updateStates();
    
    void retrieveState();
    void fullGridInterpolation();
    
private:
    int sampleAtWhichToRetrieveState;
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
    int curPercentage = 0;
    int test = 0; 
};

#endif /* OneDWave_h */
