//
//  OneDWaveDynamic.hpp
//  FullGridInterpolation
//
//  Created by Silvin Willemsen on 03/11/2020.
//  Copyright © 2020 Silvin Willemsen. All rights reserved.
//

#ifndef OneDWaveDynamic_h
#define OneDWaveDynamic_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "math.h"

enum DynamicInterpolationType
{
    dLinear,
    dCubic
};


class OneDWaveDynamic
{
public:
    OneDWaveDynamic (double startN, double endN,
              double fs, double outLength,
              double excitationLoc, double excitationWidth,
              double outputLocStart, DynamicInterpolationType dyIntType,
              double lambdaMultiplier,
              bool changeC, int whereToAddPoints);
    /*
        The whereToAddPoints variable determines Number from the right boundary (quite important, switches between different techniques)
        -1: Adding to the center alternating between left and right string.
        0: Interpolated boundary
        1: Right string has a single moving point. Using simply supported boundary condition
        2: Right string has two moving points. When trying to solve the cubic
        interpolation, w_2 is always 0 (that's why this is a bit different)
                                        >3: (Expected behaviour) Selects where to add points (to left string).
    */
    ~OneDWaveDynamic();
    
    void recalculateCoeffs (int n);
    void calculateInterpolatedPoints();
    void scheme();
    void updateStates();
    
    void retrieveStateU();
    void retrieveStateW();

    void retrieveOutput (int indexFromBoundary, bool fromRightBoundary);
    
    bool simulationStopped() { return stopSimulation; };
    
    void addRemoveInCenter();
    
private:
    int lengthSound;
    double startN, endN, startNTrue, endNTrue, NDiff;
    double fs, outLength, excitationLoc, excitationWidth, outputLocStart, lambdaMultiplier;
    double startC, endC, cDiff;
    
    int N, NPrev, outLoc, Mu, Mw;
    double k, h, c, lambdaSq, NDouble, alf, alfTick;
    double B0, B1, C;

    // calculating interpolated points at the connection
    double a11, a12, a21, a22;
    double v1, v2;
    double alpha, beta, gamma, delta, oOdet;
    
    std::vector<std::vector<double>> uVecs;
    std::vector<std::vector<double>> wVecs;
    
    std::vector<double*> u;
    std::vector<double*> w;
    double uMuP1 = 0;
    double wMin1 = 0;

    std::vector<double> customIp, customIp1;
    
    std::ofstream stateAt, plotIdx, output, cSave, NSave, NChange, lambdaSqSave, alfTickSave;
    
    int curPlotIdx = 1;
    
    bool stopSimulation = false;
    bool changeC;
    
    DynamicInterpolationType dyIntType;
    int whereToAddPoints;
    
    void (OneDWaveDynamic::*addRemovePoints)();

};

#endif /* OneDWaveDynamic_h */

