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
#include "../eigen/Eigen/Eigen"

enum DynamicInterpolationType
{
    dLinear,
    dCubic,
    dSinc
};

struct SincInterpolVals
{
    SincInterpolVals (int sincWidth = 0, double alphaBand = 0) :
    sincWidth (sincWidth), alphaBand (alphaBand)
    {}
    int sincWidth;
    double alphaBand;
};

class OneDWaveDynamic
{
public:
    OneDWaveDynamic (double startN, double endN,
              double fs, double outLength,
              double excitationLoc, double excitationWidth,
              double outputLocStart, DynamicInterpolationType dyIntType,
              SincInterpolVals& sIV,
              double lambdaMultiplier,
              bool changeC, int numFromRightBound,
              bool lpConnection, double lpExponent,
              bool LFO, double LFOfreq = 0,
              double changeS = 0, double changeE = 1);
    /*
        The numFromRightBound variable determines Number from the right boundary (quite important, switches between different techniques)
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
    void addRemoveAtPoint();

    void createCustomIp();
    void lpPoints(); // low-pass the connection

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
    int numFromRightBound;
    
    void (OneDWaveDynamic::*addRemovePoints)();

    bool LFO;
    double LFOFreq;
    double changeStart;
    double changeEnd;
    double changeDiff;

    // LowPass the connection
    bool lpConnection;
    double lpExponent;
    
    
    SincInterpolVals sIV;
    int sincWidth, iLen;
    double bmax;
    std::vector<double> xUMp1;
    Eigen::VectorXd aU, bU, aW;
    Eigen::MatrixXd AU;
};

#endif /* OneDWaveDynamic_h */

