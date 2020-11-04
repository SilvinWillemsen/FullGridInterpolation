//
//  OneDWaveDynamic.cpp
//  FullGridInterpolation
//
//  Created by Silvin Willemsen on 03/11/2020.
//  Copyright © 2020 Silvin Willemsen. All rights reserved.
//

#include "OneDWaveDynamic.h"

OneDWaveDynamic::OneDWaveDynamic (double startN, double endN,
                    double fs, double outLength,
                    double excitationLoc, double excitationWidth,
                    double outputLocStart, double lambdaMultiplier) :
startN (startN), endN (endN),
fs (fs), outLength (outLength),
excitationLoc (excitationLoc),
excitationWidth(excitationWidth),
outputLocStart (outputLocStart),
lambdaMultiplier (lambdaMultiplier)
{
    lengthSound = outLength * fs;
    
    k = 1.0 / fs;
    
    startC = 44100 / startN;
    c = startC;
    endC = 44100 / endN;
    
    cDiff = startC - endC;
    h = c * k;
    N = floor (1.0 / h);
    h = 1 / static_cast<double> (N);
    NPrev = N;
    lambdaSq = lambdaMultiplier * c * c * k * k / (h * h);
    
    Mu = ceil (N * 0.5);
    Mw = floor (N * 0.5);
    uVecs.resize (3, std::vector<double> (ceil (std::max (double(Mu), fs / 44100.0 * endN) * 0.5), 0));
    wVecs.resize (3, std::vector<double> (floor (std::max (double(Mw), fs / 44100.0 * endN) * 0.5), 0));
    
    u.resize(3, nullptr);
    w.resize(3, nullptr);

    for (int i = 0; i < uVecs.size(); ++i)
    {
        u[i] = &uVecs[i][0];
        w[i] = &wVecs[i][0];
    }
    
    int loc = excitationLoc * N;
    int width = std::max (2.0, excitationWidth * N);
    int raisedCosStart = floor (loc - width * 0.5) - 1;
    int raisedCosEnd = floor (loc + width * 0.5) - 1;
    
    // check if overlap is part of raised cos
    if (raisedCosStart <= Mu && raisedCosEnd >= Mu)
    {
        std::cout << "Overlap!" << std::endl;
        stopSimulation = true;
        return;
    }
    else if (raisedCosEnd < Mu)
    {
        for (int i = raisedCosStart; i <= raisedCosEnd; ++i)
        {
            u[1][i] = 0.5 * (1 - cos (2.0 * M_PI * (i-raisedCosStart) / width));
            u[2][i] = u[1][i];
        }
    }
    else if (raisedCosEnd > Mu)
    {
        raisedCosStart -= Mu;
        raisedCosEnd -= Mu;
        for (int i = raisedCosStart; i <= raisedCosEnd; ++i)
        {
            w[1][i] = 0.5 * (1 - cos (2.0 * M_PI * (i-raisedCosStart) / width));
            w[2][i] = w[1][i];
        }
    }
    
    outLoc = outputLocStart * N;
    
    customIp.resize (4, 0);
    customIp1.resize (4, 0);

    B0 = 2.0 - 2.0 * lambdaSq;
    B1 = lambdaSq;
    C = -1.0;
    
    a11 = 1;
    a22 = 1;
    
    
    // file writing
    stateAt.open ("stateAtD.csv");
    plotIdx.open ("plotIdxD.csv");
    output.open ("outputD.csv");
    cSave.open ("cSaveD.csv");
    NSave.open ("NSaveD.csv");
    NChange.open ("NChangeD.csv");
    lambdaSqSave.open ("lambdaSqSaveD.csv");
    alfTickSave.open ("alfTickSaveD.csv");
    
}

OneDWaveDynamic::~OneDWaveDynamic()
{
    stateAt.close();
    plotIdx.close();
    output.close();
    cSave.close();
    NSave.close();
    NChange.close();
    lambdaSqSave.close();
    alfTickSave.close();
}

void OneDWaveDynamic::recalculateCoeffs (int n)
{
    c = startC - n / static_cast<double> (lengthSound-1) * cDiff;
    h = c * k;
    NDouble = 1.0 / h;
    N = floor (NDouble);
    
    lambdaSq = lambdaMultiplier * c * c * k * k / (h * h);
    lambdaSqSave << lambdaSq << ";\n";
    cSave << c << ";\n";
    NSave << N << ";\n";
    
    alf = NDouble - N;
//    if (n % 441 == 0)

//    if (n > 10)
//    {
//        std::cout << std::endl;
//        stopSimulation = true;
//    }
    
    if (N != NPrev)
    {
        NChange << NPrev << ", " << N << ";\n";
        if (abs(N - NPrev) > 1)
            std::cout << "too fast..?" << std::endl;
        alfTick = ((1.0-Mw * h) - (Mu * h + h)) / h;
        alfTickSave << alfTick << ";\n";
        if (round(round(alfTick * 10) / 15.0) == 1)
            std::cout << alfTick - alf << std::endl;

//        alfTick = alf;
//        std::cout << alfTick - alf << std::endl;
        if (N > NPrev) // add point
        {
            customIp[0] = -alfTick * (alfTick + 1.0) / ((alfTick + 2.0) * (alfTick + 3.0));
            customIp[1] = 2.0 * alfTick / (alfTick + 2.0);
            customIp[2] = 2.0 / (alfTick + 2.0);
            customIp[3] = -2.0 * alfTick / ((alfTick + 3.0) * (alfTick + 2.0));
            
            customIp1[3] = alf * (alf - 1) * (alf - 2) / -6.0;
            customIp1[2] = (alf - 1) * (alf + 1) * (alf - 2) / 2.0;
            customIp1[1] = alf * (alf + 1) * (alf - 2) / -2.0;
            customIp1[0] = alf * (alf + 1) * (alf - 1) / 6.0;
            if (N % 2 == 1)
            {
                u[1][Mu] = customIp[0] * u[1][Mu-2]
                    + customIp[1] * u[1][Mu-1]
                    + customIp[2] * w[1][0]
                    + customIp[3] * w[1][1];
                u[2][Mu] = customIp[0] * u[2][Mu-2]
                    + customIp[1] * u[2][Mu-1]
                    + customIp[2] * w[2][0]
                    + customIp[3] * w[2][1];
                ++Mu;

            }
            else
            {
                // move w vector one up (can be optimised)
                for (int l = Mw-1; l >= 0; --l)
                {
                    w[1][l+1] = w[1][l];
                    w[2][l+1] = w[2][l];

                }
                w[1][0] = customIp[3] * u[1][Mu-2]
                    + customIp[2] * u[1][Mu-1]
                    + customIp[1] * w[1][0]
                    + customIp[0] * w[1][1];
                w[2][0] = customIp[3] * u[2][Mu-2]
                    + customIp[2] * u[2][Mu-1]
                    + customIp[1] * w[2][0]
                    + customIp[0] * w[2][1];
                ++Mw;
            }
        } else {
            if (N % 2 == 0)
            {
                u[1][Mu-1] = 0;
                u[2][Mu-1] = 0;
                --Mu;
                
            }
            else
            {
                // move w vector one down (can be optimised)
                for (int l = 0; l < Mw; ++l)
                {
                    w[1][l] = w[1][l+1];
                    w[2][l] = w[2][l+1];
                }
                w[1][Mw-1] = 0;
                w[2][Mw-1] = 0;

                --Mw;
            }
        }
    }
    NPrev = N;
    
    calculateInterpolatedPoints();

    
    B0 = 2.0 - 2.0 * lambdaSq;
    B1 = lambdaSq;
    C = -1.0;
}


void OneDWaveDynamic::calculateInterpolatedPoints()
{
    alpha = alf * (alf - 1) * (alf - 2) / -6.0;
    beta = (alf - 1) * (alf + 1) * (alf - 2) * 0.5;
    gamma = alf * (alf + 1) * (alf - 2) * -0.5;
    delta = alf * (alf + 1) * (alf - 1) / 6.0;
    
    a21 = -delta;
    a12 = -delta;
    
    oOdet = 1.0 / (a11 * a22 - a12 * a21);
    v1 = alpha * w[1][2] + beta * w[1][1] + gamma * w[1][0];
    v2 = alpha * u[1][Mu-3] + beta * u[1][Mu-2] + gamma * u[1][Mu-1];
    
//    if (v2 != 0)
//    {
//        std::cout << "wait" << std::endl;
//        stopSimulation = true;
//    }
    uMuP1 = (v1 * a22 - v2 * a12) * oOdet;
    wMin1 = (-v1 * a21 + v2 * a11) * oOdet;
}

void OneDWaveDynamic::scheme()
{
    
    //// Inner scheme ////
    for (int l = 1; l < Mu-1; ++l)
    {
        u[0][l] = B0 * u[1][l] + B1 * (u[1][l+1] + u[1][l-1]) + C * u[2][l];
    }
    
    //// Next to boundaries ////
    u[0][0] = B0 * u[1][0] + B1 * u[1][1] + C * u[2][0];
    u[0][Mu-1] = B0 * u[1][Mu-1] + B1 * (u[1][Mu-2] + uMuP1) + C * u[2][Mu-1];
    
    //// Inner scheme ////
    for (int l = 1; l < Mw-1; ++l)
    {
        w[0][l] = B0 * w[1][l] + B1 * (w[1][l+1] + w[1][l-1]) + C * w[2][l];
    }
    
    //// Next to boundaries ////
    w[0][Mw-1] = B0 * w[1][Mw-1] + B1 * w[1][Mw-2] + C * w[2][Mw-1];
    w[0][0] = B0 * w[1][0] + B1 * (w[1][1] + wMin1) + C * w[2][0];
    
}


void OneDWaveDynamic::updateStates()
{
    double* uTmp = u[2];
    u[2] = u[1];
    u[1] = u[0];
    u[0] = uTmp;
    
    double* wTmp = w[2];
    w[2] = w[1];
    w[1] = w[0];
    w[0] = wTmp;
}

void OneDWaveDynamic::retrieveOutput (int indexFromBoundary, bool fromRightBoundary)
{
//    output << w[1][fromRightBoundary ? (Mw-1-indexFromBoundary) : (indexFromBoundary - 1)] << ";\n";
    output << w[1][std::max (0, (Mw-1-indexFromBoundary))] << ";\n";
}

void OneDWaveDynamic::retrieveStateU()
{
    plotIdx << curPlotIdx << ";\n";
    for (int l = 0; l < Mu; ++l)
        stateAt << std::to_string(u[1][l]) << ";\n";
    curPlotIdx += Mu-1;
    plotIdx << curPlotIdx << ";\n";
    ++curPlotIdx;
//    std::cout << "Last index of uNext " << u[0][Mu-1] << std::endl;
    
}

void OneDWaveDynamic::retrieveStateW()
{
    plotIdx << curPlotIdx << ";\n";
    for (int l = 0; l < Mw; ++l)
        stateAt << std::to_string(w[1][l]) << ";\n";
    curPlotIdx += Mw-1;
    plotIdx << curPlotIdx << ";\n";
    ++curPlotIdx;
    
//    std::cout << "First index of wNext " << w[0][0] << std::endl;
}
