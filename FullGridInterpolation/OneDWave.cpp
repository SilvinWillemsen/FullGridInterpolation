//
//  OneDWave.cpp
//  OneDWave
//
//  Created by Silvin Willemsen on 31/10/2020.
//  Copyright © 2020 Silvin Willemsen. All rights reserved.
//

#include "OneDWave.h"

OneDWave::OneDWave (int startN, int endN,
                  double fs, double outLength,
                  double excitationLoc, double excitationWidth,
                  double outputLocStart) :   startN (startN), endN (endN),
    fs (fs), outLength (outLength),
    excitationLoc (excitationLoc),
    excitationWidth(excitationWidth),
    outputLocStart (outputLocStart)
{
    lengthSound = outLength * fs;
    
    k = 1.0 / fs;
    
    startC = 44100 / static_cast<double> (startN);
    c = startC;
    endC = 44100 / static_cast<double> (endN);
    
    cDiff = startC - endC;
    h = c*k;
    N = floor (1.0 / h);
    h = 1 / static_cast<double> (N);
    NPrev = N;
    lambdaSq = c * c * k * k / (h * h);
    
    uVecs.resize (3, std::vector<double> (N-1, 0));
    u.resize(3, nullptr);
    
    for (int i = 0; i < uVecs.size(); ++i)
        u[i] = &uVecs[i][0];
    
    
    int loc = excitationLoc * N;
    int width = excitationWidth * N;
    int raisedCosStart = floor (loc - width  / 2) - 1;
    int raisedCosEnd = floor (loc + width / 2) - 1;
    
    for (int i = raisedCosStart; i <= raisedCosEnd; ++i)
    {
        u[1][i] = 0.5 * (1 - cos (2.0 * M_PI * (i-raisedCosStart) / width));
        u[2][i] = u[1][i];
    }
    
    out.resize (lengthSound, 0);
    outLoc = outputLocStart * N;
    
    B0 = 2.0 - 2.0 * lambdaSq;
    B1 = lambdaSq;
    C = -1.0;
    
    // file writing
    stateAt.open ("stateAt.csv");
}

OneDWave::~OneDWave()
{
    stateAt.close();
}

//void OneDWave::calculate()
//{
//    
//}
void OneDWave::scheme()
{
    
    //// Inner scheme ////
    for (int l = 1; l < N-2; ++l)
    {
        u[0][l] = B0 * u[1][l] + B1 * (u[1][l+1] + u[1][l-1]) + C * u[2][l];
    }
    
    //// Boundaries ////
    u[0][0] = B0 * u[1][0] + B1 * u[1][1] + C * u[2][0];
    u[0][N-2] = B0 * u[1][N-2] + B1 * u[1][N-3] + C * u[2][N-2];

}


void OneDWave::updateStates()
{
    double* uTmp = u[2];
    u[2] = u[1];
    u[1] = u[0];
    u[0] = uTmp;
}

void OneDWave::retrieveState()
{
    for (int l = 0; l < N-1; ++l)
        stateAt << std::to_string(u[1][l]) << ";\n";
}

void OneDWave::fullGridInterpolation()
{
    
}
