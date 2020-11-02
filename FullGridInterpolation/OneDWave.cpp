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
    h = c * k;
    N = floor (1.0 / h);
    h = 1 / static_cast<double> (N);
    NPrev = N;
    lambdaSq = c * c * k * k / (h * h);

    emptyVector.resize (std::max(startN, endN)-1, -1);
    uVecs1 = std::make_shared<std::vector<std::vector<double>>> (3, emptyVector);
    uVecs2 = std::make_shared<std::vector<std::vector<double>>> (3, emptyVector);

//    uVecs1.resize (3, emptyVector);
//    uVecs2.resize (3, emptyVector);

    u.resize(3, nullptr);
//    auto test = uVecs1.get()[0][0][0];
//    auto test = uVecs1.get()[0][0];
    for (int i = 0; i < uVecs1->size(); ++i)
        u[i] = &uVecs1.get()[0][i][0];
    
    pointingAtuVecs1 = true;
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
    plotIdx.open ("plotIdx.csv");
    output.open ("output.csv");

}

OneDWave::~OneDWave()
{
    stateAt.close();
    plotIdx.close();
    output.close();

}

void OneDWave::scheme()
{
    
    //// Inner scheme ////
    for (int l = 1; l < N-2; ++l)
    {
        u[0][l] = B0 * u[1][l] + B1 * (u[1][l+1] + u[1][l-1]) + C * u[2][l];
    }
    
    //// Next to boundaries ////
    u[0][0] = B0 * u[1][0] + B1 * u[1][1] + C * u[2][0];
    u[0][N-2] = B0 * u[1][N-2] + B1 * u[1][N-3] + C * u[2][N-2];

    if (retrievingState)
        retrieveState (N);
}


void OneDWave::updateStates()
{
    double* uTmp = u[2];
    u[2] = u[1];
    u[1] = u[0];
    u[0] = uTmp;
}

void OneDWave::retrieveOutput (int indexFromBoundary, bool fromRightBoundary)
{
    output << u[1][fromRightBoundary ? N-1-indexFromBoundary : indexFromBoundary - 1] << ";\n";
}

void OneDWave::retrieveState (int end)
{
    plotIdx << curPlotIdx << ";\n";
    for (int l = 0; l < end-1; ++l)
        stateAt << std::to_string(u[1][l]) << ";\n";
    curPlotIdx += end-2;
    plotIdx << curPlotIdx << ";\n";
    ++curPlotIdx;
}

void OneDWave::fullGridInterpolation (bool isNGrowing)
{
    if (isNGrowing)
    {
        double alphSpace = NPrev / static_cast<double> (N);
        double alph = alphSpace;
        double m = alphSpace - 1;

        if (pointingAtuVecs1)
        {
            uVecs2.get()[0].resize (3, emptyVector);

            for (int l = 0; l < N; ++l)
            {
//                std::cout << m << ", " << alph << std::endl;
                uVecs2.get()[0][1][l] = cubicInterpolation (u[1], floor(m), alph);
                uVecs2.get()[0][2][l] = cubicInterpolation (u[2], floor(m), alph);
                m += alphSpace;
                alph += alphSpace;
                if (alph >= 1)
                    alph -= 1;
            }
//            uVecs2[1][0] = 0;
//            uVecs2[1][N-1] = 0;
//            uVecs2[2][0] = 0;
//            uVecs2[2][N-1] = 0;
            
            for (int i = 0; i < uVecs2->size(); ++i)
                u[i] = &uVecs2.get()[0][i][0];

        } else {
            uVecs1.get()[0].resize (3, emptyVector);

            for (int l = 0; l < N; ++l)
            {
//                std::cout << m << ", " << alph << std::endl;
                uVecs1.get()[0][1][l] = cubicInterpolation (u[1], floor(m), alph);
                uVecs1.get()[0][2][l] = cubicInterpolation (u[2], floor(m), alph);
                m += alphSpace;
                alph += alphSpace;
                if (alph >= 1)
                    alph -= 1;
            }
            for (int i = 0; i < uVecs1->size(); ++i)
                u[i] = &uVecs1.get()[0][i][0];
        }
        
        
    } else { // Add the division h1/h2 division here somewhere from stefans book
        
    }
    pointingAtuVecs1 = !pointingAtuVecs1;

}

void OneDWave::recalculateCoeffs (int n)
{
    c = startC - n / static_cast<double> (lengthSound-1) * cDiff;
    h = c * k;
    N = floor (1.0 / h);
    h = 1 / static_cast<double> (N);
    lambdaSq = c * c * k * k / (h * h);
    if (N != NPrev)
    {
        retrieveState (NPrev);
        fullGridInterpolation (N > NPrev);
        retrieveState (N);
//        stateAt.close();

    }
    NPrev = N;

    
    B0 = 2.0 - 2.0 * lambdaSq;
    B1 = lambdaSq;
    C = -1.0;
}

double OneDWave::cubicInterpolation (double* uVec, int l, double alph)
{
    // if l == 0 we're at u_1 as we're not including the boundaries in the calculations
    if (l == -1)
    {
        return -uVec[l + 1] * (alph * (alph - 1) * (alph - 2)) / -6.0 // simply supported
//        + uVec[l] * ((alph - 1) * (alph + 1) * (alph - 2)) * 0.5 // is 0 anyway
        + uVec[l + 1] * (alph * (alph + 1) * (alph - 2)) * -0.5
        + uVec[l + 2] * (alph * (alph + 1) * (alph - 1)) / 6.0;
    } else if (l == 0)
    {
        return // uVec[l - 1] * (alph * (alph - 1) * (alph - 2)) / -6.0 // is 0 anyway
        uVec[l] * ((alph - 1) * (alph + 1) * (alph - 2)) * 0.5
        + uVec[l + 1] * (alph * (alph + 1) * (alph - 2)) * -0.5
        + uVec[l + 2] * (alph * (alph + 1) * (alph - 1)) / 6.0;
    } else if (l == N-2) {
        return uVec[l - 1] * (alph * (alph - 1) * (alph - 2)) / -6.0
        + uVec[l] * ((alph - 1) * (alph + 1) * (alph - 2)) * 0.5
        + uVec[l + 1] * (alph * (alph + 1) * (alph - 2)) * -0.5;
//        + uVec[l + 2] * (alph * (alph + 1) * (alph - 1)) / 6.0; // is 0 anyway
    }  else if (l == N-1) {
            return uVec[l - 1] * (alph * (alph - 1) * (alph - 2)) / -6.0
            + uVec[l] * ((alph - 1) * (alph + 1) * (alph - 2)) * 0.5
//            + uVec[l + 1] * (alph * (alph + 1) * (alph - 2)) * -0.5;// is 0 anyway
             - uVec[l] * (alph * (alph + 1) * (alph - 1)) / 6.0; // simply supported boundary
    } else {
        return uVec[l - 1] * (alph * (alph - 1) * (alph - 2)) / -6.0
            + uVec[l] * ((alph - 1) * (alph + 1) * (alph - 2)) * 0.5
            + uVec[l + 1] * (alph * (alph + 1) * (alph - 2)) * -0.5
            + uVec[l + 2] * (alph * (alph + 1) * (alph - 1)) / 6.0;
    }
    
//    double val = 0;
//    if (bp > simplySupported ? 2 : 3)
//        val = val + uVec[bp - 1] * (alph * (alph - 1) * (alph - 2)) * -oneOverSix;
//
//    val = val + uVec[bp] * ((alph - 1) * (alph + 1) * (alph - 2)) * 0.5;
//
//    if (bp < N - (simplySupported ? 2 : 3))
//        val = val + uVec[bp + 1] * (alph * (alph + 1) * (alph - 2)) * -0.5;
//    if (bp < N - (simplySupported ? 3 : 4))
//        val = val + uVec[bp + 2] * (alph * (alph + 1) * (alph - 1)) * oneOverSix;
    //    std::cout << val - val1 << std::endl;
//    return val;
}
