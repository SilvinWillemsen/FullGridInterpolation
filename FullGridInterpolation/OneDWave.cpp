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
                  double outputLocStart, InterpolationType interpolationType) :
    startN (startN), endN (endN),
    fs (fs), outLength (outLength),
    excitationLoc (excitationLoc),
    excitationWidth(excitationWidth),
    outputLocStart (outputLocStart),
    interpolationType (interpolationType)
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

    emptyVector.resize (std::max(fs / 44100.0 * startN, fs / 44100.0 * endN)-1, 0);
    uVecs1.resize (3, emptyVector);
    uVecs2.resize (3, emptyVector);
    u.resize(3, nullptr);
    
    for (int i = 0; i < uVecs1.size(); ++i)
        u[i] = &uVecs1[i][0];
    
    pointingAtuVecs1 = true;
    int loc = excitationLoc * N;
    int width = excitationWidth * N;
    int raisedCosStart = floor (loc - width * 0.5) - 1;
    int raisedCosEnd = floor (loc + width * 0.5) - 1;
    
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
    cSave.open ("cSave.csv");
    NSave.open ("NSave.csv");
    NChange.open ("NChange.csv");
    lambdaSqSave.open ("lambdaSqSave.csv");


}

OneDWave::~OneDWave()
{
    stateAt.close();
    plotIdx.close();
    output.close();
    cSave.close();
    NSave.close();
    NChange.close();
    lambdaSqSave.close();
}

void OneDWave::recalculateCoeffs (int n)
{
    c = startC - n / static_cast<double> (lengthSound-1) * cDiff;
    if (interpolationType != none)
    {
        h = c * k;
        N = floor (1.0 / h);
        h = 1 / static_cast<double> (N);
    }
    lambdaSq = 0.999 * c * c * k * k / (h * h);
    lambdaSqSave << lambdaSq << ";\n";
    cSave << c << ";\n";
    NSave << N << ";\n";
    
    //    if (n % int (fs / 44.10) == 0)
    //        retrieveState (N);
    
    if (N != NPrev)
    {
        NChange << NPrev << ", " << N << ";\n";
        if (abs(N - NPrev) > 1)
            std::cout << "too fast..?" << std::endl;
        
        fullGridInterpolation (N > NPrev);
    }
    NPrev = N;
    
    B0 = 2.0 - 2.0 * lambdaSq;
    B1 = lambdaSq;
    C = -1.0;
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

}


void OneDWave::updateStates()
{
    double* uTmp = u[2];
    u[2] = u[1];
    u[1] = u[0];
    u[0] = uTmp;
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
            uVecs2.resize (3, std::vector<double> (std::max(fs / 44100.0 * startN, fs / 44100.0 * endN)-1, 0));
            
            if (interpolationType == linear)
            {
                for (int l = 0; l < N-1; ++l)
                {
                    uVecs2[1][l] = linearInterpolation (u[1], floor(m), alph);
                    uVecs2[2][l] = linearInterpolation (u[2], floor(m), alph);
                    m += alphSpace;
                    alph += alphSpace;
                    if (alph >= 1)
                        alph -= 1;
                }
            }
            else if (interpolationType == cubic)
            {
                for (int l = 0; l < N-1; ++l)
                {
                    uVecs2[1][l] = cubicInterpolation (u[1], floor(m), alph);
                    uVecs2[2][l] = cubicInterpolation (u[2], floor(m), alph);
                    m += alphSpace;
                    alph += alphSpace;
                    if (alph >= 1)
                        alph -= 1;
                }
            }
            
            for (int i = 0; i < uVecs2.size(); ++i)
                u[i] = &uVecs2[i][0];

        } else {
            uVecs1.resize (3, std::vector<double> (std::max(fs / 44100.0 * startN, fs / 44100.0 * endN)-1, 0));
            
            if (interpolationType == linear)
            {
            for (int l = 0; l < N-1; ++l)
            {
                uVecs1[1][l] = linearInterpolation (u[1], floor(m), alph);
                uVecs1[2][l] = linearInterpolation (u[2], floor(m), alph);
                m += alphSpace;
                alph += alphSpace;
                if (alph >= 1)
                    alph -= 1;
            }
            }
            
            if (interpolationType == cubic)
            {
                for (int l = 0; l < N-1; ++l)
                {
                    uVecs1[1][l] = cubicInterpolation (u[1], floor(m), alph);
                    uVecs1[2][l] = cubicInterpolation (u[2], floor(m), alph);
                    m += alphSpace;
                    alph += alphSpace;
                    if (alph >= 1)
                        alph -= 1;
                }
            }
            for (int i = 0; i < uVecs1.size(); ++i)
                u[i] = &uVecs1[i][0];
        }
        
        
    } else { // Add the division h1/h2 division here somewhere from stefans book
        
    }
    pointingAtuVecs1 = !pointingAtuVecs1;

}

double OneDWave::linearInterpolation (double* uVec, int l, double alph)
{
    if (l == -1)
        return alph * uVec[l+1];
    else if (l == NPrev-2)
        return (1.0 - alph) * uVec[l];
    
    return (1.0 - alph) * uVec[l] + alph * uVec[l+1];
    
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
    }
    else if (l == 0)
    {
        return // uVec[l - 1] * (alph * (alph - 1) * (alph - 2)) / -6.0 // is 0 anyway
            uVec[l] * ((alph - 1) * (alph + 1) * (alph - 2)) * 0.5
            + uVec[l + 1] * (alph * (alph + 1) * (alph - 2)) * -0.5
            + uVec[l + 2] * (alph * (alph + 1) * (alph - 1)) / 6.0;
    }
    else if (l == NPrev-3)
    {
        return uVec[l - 1] * (alph * (alph - 1) * (alph - 2)) / -6.0
            + uVec[l] * ((alph - 1) * (alph + 1) * (alph - 2)) * 0.5
            + uVec[l + 1] * (alph * (alph + 1) * (alph - 2)) * -0.5;
//        + uVec[l + 2] * (alph * (alph + 1) * (alph - 1)) / 6.0; // is 0 anyway
    }
    else if (l == NPrev-2)
    {
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


void OneDWave::retrieveOutput (int indexFromBoundary, bool fromRightBoundary)
{
    output << u[1][fromRightBoundary ? (N-1-indexFromBoundary) : (indexFromBoundary - 1)] << ";\n";
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

