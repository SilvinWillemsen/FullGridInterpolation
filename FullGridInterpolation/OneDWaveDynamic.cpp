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
                    double outputLocStart, DynamicInterpolationType dyIntType, bool even, bool shiftedIn,
                    SincInterpolVals& sIV, double lambdaMultiplier,
                    bool changeC, int numFromRightBound,
                    bool lpConnection, double lpExponent, double lpDelay,
                    bool LFO, double LFOFreq,
                    double changeS, double changeE) :
lengthSound (outLength * fs),
startN (startN), endN (endN),
fs (fs), outLength (outLength),
excitationLoc (excitationLoc),
excitationWidth(excitationWidth),
outputLocStart (outputLocStart),
dyIntType (dyIntType),
even (even),
shifted (even ? false : shiftedIn),
lambdaMultiplier (lambdaMultiplier),
changeC (changeC), numFromRightBound (numFromRightBound),
LFO (LFO), LFOFreq (LFOFreq),
changeStart (round(changeS * lengthSound)),
changeEnd (round(changeE * lengthSound)),
lpConnection (lpConnection),
lpExponent (lpExponent),
lpDelay (lpDelay),
sIV (sIV),
sincWidth (sIV.sincWidth)
{
    bool initWithDiffAtConn = false;
    
    lengthSound = outLength * fs;
    
    k = 1.0 / fs;
    
    startC = 44100 / startN;
    c = startC;
    endC = 44100 / endN;
    if (!LFO)
    {
        if (changeC)
        {
            cDiff = startC - endC;
        } else {
            startNTrue = startN * fs / 44100.0;
            endNTrue = endN * fs / 44100.0;
            NDiff = endNTrue - startNTrue;
        }
    }
    
    h = c * k;
    N = floor (1.0 / h);
    h = 1 / static_cast<double> (N);
    NPrev = N;
    lambdaSq = lambdaMultiplier * c * c * k * k / (h * h);
    
    
//    if (numFromRightBound > N*0.5) // only add at right boundary
//    {
//        stopSimulation = true;
//        std::cout << "NumFromRightBound is too large" << std::endl;
//        return;
//    }
    
    // Resize u and w vecs according to where points are added
    switch (numFromRightBound)
    {
        case -1:
        {
            double maxM = std::max (double(N), fs / 44100.0 * endN) * 0.5;
            Mu = ceil (N * 0.5);
            Mw = floor (N * 0.5);
            uVecs.resize (3, std::vector<double> (ceil (maxM), 0));
            wVecs.resize (3, std::vector<double> (floor (maxM), 0));
            addRemovePoints = &OneDWaveDynamic::addRemoveInCenter;
            break;
        }
        default:
        {
            Mu = N - numFromRightBound;
            Mw = numFromRightBound;
            double maxM = std::max (double(N - numFromRightBound), fs / 44100.0 * (endN - numFromRightBound));
            uVecs.resize (3, std::vector<double> (maxM, 0));
            wVecs.resize (3, std::vector<double> (Mw, 0));
            addRemovePoints = &OneDWaveDynamic::addRemoveAtPoint;
            break;

        }
    }
    u.resize(3, nullptr);
    w.resize(3, nullptr);

    for (int i = 0; i < uVecs.size(); ++i)
    {
        u[i] = &uVecs[i][0];
        w[i] = &wVecs[i][0];
    }
    if (initWithDiffAtConn)
    {
        uVecs[1][Mu-1] = 0.5;
        wVecs[1][0] = -0.5;
//        uVecs[2][Mu-1] = 0.5;
//        wVecs[2][0] = -0.5;
    } else {
        double loc = excitationLoc * N;
        int width = std::max (4.0, excitationWidth * N);
        int raisedCosStart = floor (loc - width * 0.5) - 1;
        
        int raisedCosEnd = floor (loc + width * 0.5) - 1 - std::min(raisedCosStart, 0);
        raisedCosStart = raisedCosStart - std::min(raisedCosStart, 0);
        // should probably also do this for the right bound
        
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
        else if (raisedCosStart > Mu)
        {
            raisedCosStart -= Mu;
            raisedCosEnd -= Mu;
            for (int i = raisedCosStart; i <= raisedCosEnd; ++i)
            {
                w[1][i] = 0.5 * (1 - cos (2.0 * M_PI * (i-raisedCosStart) / width));
                w[2][i] = w[1][i];
            }
        }
    }
    
////    initialise with nyquist
//    for (int i = 0; i < uVecs[0].size(); ++i)
//    {
//        u[1][i] = (i % 2) ? 0.5 : -0.5;
//        u[2][i] = (i % 2) ? 0.5 : -0.5;
//    }
//    w[1][0] = u[1][Mu-1];
//    w[2][0] = u[2][Mu-1];
    
    outLoc = outputLocStart * N;
    
    customIp.resize (4, 0);
    customIp1.resize (4, 0);

    B0 = 2.0 - 2.0 * lambdaSq;
    B1 = lambdaSq;
    C = -1.0;
    
    a11 = 1;
    a22 = 1;

    changeDiff = changeEnd - changeStart;
    
    // sinc
    bmax = M_PI * sIV.alphaBand;
    if (sIV.sincWidth == -1)
        sincWidth = numFromRightBound + 1;
    
    iLen = sincWidth * 2 + (even ? 1 : 0);
    xUMp1.resize (iLen, 0);
    
    bU = Eigen::VectorXd (iLen, 1);
    aW = Eigen::VectorXd (iLen, 1);

    AU = Eigen::MatrixXd (iLen, iLen);
    
    // file writing
    stateAt.open ("stateAtD.csv");
    plotIdx.open ("plotIdxD.csv");
    output.open ("outputD.csv");
    cSave.open ("cSaveD.csv");
    NSave.open ("NSaveD.csv");
    NChange.open ("NChangeD.csv");
    lambdaSqSave.open ("lambdaSqSaveD.csv");
    alfTickSave.open ("alfTickSaveD.csv");
    aUSave.open ("aUSaveD.csv");

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
    aUSave.close();
}

void OneDWaveDynamic::recalculateCoeffs (int n)
{
    if (changeC)
    {
        if (LFO)
            c = startC + (cos (2.0 * M_PI * n * LFOFreq / fs + M_PI) + 1) * 0.5 * (endC - startC);
        else if (n > changeStart && n <= changeEnd)
            c = startC - (n-changeStart) / static_cast<double> (changeDiff-1) * cDiff;
        h = c * k;
        NDouble = 1.0 / h;
        N = floor (NDouble);
    } else {
        if (LFO)
            NDouble = startN + (cos (2.0 * M_PI * n * LFOFreq / fs) + 1) * 0.5 * (endN - startN);
        else if (n >= changeStart && n <= changeEnd) // increase in N
            NDouble = startNTrue + (n-changeStart) / static_cast<double> (changeDiff) * NDiff;
        else
        {
            h = c * k;
            NDouble = 1.0 / h;
        }
        N = floor (NDouble);
        h = 1.0 / NDouble;
        c = h / k;
    }
//    h = c * k;
//    NDouble = 1.0 / h;
//    N = floor (NDouble);
    
    lambdaSq = lambdaMultiplier * c * c * k * k / (h * h);
    lambdaSqSave << lambdaSq << ";\n";
    NSave << NDouble << ";\n";
    
    alf = NDouble - N;
    alfTickSave << alf << ";\n"; // using for alf instead of alfTick now
    if (N != NPrev)
    {
        NChange << NPrev << ", " << N << ";\n";
        if (abs(N - NPrev) > 1)
            std::cout << "too fast..?" << std::endl;
        alfTick = ((1.0-Mw * h) - ((Mu + 1) * h)) / h;
//        alfTickSave << alfTick << ";\n";
        if (round(round(alfTick * 10) / 15.0) == 1)
            std::cout << alfTick - alf << std::endl;
        
        if (N > NPrev)
        {
            createCustomIp();
        }
        
        // call function
        (this->*addRemovePoints)();
    }
//    if (n % 441 == 0)

//    if (n > 10)
//    {
//        std::cout << std::endl;
//        stopSimulation = true;
//    }    
    NPrev = N;
    
    calculateInterpolatedPoints();
    
    B0 = 2.0 - 2.0 * lambdaSq;
    B1 = lambdaSq;
    C = -1.0;
}


void OneDWaveDynamic::calculateInterpolatedPoints()
{
    double w1 = 0;
    double w2 = 0;
    
    switch (dyIntType)
    {
        case dSinc:
        {
            if (sIV.sincWidth == -1) // make use of the full range
            {
                if (numFromRightBound == -1)
                {
                    sincWidth = floor (N * 0.5) - 1;
                    iLen = sincWidth * 2 + 1;
                    xUMp1.resize (iLen, 0);
                    bU = Eigen::VectorXd (iLen, 1);
                    aW = Eigen::VectorXd (iLen, 1);
                    AU = Eigen::MatrixXd (iLen, iLen);
                
                } else {
                    sincWidth = numFromRightBound + 1;
                }
            }
            bU.setZero();
            aU.setZero();
            aW.setZero();
            AU.setZero();
            
            // set alpha to 1e-6, otherwise the system can't solve
            if (alf < 1e-6)
                alf = 1e-6;
            

            for (int i = 0; i < sincWidth - (even ? 0 : 1); ++i)
                xUMp1[i] = i-sincWidth + (even ? 0 : 1);
            for (int i = sincWidth - (even ? 0 : 1); i < iLen; ++i)
                xUMp1[i] = i - sincWidth + alf - (even ? 1 : 0);
            
            for (int i = 0; i < iLen; ++i)
            {
                if (xUMp1[i] == 0)
                    bU[i] = 1;
                else
                    bU[i] = sin(bmax*xUMp1[i]) / xUMp1[i];
            }
            double dist = 0;
            for (int i = 0; i < iLen; ++i)
            {
                for (int j = 0; j < iLen; ++j)
                {
                    if (i == j)
                    {
                        AU(i, j) = bmax;
                        continue;
                    }
           
                    dist = xUMp1[i] - xUMp1[j];
                    AU(i, j) = sin (bmax*dist) / dist;
                }
            }
            aU = AU.bdcSvd (Eigen::ComputeThinU | Eigen::ComputeThinV).solve(bU); // aU = AU\bU
//            for (int i = 0; i < iLen; ++i)
//                aUSave << aU[i] << ", ";
//            aUSave << ";\n";

            for (int i = 0; i < iLen; ++i) // fliplr
                aW[iLen - 1 - i] = aU[i];
            
            // Calculate points
            uMuP1 = 0;
            wMin1 = 0;
            int idx = 0;
            int rangeEnd = iLen;
            bool useBoundary = false;
            
            // set range in different cases
            if (Mu + sincWidth == N+1) // use boundary condition
            {
                rangeEnd = iLen - 2;
                useBoundary = true;
            }
            else if (Mu + sincWidth == N) // use boundary condition
            {
                rangeEnd = iLen - 1;
            } else {
                rangeEnd = iLen;
            }
            
            // Calculate virtual grid points
            for (int i = 0; i < rangeEnd; ++i)
            {
                idx = Mu-sincWidth + i + (even || shifted ? 0 : 1);
                if (idx < Mu)
                    uMuP1 += aU[i] * u[1][idx];
                else
                    uMuP1 += aU[i] * w[1][idx-Mu];
            }
            
            if (useBoundary)
                uMuP1 -= aU[iLen-1] * w[1][Mw - 1];
            
            for (int i = 0; i < fmin (rangeEnd + 1, iLen); ++i)
            {
                idx = Mu-sincWidth-(shifted ? 0 : 1) + i;
                if (idx < Mu)
                    wMin1 += aW[i] * u[1][idx];
                else
                    wMin1 += aW[i] * w[1][idx-Mu];
            }
            break;
        }
        case dQuadratic:
        {
            alpha = (alf - 1.0) / (alf + 1.0);
            beta = 1.0;
            gamma = -(alf - 1.0) / (alf + 1.0);
            
            if (numFromRightBound != 1)
                w1 = w[1][1];
            
            uMuP1 = alpha * u[1][Mu-1] + beta * w[1][0] + gamma * w1;
            wMin1 = gamma * u[1][Mu-2] + beta * u[1][Mu-1] + alpha * w[1][0];

            break;
        }
        case dCubic:
        {
            alpha = alf * (alf - 1) * (alf - 2) / -6.0;
            beta = (alf - 1) * (alf + 1) * (alf - 2) * 0.5;
            gamma = alf * (alf + 1) * (alf - 2) * -0.5;
            delta = alf * (alf + 1) * (alf - 1) / 6.0;
            
            a21 = -delta;
            a12 = -delta;
            
            oOdet = 1.0 / (a11 * a22 - a12 * a21);
            
            if (numFromRightBound > 2 || numFromRightBound == -1)
            {
                w2 = w[1][2];
                w1 = w[1][1];
            }
            else if (numFromRightBound == 2)
            {
                w1 = w[1][1];
            }
            else if (numFromRightBound == 1)
            {
                w2 = -w[1][0];
            }
            
            v1 = alpha * w2 + beta * w1 + gamma * w[1][0];
            v2 = alpha * u[1][Mu-3] + beta * u[1][Mu-2] + gamma * u[1][Mu-1];
            
            uMuP1 = (v1 * a22 - v2 * a12) * oOdet;
            wMin1 = (-v1 * a21 + v2 * a11) * oOdet;
            break;
        }
        case dAltCubic:
        {
            if (shifted)
            {
                alpha = -(alf*(alf - 1.0))/((alf + 1.0)*(alf + 2.0));
                beta = (2.0*(alf - 1.0))/(alf + 1.0);
                gamma = 2.0/(alf + 1.0);
                delta = -(2.0*(alf - 1.0))/((alf + 1.0)*(alf + 2.0));
            } else {
                alpha = (alf - 1.0) / (alf + 2.0);
                beta =  (alf + 1.0) / 2.0;
                gamma = 1.0 - alf;
                delta =  alf * (alf - 1.0) / (2.0 * (alf + 2.0));
            }
            if (numFromRightBound > 2 || numFromRightBound == -1)
            {
                w2 = w[1][2];
                w1 = w[1][1];
            }
            else if (numFromRightBound == 2)
            {
                w1 = w[1][1];
            }
            else if (numFromRightBound == 1)
            {
                w2 = -w[1][0];
            }
            
            if (shifted) // use the same points for both virtual grid points
            {
                uMuP1 = alpha * u[1][Mu-2] + beta * u[1][Mu-1] + gamma * w[1][0] + delta * w1;
                wMin1 = delta * u[1][Mu-2] + gamma * u[1][Mu-1] + beta * w[1][0] + alpha * w1;
            } else {
                uMuP1 = alpha * u[1][Mu-1] + beta * w[1][0] + gamma * w1 + delta * w2;
                wMin1 = delta * u[1][Mu-3] + gamma * u[1][Mu-2] + beta * u[1][Mu-1] + alpha * w[1][0];
            }
            break;
        }
        case dQuartic:
        {
            alpha = -(alf*(alf - 1.0))/((alf + 2.0)*(alf + 3.0));
            beta =  (2.0*(alf - 1.0))/(alf + 2.0);
            gamma = 1.0;
            delta =  -(2.0*(alf - 1.0))/(alf + 2.0);
            epsilon = (alf*(alf - 1.0))/((alf + 2.0)*(alf + 3.0));
            
            if (numFromRightBound > 2 || numFromRightBound == -1)
            {
                w2 = w[1][2];
                w1 = w[1][1];
            }
            else if (numFromRightBound == 2)
            {
                w1 = w[1][1];
            }
            else if (numFromRightBound == 1)
            {
                w2 = -w[1][0];
            }
            
            uMuP1 = alpha * u[1][Mu-2] + beta * u[1][Mu-1] + gamma * w[1][0] + delta * w1 + epsilon * w2;
            wMin1 = epsilon * u[1][Mu-3] + delta * u[1][Mu-2] + gamma * u[1][Mu-1] + beta * w[1][0] + alpha * w1;
            break;
        }
        case dLinear:
        {
            if (numFromRightBound != 1)
            {
                w1 = (1.0 - alf) * w[1][1];
            }
            uMuP1 = w1 + alf * w[1][0];
            wMin1 = (1.0 - alf) * u[1][Mu-2] + alf * u[1][Mu-1];
            break;
        }
    }
}

void OneDWaveDynamic::scheme()
{
    double w1 = 0;
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
    if (numFromRightBound != 1)
        w1 = w[1][1];
    
    w[0][0] = B0 * w[1][0] + B1 * (w1 + wMin1) + C * w[2][0];
    if (numFromRightBound != 1)
        w[0][Mw-1] = B0 * w[1][Mw-1] + B1 * w[1][Mw-2] + C * w[2][Mw-1];

    if (lpConnection)
        lpPoints();
    
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
    if (fromRightBoundary)
    {
        if (Mw-indexFromBoundary < 0)
            output << u[1][Mu-indexFromBoundary + Mw] << ";\n";
        else
            output << w[1][Mw-indexFromBoundary] << ";\n";
    //    output << w[1][0] << ";\n";
    } else {
        output << u[1][indexFromBoundary - 1] << ";\n";
    }
}

void OneDWaveDynamic::retrieveStateU()
{
    plotIdx << curPlotIdx << ";\n";
    for (int l = 0; l < Mu; ++l)
        stateAt << std::to_string(u[1][l]) << ";\n";
    curPlotIdx += Mu-1;
    plotIdx << curPlotIdx << ";\n";
    ++curPlotIdx;
    cSave << c << ";\n";
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

void OneDWaveDynamic::addRemoveInCenter()
{
    if (N > NPrev) // add point
    {
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
            // save w0 (and prev) beforehand, otherwise things will be overwritten
            double w0 = customIp[3] * u[1][Mu-2]
                + customIp[2] * u[1][Mu-1]
                + customIp[1] * w[1][0]
                + customIp[0] * w[1][1];
            double w0Prev = customIp[3] * u[2][Mu-2]
                + customIp[2] * u[2][Mu-1]
                + customIp[1] * w[2][0]
                + customIp[0] * w[2][1];

            // move w vector one up (can be optimised)
            for (int l = Mw-1; l >= 0; --l)
            {
                w[1][l+1] = w[1][l];
                w[2][l+1] = w[2][l];
                
            }
            w[1][0] = w0;
            w[2][0] = w0Prev;
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

void OneDWaveDynamic::addRemoveAtPoint()
{
    double w1 = 0;
    double w1Prev = 0;
    if (N > NPrev) // add point
    {
        if (numFromRightBound != 1)
        {
            w1 = w[1][1];
            w1Prev = w[2][1];
        }
    
        u[1][Mu] = customIp[0] * u[1][Mu-2]
        + customIp[1] * u[1][Mu-1]
        + customIp[2] * w[1][0]
        + customIp[3] * w1;
        u[2][Mu] = customIp[0] * u[2][Mu-2]
        + customIp[1] * u[2][Mu-1]
        + customIp[2] * w[2][0]
        + customIp[3] * w1Prev;
        ++Mu;
        
    } else {
        u[1][Mu-1] = 0;
        u[2][Mu-1] = 0;
//        double avg_uMw0 = (u[1][Mu-2] + w[1][0]) * 0.5;
//        double avg_uMw0m1 = (u[2][Mu-2] + w[2][0]) * 0.5;
//        u[1][Mu-2] = (1 - alf) * avg_uMw0 + alf * u[1][Mu-2];
//        u[2][Mu-2] = (1 - alf) * avg_uMw0m1 + alf * u[2][Mu-2];
//        w[1][0] = (1 - alf) * avg_uMw0 + alf * w[1][0];
//        w[2][0] = (1 - alf) * avg_uMw0m1 + alf * w[2][0];
        --Mu;
    }
}

void OneDWaveDynamic::createCustomIp()
{
    switch (dyIntType)
    {
        case dLinear:
        {
            customIp[0] = 0;
            customIp[1] = alfTick / (alfTick + 1.0);
            customIp[2] = 1.0 / (alfTick + 1.0);
            customIp[3] = 0;
            break;
        }
        case dQuadratic:
        case dCubic:
        case dAltCubic:
        case dQuartic:
//        case dSinc:
        {
            customIp[0] = -alfTick * (alfTick + 1.0) / ((alfTick + 2.0) * (alfTick + 3.0));
            customIp[1] = 2.0 * alfTick / (alfTick + 2.0);
            customIp[2] = 2.0 / (alfTick + 2.0);
            customIp[3] = -2.0 * alfTick / ((alfTick + 3.0) * (alfTick + 2.0));
//            if (alf-alfTick != 0)
//                std::cout << "alf: " << alf << " alfTick " << alfTick << std::endl;
            //            customIp1[3] = alf * (alf - 1) * (alf - 2) / -6.0;
            //            customIp1[2] = (alf - 1) * (alf + 1) * (alf - 2) / 2.0;
            //            customIp1[1] = alf * (alf + 1) * (alf - 2) / -2.0;
            //            customIp1[0] = alf * (alf + 1) * (alf - 1) / 6.0;
            break;
        }
            
        case dSinc:
        {
            int custSincWidth = 2;
            double custBMax = M_PI;
            std::vector <double> xPosCustIp (4, 0);

            for (int i = 0; i < custSincWidth; ++i)
                xPosCustIp[i] = i-custSincWidth;
            for (int i = custSincWidth; i < custSincWidth * 2; ++i)
                xPosCustIp[i] = i - custSincWidth + alf;

            for (int i = 0; i < custSincWidth * 2; ++i)
                customIp[i] = sin(custBMax * xPosCustIp[i]) / (xPosCustIp[i] * custBMax);

            if (alf == 0)
                customIp[custSincWidth] = 1;
            break;
        }
    }
}

void OneDWaveDynamic::lpPoints()
{
//    double diffAtConn = w[1][0] - u[1][Mu-1];
//    double lpCoeff = pow (1-alf, lpExponent);
//    u[1][Mu-1] = u[1][Mu-1] + lpCoeff * diffAtConn * 0.5;
//    w[1][0] = w[1][0] - lpCoeff * diffAtConn * 0.5;
    
    double etaPrev = (w[2][0] - u[2][Mu-1]) * 0.5;
    double sig0 = 2.0;
    double rForce = (1.0 - sig0 / k) / (1.0 + sig0 / k);
    double oOP = (h * (1.0 + sig0 / k) * (1.0-alf)) / (2.0 * h * alf + k * k * (1.0 + sig0 / k) * (1.0-alf));
    
    double F = ((w[0][0] - u[0][Mu-1]) * 0.5 + rForce * etaPrev) * oOP;
    
    u[0][Mu-1] += k*k/h * F;
    w[0][0] -= k*k/h * F;

}
