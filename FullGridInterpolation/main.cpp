//
//  main.cpp
//  OneDWave
//
//  Created by Silvin Willemsen on 31/10/2020.
//  Copyright © 2020 Silvin Willemsen. All rights reserved.
//

#include "OneDWave.h"
#include "OneDWaveDynamic.h"

#include <iostream>
#include <fstream>

void loop (OneDWave& oneDWave, int start, int end, int lengthSound, double fs, int outputLocFromRightBoundary)
{
    int curPercentage = start / static_cast<double> (lengthSound - 1);
    double test = 0;
    for (int n = start; n < end; ++n)
    {
        
        oneDWave.recalculateCoeffs (n);
        oneDWave.scheme();
        oneDWave.retrieveOutput (fs / 44100.0 * outputLocFromRightBoundary, false);
        if (n % int (fs / 44100.0) == 0 && n < 4000)
            oneDWave.retrieveState (-1);
        oneDWave.updateStates();
        
        test = n * 100 / static_cast<double>(lengthSound - 1);
        if (test > curPercentage)
        {
            std::cout << curPercentage << "% done" << std::endl;
            ++curPercentage;
        }
    }
}

void loop (OneDWaveDynamic& oneDWave, int start, int end, int lengthSound, double fs, int outputLocFromBoundary)
{
    int curPercentage = start / static_cast<double> (lengthSound - 1);
    double test = 0;
    int counter = 0;
    for (int n = start; n < end; ++n)
    {
        if (oneDWave.simulationStopped())
            return;
        oneDWave.recalculateCoeffs (n);
        oneDWave.scheme();
        oneDWave.retrieveOutput (fs / 44100.0 * outputLocFromBoundary, false);
//        if (int(n / 100.0) % int (fs) == 0)
        if (n < 1000)
        {
            oneDWave.retrieveStateU();
            oneDWave.retrieveStateW();
            ++counter;
//            std::cout << counter << std::endl;
        }
        oneDWave.updateStates();
        
        test = n * 100 / static_cast<double>(lengthSound - 1);
        if (test > curPercentage)
        {
            std::cout << curPercentage << "% done" <<  std::endl;
            ++curPercentage;
        }
        
    }
}

int adder(int a , int b){
    return a + b;
}

int main(int argc, const char * argv[]) {
    // insert code here...

#if DEBUG == 1
    std::cout << "I'm in DEBUG!" << std::endl;
#endif
    double fs = 44100;
    double outLength = 10;
    double lengthSound = fs * outLength;
    double NStart = 300;
    double NEnd = 15;
    double lambdaMultiplier = 1;
    int outputLocFromBoundary = 1;
    double excitationWidth = 0.25; //0.2
    double excitationLoc = 0.2; //1.0 / M_PI;
    double outputLocStart = 0;
    
    std::string version = "D"; // [I]nterpolation, [D]ynamic or [B]oth
    
    /*
         The numFromRightBound variable determines Number from the right boundary (quite important, switches between different techniques)
         -1: Adding to the center alternating between left and right string.
         0: Interpolated boundary
         1: Right string has a single moving point. Using simply supported boundary condition (if cubic)
         2: Right string has two moving points. When trying to solve the cubic
         interpolation, w_2 is always 0 (that's why this is a bit different)
         >3: (Expected behaviour) Selects where to add points (to left string).
     */
    
    int numFromRightBound = -1;
    
    // Set up different data
    std::ofstream curFs, curVersion;
    curFs.open ("curFs.csv");
    curFs << fs;
    curFs.close();
    
    curVersion.open ("curVersion.txt");
    curVersion << version;
    curVersion.close();
    
    SincInterpolVals sIV (2, 1); // sincWidth, alphaBandwidth. If sincWidth == -1, this results in using numFromRightBound + 1. If also numFromRightBound == -1, it uses the entire range.
    
    if (version == "I" || version == "B")
    {
        std::cout << "Interpolated version" << std::endl;
        OneDWave oneDWaveInterpol (NStart, NEnd,
                                   fs, outLength,
                                   excitationLoc, excitationWidth,
                                   outputLocStart, linear,
                                   lambdaMultiplier);
        loop (oneDWaveInterpol, 0, lengthSound, lengthSound, fs, outputLocFromBoundary);
        std::cout << "Done with interpolated version" << std::endl;

    }
    if (version == "D" || version == "B")
    {
        std::cout << "Dynamic version" << std::endl;
        OneDWaveDynamic oneDWaveDynamic (NStart, NEnd,
                                         fs, outLength,
                                         excitationLoc, excitationWidth,
                                         outputLocStart,
                                         dQuadratic, false, false, // interpolation type, even, shifted?
                                         sIV,
                                         lambdaMultiplier,
                                         false, numFromRightBound, //change C (or N)
                                         true, 30, 0,  // Low-pass the connection? And with what exponent?, And after how many seconds?
                                         false, 5); // LFO wavespeed/number of points? Freq
//        ,0.5 * fs / lengthSound, (0.5 + 3*abs(NEnd - NStart) / fs) * fs / lengthSound);
        loop (oneDWaveDynamic, 0, lengthSound, lengthSound, fs, outputLocFromBoundary);
        std::cout << "Done with dynamic version" << std::endl;

    }
#if DEBUG == 1
    std::cout << "I'm in DEBUG!" << std::endl;
#endif
    
    return 0;
}
