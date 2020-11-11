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
        oneDWave.retrieveOutput (fs / 44100.0 * outputLocFromRightBoundary, true);
        if (n % int (fs / 44100.0) == 0 && n < 1000)
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

void loop (OneDWaveDynamic& oneDWave, int start, int end, int lengthSound, double fs, int outputLocFromRightBoundary)
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
        oneDWave.retrieveOutput (fs / 44100.0 * outputLocFromRightBoundary, true);
        if (n % int (fs) == 0)
        {
            oneDWave.retrieveStateU();
            oneDWave.retrieveStateW();
            ++counter;
            std::cout << counter << std::endl;
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
    double NStart = 15;
    double NEnd = 20;
    double lambdaMultiplier = 1;
    int outputLocFromRightBoundary = 2;
    double excitationWidth = 0.2;
    double excitationLoc = 1.0 / M_PI;
    double outputLocStart = 0.9;

    std::string version = "D"; // [I]nterpolation, [D]ynamic or [B]oth
    int whereToAddPoints = -1;
    
    std::ofstream curFs, curVersion;
    curFs.open ("curFs.csv");
    curFs << fs;
    curFs.close();
    
    curVersion.open ("curVersion.txt");
    curVersion << version;
    curVersion.close();
    
    if (version == "I" || version == "B")
    {
        std::cout << "Interpolated version" << std::endl;
        OneDWave oneDWaveInterpol (NStart, NEnd,
                                   fs, outLength,
                                   excitationLoc, excitationWidth,
                                   outputLocStart, cubic,
                                   lambdaMultiplier);
        loop (oneDWaveInterpol, 0, lengthSound, lengthSound, fs, outputLocFromRightBoundary);
        std::cout << "Done with interpolated version" << std::endl;

    }
    if (version == "D" || version == "B")
    {
        std::cout << "Dynamic version" << std::endl;
        OneDWaveDynamic oneDWaveDynamic (NStart, NEnd,
                                         fs, outLength,
                                         excitationLoc, excitationWidth,
                                         outputLocStart, dCubic,
                                         lambdaMultiplier,
                                         true, whereToAddPoints);
        loop (oneDWaveDynamic, 0, lengthSound, lengthSound, fs,  outputLocFromRightBoundary);
        std::cout << "Done with dynamic version" << std::endl;

    }
#if DEBUG == 1
    std::cout << "I'm in DEBUG!" << std::endl;
#endif
    
    return 0;
}
