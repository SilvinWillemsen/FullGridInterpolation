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
    for (int n = start; n < end; ++n)
    {
        if (oneDWave.simulationStopped())
            return;
        oneDWave.recalculateCoeffs (n);
        oneDWave.scheme();
        oneDWave.retrieveOutput (fs / 44100.0 * outputLocFromRightBoundary, true);
        if (n % int (fs / 44100.0) == 0 && n < 1000)
        {
            oneDWave.retrieveStateU();
            oneDWave.retrieveStateW();
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

int main(int argc, const char * argv[]) {
    // insert code here...
#if DEBUG == 1
    std::cout << "I'm in DEBUG!" << std::endl;
#endif
    double fs = 44100;
    double outLength = 10;
    double lengthSound = fs * outLength;
    double NStart = 35.0;
    double NEnd = 65.0;
    double lambdaMultiplier = 1;
    int outputLocFromRightBoundary = 2;
    double excitationWidth = 0.2;
    double excitationLoc = 0.2;
    double outputLocStart = 0.9;

    std::string version = "D"; // [I]nterpolation, [D]ynamic
    
    std::ofstream curFs, curVersion;
    curFs.open ("curFs.csv");
    curFs << fs;
    curFs.close();
    
    curVersion.open ("curVersion.txt");
    curVersion << version;
    curVersion.close();
    
    if (version == "I")
    {
        std::cout << "Interpolated version" << std::endl;
        OneDWave oneDWaveInterpol (NStart, NEnd, fs, outLength, excitationWidth, excitationLoc, outputLocStart, cubic, lambdaMultiplier);
        loop (oneDWaveInterpol, 0, lengthSound, lengthSound, fs, outputLocFromRightBoundary);
    }
    else if (version == "D")
    {
        std::cout << "Dynamic version" << std::endl;
        OneDWaveDynamic oneDWaveDynamic (NStart, NEnd, fs, outLength, excitationWidth, excitationLoc, outputLocStart, lambdaMultiplier);
        loop (oneDWaveDynamic, 0, lengthSound, lengthSound, fs,  outputLocFromRightBoundary);
    }
#if DEBUG == 1
    std::cout << "I'm in DEBUG!" << std::endl;
#endif
    
    return 0;
}
