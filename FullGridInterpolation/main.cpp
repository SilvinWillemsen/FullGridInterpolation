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

void loop (OneDWave& oneDWave, int start, int end, int lengthSound, double fs)
{
    int curPercentage = start / static_cast<double> (lengthSound - 1);
    double test = 0;
    for (int n = start; n < end; ++n)
    {
        
        oneDWave.recalculateCoeffs (n);
        oneDWave.scheme();
        oneDWave.retrieveOutput (fs / 44100.0 * 5, true);
        oneDWave.updateStates();
        
        test = n * 100 / static_cast<double>(lengthSound - 1);
        if (test > curPercentage)
        {
            std::cout << curPercentage << std::endl;
            ++curPercentage;
        }
    }
}

void loop (OneDWaveDynamic& oneDWave, int start, int end, int lengthSound, double fs)
{
    int curPercentage = start / static_cast<double> (lengthSound - 1);
    double test = 0;
    for (int n = start; n < end; ++n)
    {
        if (oneDWave.simulationStopped())
            return;
        oneDWave.recalculateCoeffs (n);
        oneDWave.scheme();
        oneDWave.retrieveOutput (fs / 44100.0 * 5, true);
        oneDWave.updateStates();
        
        test = n * 100 / static_cast<double>(lengthSound - 1);
        if (test > curPercentage)
        {
            std::cout << curPercentage << std::endl;
            ++curPercentage;
        }
        
    }
}

int main(int argc, const char * argv[]) {
    // insert code here...
    
    double fs = 44100;
    double outLength = 3;
    double lengthSound = fs * outLength;
    double NStart = 30.0;
    double NEnd = 30.99;
    std::ofstream curFs;
    curFs.open ("curFs.csv");
    curFs << fs;
    curFs.close();
    OneDWave oneDWaveInterpol (NStart, NEnd, fs, outLength, 0.25, 0.2, 0.9, linear);
    
    OneDWaveDynamic oneDWaveDynamic (NStart, NEnd, fs, outLength, 0.25, 0.2, 0.9);

//    loop (oneDWaveInterpol, 0, lengthSound, lengthSound, fs);
    loop (oneDWaveDynamic, 0, lengthSound, lengthSound, fs);

    return 0;
}
