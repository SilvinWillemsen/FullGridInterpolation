//
//  main.cpp
//  OneDWave
//
//  Created by Silvin Willemsen on 31/10/2020.
//  Copyright © 2020 Silvin Willemsen. All rights reserved.
//

#include "OneDWave.h"
#include <iostream>
#include <fstream>

void loop (OneDWave& oneDWave, int start, int end, int lengthSound)
{
    int curPercentage = start / static_cast<double> (lengthSound - 1);
    double test = 0;
    for (int n = start; n < end; ++n)
    {
        
        
        oneDWave.scheme();
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
    
    double fs = 4410000;
    double outLength = 0.5;
    double lengthSound = fs * outLength;
    
    OneDWave oneDWave (30, 30, fs, outLength, 0.25, 0.2, 0.9);
    int sampleAtWhichToRetrieveState = 1000;
    
    loop (oneDWave, 0, sampleAtWhichToRetrieveState, lengthSound);
    oneDWave.retrieveState();
    loop (oneDWave, sampleAtWhichToRetrieveState, lengthSound, lengthSound);

    
    return 0;
}
