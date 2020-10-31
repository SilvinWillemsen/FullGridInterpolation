//
//  main.cpp
//  FullGridInterpolation
//
//  Created by Silvin Willemsen on 31/10/2020.
//  Copyright © 2020 Silvin Willemsen. All rights reserved.
//

#include "FullGridInterpolation.h"
#include <iostream>
#include <fstream>

int main(int argc, const char * argv[]) {
    // insert code here...
    
    FullGridInterpolation fgi (30, 30, 44100, 0.5, 0.25, 0.2, 0.9);
    
    fgi.calculate();
    fgi.updateStates();
    return 0;
}
