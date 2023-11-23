#ifndef RING_H
#define RING_H

#include "InputData.h"
#include "Beam.h"
#include "Cavity.h"
#include <iostream>

class Ring{
public:
    double eta;
    double T0;
    Ring(InputData& inputData); // Constructor
    
    void oneTurn(InputData& inputData, Cavity& cavity, Beam& beam);

    //OutputData getResults();
};

#endif