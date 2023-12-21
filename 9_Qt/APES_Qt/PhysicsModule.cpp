#include "PhysicsModule.h"
#include "InputData.h"
#include "Beam.h"
#include "Cavity.h"
#include "Ring.h"

#include <iostream>
PhysicsModule::PhysicsModule() {
    // Initialization code
}

void PhysicsModule::oneTurn(InputData& inputData, Cavity& cavity, Ring& ring, Beam& beam) {
    // Code to perform a full turn tracking of the beam
    // First we need to update the vAdd based on the charge per macro particle incase there is 
    // a change in the charge per macro particle for example we need to artificially turn off some bunch to 
    // simulate the bunch gap in the on-axis injectrion operation.
    cavity.fullMap(beam);
    beam.calcualteCentroids();
    beam.storeParCentroids(cavity.iTurn-1);
    ring.oneTurn(inputData, cavity, beam);

}
