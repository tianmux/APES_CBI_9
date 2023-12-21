#include "InputData.h"
#include "Beam.h"
#include "Cavity.h"
#include "Ring.h"
class PhysicsModule {
public:
    PhysicsModule(); // Constructor
    
    void oneTurn(InputData& inputData, Cavity& cavity, Ring& ring, Beam& beam);

    //OutputData getResults();

private:
    // Internal state variables
    //FieldData fieldData;
    //CavityProperties cavityProperties;
    //BeamProperties beamProperties;

    void calculateFieldInteractions();
    // Other private methods for internal calculations
};
