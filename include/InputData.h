// InputData.h

#ifndef INPUTDATA_H
#define INPUTDATA_H

#include <vector>
#include <string>
#include <complex>


// Define the imaginary unit 'i'
const std::complex<double> ii(0.0, 1.0);


// Define additional structures if needed
struct CavityProperties {
    int numCavities; //number of cavities
    int h;   //hamonic number of the cavity
    double RoQ; // R over Q of the cavity
    double QL;  // Loaded Q of the cavity
    double Vsynch; // Synchronous voltage of the cavity
    double Vquad;   // Quadrature voltage of the cavity
    double df;  // Frequency shift of the resonance frequency of the cavity from the RF frequency
    double fRF; // RF frequency
    double TRF; // RF period
    double wRF; // RF angular frequency
    double fCav;    // cavity frequency
    double TCav;    // cavity period    
    double wCav;    // cavity angular frequency

    // derived quantities for circuit model
    double L;   // inductance of the cavity
    double C;   // capacitance of the cavity
    double R;   // resistance of the cavity
    double alpha;   // damping factor of the cavity
    std::complex<double> Z; // impedance of the cavity
    std::complex<double> Vc;    // cavity voltage
    std::complex<double> Ig;    // generator current
    std::complex<double> IbDC;  // DC current of the beam
    std::complex<double> Vb_ref; // reference voltage from the beam


    // Feedback parameters
    int feed_step; // feedback step in unit of bucket
    double gp;  // Proportional gain of the feedback
    double gi;  // Integral gain of the feedback
    double gd;  // Derivative gain of the feedback
    double gc;  // Gain of the comb
    double epsilon; // epsilon of the comb feedback
    double delay;   // delay of the feedback loop
    double phis;   // synchronous phase of the cavity
    // Add other cavity-specific properties
};

struct BeamProperties {
    std::vector<int>  pattern; // bunch pattern
    int nBunch; // total number of bunches 
    int fill_step; // filling step in unit of bucket
    int nPar;   // number of macro particles per bunch
    int dynamicOn;  // flag to turn on the dynamic simulation
    int parStoreStep;   // step to store the particle data
    double Ek;   // kinetic energy of the beam
    double Vrad; // radiation voltage
    double qPb; // charge per bunch
    double nTrain;  // number of bunch trains
    double gamma0;  // Lorentz factor of the beam
    double gMTSQ;   // gamma_t square
    double tau; // longitudinal damping time in unit of second
    double sig_gamma;   // relative rms energy spread
    double sig_time;    // rms bunch length in unit of second


    // Add other beam-specific properties
};

struct LatticeProperties {
    // Define Lattice or ring.
    double R;   // radius of the ring
    int nRamp;  // number of turns the code use to ramp the vb
    int nTrack; // number of turns the code use to track the beam
    int sampRate;   // sampling rate of the code for M1 and M2
    int SampRate_V; // sampling rate of the code for Voltages
    // Add other initial conditions
};

struct GeneralProperties{
    double t0;  // initial time in terms of number of RF cycles. 


};

// Main structure for input data
struct InputData {
    // Physical constants
    double c = 299792458;   // speed of light
    double e = 1.60217662e-19;  // elementary charge
    double m = 9.10938356e-31;  // electron mass
    double E0 = m*c*c/e;    // rest energy of electron
    double hbar = 1.0545718e-34;    // reduced Planck constant
    double pi = 3.14159265358979323846; // pi


    CavityProperties cavityProps;
    BeamProperties beamProps;
    LatticeProperties latticeProps;
    GeneralProperties generalProps;
    // Include other relevant parameters
    
    // Dirived parameters
    double beta;    // relativistic beta
    double tRev;    // revolution period
    double fRev;    // revolution frequency
    double damp_coeff;  // damping coefficient
    double excite_coeff;    // excitation coefficient
    double eta;    // momentum compaction factor
    double Qs;    // synchrotron tune
    double df_opt; // optimal frequency shift

    // Constructor for initializing with default values
    InputData(); 
    
    // Function to load data from a file or another source
    void loadDataFromFile(const std::string& filename);
    // Function to print input data
    void printInputData() const; 
};

#endif // INPUTDATA_H
