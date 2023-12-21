// InputData.cpp

#include "InputData.h"
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

InputData::InputData() {
    // Initialize CavityProperties with default values
    cavityProps.numCavities = 1; // Example default value
    cavityProps.h = 400;
    cavityProps.RoQ = 50;
    cavityProps.QL = 2e5;
    cavityProps.Vsynch = 1.2e4;
    cavityProps.Vquad = 1.5e6;
    cavityProps.df = 0.0;
    cavityProps.feed_step = 400;
    cavityProps.gp = 0.0;
    cavityProps.gi = 0.0;
    cavityProps.gd = 0.0;
    cavityProps.gc = 0.0;
    cavityProps.epsilon = 0.0;
    cavityProps.delay = 0.0;

    // Initialize BeamProperties with default values
    beamProps.qPb = 2e-9;
    beamProps.nTrain = 1;
    beamProps.pattern.clear(); // Clearing the vector
    beamProps.nBunch = 400;
    beamProps.fill_step = 1;
    beamProps.Ek = 1890189000;
    beamProps.gamma0 = 3.699e+03;
    beamProps.gMTSQ = 5.882e+01;
    beamProps.nPar = 1;
    beamProps.tau = 0.0222e11;
    beamProps.sig_gamma = 220.7440935054759e-22;
    beamProps.sig_time = 2.6e-111;
    beamProps.dynamicOn = 0;

    // Initialize LatticeProperties with default values
    latticeProps.R = 0.0;
    latticeProps.nRamp = 0;
    latticeProps.nTrack = 0;
    latticeProps.sampRate = 0;
    latticeProps.SampRate_V = 0;
}


void InputData::loadDataFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    nlohmann::json json;
    file >> json;

    std::cout<<"Reading Cavity Properties"<<std::endl;
    // Parsing CavityProperties
    auto& cavityJson = json["cavityProps"];
    cavityProps.numCavities = cavityJson["numCavities"];
    cavityProps.h = cavityJson["h"];
    cavityProps.RoQ = cavityJson["RoQ"];
    cavityProps.QL = cavityJson["QL"];
    cavityProps.Vsynch = cavityJson["Vsynch"];
    cavityProps.Vquad = cavityJson["Vquad"];
    cavityProps.df = cavityJson["df"];
    
    cavityProps.feed_step = cavityJson["feed_step"];
    cavityProps.gp = cavityJson["gp"];
    cavityProps.gi = cavityJson["gi"];
    cavityProps.gd = cavityJson["gd"];
    cavityProps.gc = cavityJson["gc"];
    cavityProps.epsilon = cavityJson["epsilon"];
    cavityProps.delay = cavityJson["delay"];

    std::cout<<"Reading Beam Properties"<<std::endl;
    // Parsing BeamProperties
    auto& beamJson = json["beamProps"];
    beamProps.qPb = beamJson["qPb"];
    beamProps.nTrain = beamJson["nTrain"];
    beamProps.pattern.clear();
    for (const auto& val : beamJson["pattern"]) {
        beamProps.pattern.push_back(val);
    }
    beamProps.dynamicOn = beamJson["dynamicOn"];
    beamProps.parStoreStep = beamJson["parStoreStep"];

    beamProps.nBunch = beamJson["nBunch"];
    beamProps.fill_step = beamJson["fill_step"];
    beamProps.Ek = beamJson["Ek"];
    beamProps.Vrad = beamJson["Vrad"];
    beamProps.gamma0 = beamJson["gamma0"];
    beamProps.gMTSQ = beamJson["gMTSQ"];
    beamProps.nPar = beamJson["nPar"];
    beamProps.tau = beamJson["tau"];
    beamProps.sig_gamma = beamJson["sig_gamma"];
    beamProps.sig_time = beamJson["sig_time"];


    std::cout<<"Reading Latteice Properties"<<std::endl;
    // Parsing LatticeProperties
    auto& latticeJson = json["latticeProps"];
    latticeProps.R = latticeJson["R"];
    latticeProps.nRamp = latticeJson["nRamp"];
    latticeProps.nTrack = latticeJson["nTrack"];
    latticeProps.sampRate = latticeJson["sampRate"];
    latticeProps.SampRate_V = latticeJson["sampRate_V"];

    std::cout<<"Reading General Properties"<<std::endl;
    // Parsing GeneralProperties
    auto& generalJson = json["generalProps"];

    generalProps.t0 = generalJson["t0"];
    generalProps.n_thread = generalJson["n_thread"];
    generalProps.projectName = generalJson["projectName"];

    file.close();

    // Calculate derived quantities
    std::cout<<"Calculating Derived Quantities"<<std::endl;
    this->beta = sqrt(1 - 1 / pow(beamProps.gamma0, 2));
    this->tRev = 2 * M_PI * latticeProps.R / beta/c; 
    this->fRev = 1 / tRev;
    this->cavityProps.fRF = cavityProps.h * fRev;
    this->cavityProps.TRF = 1 / this->cavityProps.fRF;
    this->cavityProps.wRF = 2 * M_PI * this->cavityProps.fRF;
    this->damp_coeff = 2*this->tRev/this->beamProps.tau;
    this->excite_coeff = this->beamProps.sig_gamma*sqrt(3.0)*sqrt(4*this->tRev/this->beamProps.tau);
    this->cavityProps.fCav = this->cavityProps.fRF+this->cavityProps.df;
    this->cavityProps.wCav = 2*pi*this->cavityProps.fCav;
    this->cavityProps.TCav = 1/this->cavityProps.fCav;
    this->cavityProps.R = this->cavityProps.RoQ*this->cavityProps.QL;
    this->cavityProps.C = 1.0/this->cavityProps.RoQ/this->cavityProps.wCav;
    this->cavityProps.L = pow(this->cavityProps.RoQ,2)*this->cavityProps.C;
    this->cavityProps.alpha = this->cavityProps.wCav/2.0/this->cavityProps.QL;
    this->cavityProps.Z = this->cavityProps.RoQ*this->cavityProps.QL/(1.0+ii*this->cavityProps.QL*(this->cavityProps.fRF/this->cavityProps.fCav-this->cavityProps.fCav/this->cavityProps.fRF));
    this->cavityProps.IbDC = -this->beamProps.qPb*this->beamProps.nBunch/this->tRev;
    this->cavityProps.Vc = std::complex<double>(this->cavityProps.Vsynch,this->cavityProps.Vquad);
    this->cavityProps.Ig = (this->cavityProps.Vc/this->cavityProps.Z - this->cavityProps.IbDC*2.0);
    this->eta = 1/this->beamProps.gMTSQ - 1 / pow(this->beamProps.gamma0,2);
    this->Qs = sqrt(this->cavityProps.h*this->cavityProps.Vquad*this->eta/(2*pi*this->beamProps.Ek));
    this->cavityProps.phis = acos(-this->cavityProps.Vsynch/this->cavityProps.Vquad);
    this->df_opt = this->cavityProps.IbDC.real()*sin(this->cavityProps.phis)*this->cavityProps.RoQ/abs(this->cavityProps.Vc)*this->cavityProps.fRF;
    this->cavityProps.Vb_ref = this->cavityProps.IbDC*this->cavityProps.Z*2.0;
    std::cout<<"Done"<<std::endl;


}

void InputData::printInputData() const {

    // Print GeneralProperties
    std::cout << "\nGeneral Properties:" << std::endl;
    std::cout << "  Project Name: " << generalProps.projectName << std::endl;
    std::cout << "  Initial Time in unit of RF Cycle: " << generalProps.t0 << std::endl;
    std::cout << "  Number of Threads: " << generalProps.n_thread << std::endl;
    // Print CavityProperties
    std::cout << "Cavity Properties:" << std::endl;
    std::cout << "  Num Cavities: " << cavityProps.numCavities << std::endl;
    std::cout << "  Harmonic Number: " << cavityProps.h << std::endl;
    std::cout << "  R/Q: " << cavityProps.RoQ << std::endl;
    std::cout << "  Loaded Q: " << cavityProps.QL << std::endl;
    std::cout << "  Synchronous Voltage: " << cavityProps.Vsynch << std::endl;
    std::cout << "  Quadrature Voltage: " << cavityProps.Vquad << std::endl;
    std::cout << "  Frequency Shift: " << cavityProps.df << std::endl;
    std::cout << "  L: " << cavityProps.L << std::endl;
    std::cout << "  C: " << cavityProps.C << std::endl;
    std::cout << "  R: " << cavityProps.R << std::endl;
    std::cout << "  Alpha: " << cavityProps.alpha << std::endl;
    std::cout << "  wCav: " << cavityProps.wCav << std::endl;
    std::cout << "  TRF: " << cavityProps.TRF << std::endl;
    std::cout << "  TCav: " << cavityProps.TCav << std::endl;
    std::cout << "  fRF: " << cavityProps.fRF << std::endl;
    std::cout << "  fCav: " << cavityProps.fCav << std::endl;
    std::cout << "  phis: " << cavityProps.phis << std::endl;
    std::cout << "  Feedback Step: " << cavityProps.feed_step << std::endl;
    std::cout << "  Proportional Gain: " << cavityProps.gp << std::endl;
    std::cout << "  Integral Gain: " << cavityProps.gi << std::endl;
    std::cout << "  Derivative Gain: " << cavityProps.gd << std::endl;
    std::cout << "  Comb Gain: " << cavityProps.gc << std::endl;
    std::cout << "  Epsilon: " << cavityProps.epsilon << std::endl;
    std::cout << "  Feedback Delay: " << cavityProps.delay << std::endl;

    // Print BeamProperties
    std::cout << "\nBeam Properties:" << std::endl;
    std::cout << "  Charge per Bunch: " << beamProps.qPb << std::endl;
    std::cout << "  Number of Bunch Trains: " << beamProps.nTrain << std::endl;
    std::cout << "  Bunch pattern: ";
    for (const auto& val : beamProps.pattern) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    std::cout << "  Total Number of Bunches: " << beamProps.nBunch << std::endl;
    std::cout << "  Filling Step: " << beamProps.fill_step << std::endl;
    std::cout << "  Kinetic Energy: " << beamProps.Ek << std::endl;
    std::cout << "  Lorentz Factor: " << beamProps.gamma0 << std::endl;
    std::cout << "  Gamma_t Square: " << beamProps.gMTSQ << std::endl;
    std::cout << "  Number of Macro Particles per Bunch: " << beamProps.nPar << std::endl;
    std::cout << "  Longitudinal Damping Time: " << beamProps.tau << std::endl;
    std::cout << "  Relative RMS Energy Spread: " << beamProps.sig_gamma << std::endl;
    std::cout << "  RMS Bunch Length: " << beamProps.sig_time << std::endl;
    std::cout << "  Dynamic Simulation On: " << beamProps.dynamicOn << std::endl;
    std::cout << "  Particle Storage Step: " << beamProps.parStoreStep << std::endl;

    

    // Print LatticeProperties
    std::cout << "\nLattice Properties:" << std::endl;
    std::cout << "  Radius of the Ring: " << latticeProps.R << std::endl;
    std::cout << "  Number of Turns for Ramping: " << latticeProps.nRamp << std::endl;
    std::cout << "  Number of Turns for Tracking: " << latticeProps.nTrack << std::endl;
    std::cout << "  Sampling Rate: " << latticeProps.sampRate << std::endl;
    std::cout << "  Sampling Rate for Voltages: " << latticeProps.SampRate_V << std::endl;

    // Print Derived Quantities
    std::cout << "\nDerived Quantities:" << std::endl;
    std::cout << "  Beta: " << beta << std::endl;
    std::cout << "  Revolution Time: " << tRev << std::endl;
    std::cout << "  Revolution Frequency: " << fRev << std::endl;
    std::cout << "  Damping Coefficient: " << damp_coeff << std::endl;
    std::cout << "  Excitation Coefficient: " << excite_coeff << std::endl;
    std::cout << "  DC Beam Current: " << cavityProps.IbDC << std::endl;
    std::cout << "  Cavity Impedance: " << cavityProps.Z << std::endl;
    std::cout << "  Cavity Voltage: " << cavityProps.Vc << std::endl;
    std::cout << "  Generator Current: " << cavityProps.Ig << std::endl;
    std::cout << "  Synchrotron Tune: " << Qs << std::endl;
    std::cout << "  Optimal Frequency Shift: " << df_opt << std::endl;
    std::cout << "  Reference Beam Induced Voltage: " << cavityProps.Vb_ref << std::endl;
    std::cout << "  Momentum Compaction Factor: " << eta << std::endl;

}