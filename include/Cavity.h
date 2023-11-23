#ifndef CAVITY_H 
#define CAVITY_H

#include <vector>
#include "InputData.h"
#include "Beam.h"
#include <complex>
class Cavity {
public:
    CavityProperties cavityProps;
    int iTurn; // the turn number of the cavity.
    // we have two sets of time point that we need to keep track of.
    std::vector<double> time;     // time points for each particle.
    std::vector<double> feed_time;  // time points for feedback.
    // For particle, we need to keep track of the generator voltage and beam voltage at each arrival time point.
    std::vector<std::complex<double>> vGen;    // generator voltage at the arrival time point of each particle.
    std::vector<std::complex<double>> vBeam;   // beam voltage at the arrival time point of each particle.
    std::vector<std::complex<double>> vRef;    // reference voltage at each feedback time point.
    
    // for feedback, we only need to keep record of the key values at the last feedback time point.
    int n_feed;  // number of feedback time points.
    double tlast_feed;  // the time of the last feedback.
    double tlast_Vb; // the time of the last beam voltage update.
    std::complex<double> Ig; // generator current at last feedback time point.
    // these two variables are use to calculate the genrator voltage for every time point after the last feedback time point.
    std::complex<double> Vgm1; // the generator voltage at the last time point where the value of Ig is updated.
    std::complex<double> Ugm1; // The integral of generator voltage at last time point where the value of Ig is updated.
    // Then, at each feedback point, we need to calculate the generator and beam voltage and figure out the feedback we need to apply to the generator current.
    std::complex<double> Vg; // the generator voltage at the current feedback time point.
    std::complex<double> Vgp; // the derivative of the generator voltage at the current feedback time point.
    std::complex<double> Ug; // the integral of the generator voltage at the current feedback time point.
    std::complex<double> Vb; // the beam voltage at the current feedback time point.

    Cavity(const InputData& inputData, const Beam& beam);

    void printCavityProps() const;
    // these two functions are for the beam.
    void calVgen(Beam& beam, int& istart, int& iend,int& ilast_feed);
    void calVbeam(Beam& beam, int& istart, int& iend,int& ilast_feed);
    // these are for feedback.
    // we can reuse the calVgen but to me putting these two functionalities 
    // this way is more clear.
    void updateVadd(Beam& beam);
    void updateVg(const int& ifeed);
    void updateVb(const int& ifeed);
    void updateIg(const int& ifeed);
    void updateLastTimes(const int& ifeed);
    // the feedback function should update the generator current and Vg Vgm Ug Ugm and Vb. 
    void feedback(const int& ifeed);

    // kick the particles 
    void kickPar(Beam& beam);

    void findIstartAndIend(const int& ifeed, int& istart, int& iend);
    // full map for the cavity
    void fullMap(Beam& beam);
};

#endif // CAVITY.H