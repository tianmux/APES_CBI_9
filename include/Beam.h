#ifndef BEAM_H
#define BEAM_H

#include <vector>
#include "InputData.h"

class Beam {
public:
    int nPar;
    double qPb;
    std::vector<int> pattern;
    int fill_step;
    int nBunch;
    int nTrain;
    int dynamicOn;
    int parStoreStep;
    int n_parStore;
    double Ek;
    double E0;
    double Vrad;
    double Gamma0;
    double GMTSQ;
    double tau;
    double sig_gamma;
    double sig_time;

    std::vector<double> time; // time coordinates of each particle
    std::vector<double> gamma; // gamma coordinates of each particle
    std::vector<double> q; // charge of each particle
    std::vector<double> centroid_time; // centroid time of each bunch
    std::vector<double> bucket_centers; // bucket centers of each bunch
    std::vector<double> centroid_gamma; // centroid gamma of each bunch
    std::vector<std::vector<double>> t_Cen_Hist; // time centroid history of the particles
    std::vector<std::vector<double>> g_Cen_Hist; // gamma centroid history of the particles
    std::vector<double> vAdd; // voltage from the single macro particle
    
    // Statistical properties
    double meanTime, meanGamma;
    double varianceTime, varianceGamma;
    // ... other statistical properties as needed

    Beam(const InputData& inputData);

    void updateTimes(const std::vector<double>& newTimes);
    void updateGammas(const std::vector<double>& newGammas);
    void updateCharges(const std::vector<double>& newCharges);

    void calculateStatistics();

    // Getters for statistical moments
    void calcualteCentroids();
    void getVarianceTime() const;
    void getVarianceGamma() const;
    // ... other statistical measures as needed

    void computeMean();
    void computeVariance();
    // ... other computation methods as needed

    void storeParCentroids(int i);
    void printBeamStatistics() const;
    void printCentroids() const;
    void printVoltages() const;
    // output 
    void writeCentroidsToFile(std::string filename, int istart, int iend);
};

#endif // BEAM_H
