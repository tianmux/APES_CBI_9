#include "Beam.h"
#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <iomanip>

// ... other includes ...

Beam::Beam(const InputData& inputData) {
    this->qPb = inputData.beamProps.qPb;
    this->nTrain = inputData.beamProps.nTrain;
    this->pattern = inputData.beamProps.pattern;
    this->nBunch = inputData.beamProps.nBunch;
    this->fill_step = inputData.beamProps.fill_step;
    this->Ek = inputData.beamProps.Ek;
    this->E0 = inputData.E0;
    this->Vrad = inputData.beamProps.Vrad;
    this->Gamma0 = inputData.beamProps.gamma0;
    this->GMTSQ = inputData.beamProps.gMTSQ;
    this->tau = inputData.beamProps.tau;
    this->sig_gamma = inputData.beamProps.sig_gamma;
    this->sig_time = inputData.beamProps.sig_time;
    this->nPar = inputData.beamProps.nPar;
    this->dynamicOn = inputData.beamProps.dynamicOn;
    this->parStoreStep = inputData.beamProps.parStoreStep;
    this->n_parStore = ceil(inputData.latticeProps.nTrack/this->parStoreStep);
    // We only allocate the memory for the vectors here
    // the values of the time vector can be updated later
    this->time.resize(nPar*nBunch);
    this->gamma.resize(nPar*nBunch);
    this->q.resize(nPar*nBunch);
    this->vAdd.resize(nPar*nBunch);
    
    this->t_Cen_Hist = std::vector<std::vector<double>>(n_parStore, std::vector<double>(nBunch));
    this->g_Cen_Hist = std::vector<std::vector<double>>(n_parStore, std::vector<double>(nBunch));

    this->centroid_time.resize(nBunch);
    this->bucket_centers.resize(nBunch);
    this->centroid_gamma.resize(nBunch);
    
    // ... other initializations ...

    // Get the initial centroids of each bunch
    int tmp_ibunch = 0;
    int tmp_ibucket = 0;
    for(int i = 0; i < nTrain; ++i){
        for (int j = 0; j < pattern[i*2]; ++j){
            this->centroid_time[tmp_ibunch] = inputData.generalProps.t0+inputData.cavityProps.TRF*tmp_ibucket;
            this->bucket_centers[tmp_ibunch] = this->centroid_time[tmp_ibunch];
            tmp_ibunch++;
            tmp_ibucket += inputData.beamProps.fill_step;
        }
        tmp_ibucket += inputData.beamProps.fill_step*this->pattern[i*2+1];
    }

    // Get the initial gamma of each bunch
    tmp_ibunch = 0;
    for(int i = 0; i < nTrain; ++i){
        for (int j = 0; j < pattern[i*2]; ++j){
            this->centroid_gamma[tmp_ibunch] = this->Gamma0;
            tmp_ibunch++;
        }
    }
    // Initialize the random number generator
    std::default_random_engine generator;
    // Initialize time and gamma of each macro particle
    tmp_ibunch = 0;
    tmp_ibucket = 0;
    for(int i = 0; i < nTrain; ++i){
        for (int j = 0; j < pattern[i*2]; ++j){
            for(int k = 0; k < nPar; ++k){
            // Initialize the normal distribution with mean and standard deviation
                std::normal_distribution<double> distribution(this->centroid_time[tmp_ibunch],this->sig_time);
                this->time[tmp_ibunch*nPar+k] = distribution(generator)+1e-12;
                distribution = std::normal_distribution<double>(this->centroid_gamma[tmp_ibunch],this->sig_gamma);
                this->gamma[tmp_ibunch*nPar+k] = distribution(generator);
            }
            tmp_ibunch++;
            tmp_ibucket += inputData.beamProps.fill_step;
        }
        tmp_ibucket += inputData.beamProps.fill_step*this->pattern[i*2+1];
    }
    // Initialize charge per macro particle 
    for(int i = 0; i < nBunch; ++i){
        for(int j = 0; j < nPar; ++j){
            this->q[i*nPar+j] = this->qPb/this->nPar;
        }
    }
}
void Beam::calcualteCentroids(){
    // Calculate the centroids of each bunch
    // ...
    for(int i = 0; i < this->nBunch; ++i){
        double tmp_time = 0;
        double tmp_gamma = 0;
        for(int j = 0; j < this->nPar; ++j){
            tmp_time += this->time[i*this->nPar+j];
            tmp_gamma += this->gamma[i*this->nPar+j];
        }
        this->centroid_time[i] = tmp_time/this->nPar;
        this->centroid_gamma[i] = tmp_gamma/this->nPar;
    }
    
}
void Beam::storeParCentroids(int i){
    for(int j = 0; j < this->nBunch; ++j){
        this->t_Cen_Hist[i][j] = this->centroid_time[j]-this->bucket_centers[j];
        this->g_Cen_Hist[i][j] = this->centroid_gamma[j]-this->Gamma0;
    }
}
void Beam::updateTimes(const std::vector<double>& newTimes){
    this->time = newTimes;
}

void Beam::updateGammas(const std::vector<double>& newGammas){
    this->gamma = newGammas;
}

void Beam::updateCharges(const std::vector<double>& newCharges){
    this->q = newCharges;
}

void Beam::calculateStatistics() {
    // Calculate the mean and variance of time and gamma
    // ...
}

void Beam::printBeamStatistics() const {
    // Print the mean and variance of time and gamma
    // ...
}

void Beam::printCentroids() const {
    // Print the centroids of each bunch
    // ...
    for(const auto& val : this->centroid_time){
        std::cout<<val<<std::endl;
    }
    for(const auto& val : this->centroid_gamma){
        std::cout<<val<<std::endl;
    }
}

void Beam::printVoltages() const {
    // Print the voltages
    // ...
    for(const auto& val : this->vAdd){
        std::cout<<val<<std::endl;
    }
}

void Beam::writeCentroidsToFile(std::string filename, int istart, int iend){
    std::ofstream outfile;
    outfile.open (filename);
    std::cout << std::fixed << std::setprecision(10);    
    for(int i = 0; i < this->n_parStore; ++i){
        for (int j = istart; j < iend; ++j){
            outfile << std::scientific << std::setprecision(10)<<this->t_Cen_Hist[i][j] << " " << std::setprecision(10)<<this->g_Cen_Hist[i][j] << "\n";
        }
        //outfile << "\n";
    }
    outfile.close();
}