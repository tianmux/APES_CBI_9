#include "Ring.h"
#include "InputData.h"
#include "Beam.h"
#include "Cavity.h"
#include <iostream>
#include <random>
#include <vector>

Ring::Ring(InputData& inputData) {
    // Initialization code
    std::cout<<"Constructing Ring"<<std::endl;
    this->eta = inputData.eta;
    this->T0 = inputData.tRev;
    std::cout<<" eta: "<< this->eta << std::endl;
    std::cout<<" T0: "<< this->T0 << std::endl;
    std::cout<< "Ring constructed" << std::endl;
}

void Ring::oneTurn(InputData& inputData, Cavity& cavity, Beam& beam){
    std::vector<double> delta;
    delta.resize(beam.nBunch*beam.nPar);

    // calculate the delta
    for(int i =0; i< beam.nBunch*beam.nPar;++i){
        delta[i] = beam.gamma[i]/beam.Gamma0 - 1;
    } 

    // update the time
    double dynamicOn = 0.0;
    if(cavity.iTurn>=beam.dynamicOn){
        dynamicOn=1.0;
    }
    for(int i = 0; i < beam.nBunch*beam.nPar; i++){
        beam.time[i] += this->T0*(1+this->eta*delta[i]*dynamicOn);
    }

    // update the feedback time point
    for(int i = 0;i < cavity.n_feed;++i){
        cavity.feed_time[i] += this->T0;
    }

    // update the bucket centers
    for(int i = 0; i < beam.nBunch; ++i){
        beam.bucket_centers[i] += this->T0;
    }
    // radiation effect and damping effect
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<> distr(-3.0, 1.0);

    const int numberOfElements = beam.nBunch*beam.nPar;
    std::vector<double> randomNumbers;

    // Reserve space for the elements
    randomNumbers.reserve(numberOfElements);

    for (int i = 0; i < numberOfElements; ++i) {
        randomNumbers.push_back(distr(eng));
    }

    for (int i = 0; i< beam.nBunch*beam.nPar; ++i){
        beam.gamma[i] -= beam.Vrad/beam.E0*dynamicOn;
        //beam.gamma[i] -= (beam.gamma[i]-beam.Gamma0)*inputData.damp_coeff;
        //beam.gamma[i] -= inputData.excite_coeff*randomNumbers[i];
    }
}