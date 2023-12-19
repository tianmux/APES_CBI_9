#include "InputData.h"
#include "Beam.h"
#include "Cavity.h"
#include "PhysicsModule.h"
#include <iostream>
#include <cassert>
#include <string>
#include <iomanip>
#include <chrono>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_input_file>" << std::endl;
        return 1;
    }
    std::cout << std::fixed << std::setprecision(10);    
    std::string outfile = "centroids.txt";

    std::string inputFilePath = argv[1];
    std::string printFlag = argv[2];
    std::cout << "Input file path: " << inputFilePath << std::endl;
    std::cout << "Print flag: " << printFlag << std::endl;

    InputData inputData;
    // Load data from a test JSON file
    try {
        inputData.loadDataFromFile(inputFilePath);
    } catch (const std::exception& e) {
        std::cerr << "Failed to load input data: " << e.what() << std::endl;
        return 1;
    }

    // Initialize the Beam object
    Beam beam(inputData);

    // Initialize the Cavity object
    Cavity cavity(inputData,beam);

    // Initialize the Ring object
    Ring ring(inputData);

    // Print the loaded data
    std::cout<< "Printing the loaded data" << std::endl;
    inputData.printInputData();
    PhysicsModule physicsModule;

    // Perform one turn tracking
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < inputData.latticeProps.nTrack; i++){
        if(i%int(inputData.latticeProps.nTrack/10)==0){
            end = std::chrono::high_resolution_clock::now();
            std::cout<<"Turn: "<< i << std::endl;
            std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;
        }
        physicsModule.oneTurn(inputData, cavity, ring, beam);    
        //beam.storeParCentroids(i);
        //cavity.printCavityProps();
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken to finish full simulation: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;
    // write the centroids to file
    beam.writeCentroidsToFile(inputData.generalProps.projectName+"_centroid.txt",0,1);


    if(printFlag == "1"){
        // Print the beam statistics
        std::cout<< "Printing the beam statistics" << std::endl;
        beam.printBeamStatistics();

        // Print the centroids of the beam
        std::cout<< "Printing the centroids of the beam" << std::endl;
        beam.printCentroids();

        // Print the voltages
        std::cout<< "Printing the voltages" << std::endl;
        beam.printVoltages();
    }

    

    std::cout << "All tests passed." << std::endl;
    return 0;
}
