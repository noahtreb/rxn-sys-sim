#ifndef MAIN_H
#define MAIN_H

#include <mutex>
#include <random>
#include <string>
#include <vector>

class System;
class FileInterface;
class Distribution;

class RerunTrial {
public:
    int trialId;
    int startTimePt;
    int endTimePt;
    int boundBreachSpeciesId;
    
    int numPrevTrials;
    bool threadAssigned;
    int threadId;
    
    RerunTrial(int trialId, int startTimePt, int endTimePt, int boundBreachSpeciesId);
};

void rerunTrial(FileInterface* fi, System* sys, Distribution** revDists, Distribution* dist, std::vector<RerunTrial*>* rerunTrials, std::mt19937* rng, std::mutex* vecLock, std::mutex* fileLock, std::mutex* absCurrLock, std::string varName, double** boundStatePts, double** state, double* time, int** absCurr, int* speciesDistKey, int threadId, int boundHandlingMethod, int numTrials, int numTimePts);
bool rerunTrialSort(RerunTrial left, RerunTrial right);
bool simFwd(System* sys, int startTimePt, int endTimePt, double* time, double** state, Distribution** dists, int* speciesDistKey, int boundHandlingMethod, int& boundBreachSpeciesId, int& breachTimePt);
bool simRev(System* sys, int startTimePt, int endTimePt, double* time, double** state, Distribution** dists, int* speciesDistKey, int boundHandlingMethod);
void writeStateData(int dataSavePtId, System* sys, FileInterface* fi, std::string varName, int trial, double** state, int startTimePt, int numTimePts, omp_lock_t& fileLock);       
void writeStateData(int dataSavePtId, System* sys, FileInterface* fi, std::string varName, int trial, double** state, int startTimePt, int numTimePts, std::mutex* fileLock);
void averageStateData(int numTrials, int numTimePts, int numSpecies, double*** stateData, double** avStateData, double** lastAvStateData);
double calcDist(double** avStateData1, double** avStateData2);

#endif