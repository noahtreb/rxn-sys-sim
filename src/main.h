#ifndef MAIN_H
#define MAIN_H

#include <string>
#include <utility>

class System;
class FileInterface;
class Distribution;

bool pairSort(std::pair<int, int> left, std::pair<int, int> right);
int simFwd(System* sys, int startTimePt, int endTimePt, double* time, double** state, Distribution** dists, int* speciesDistKey, int boundHandlingMethod, int& boundBreachSpeciesId, int& breachTimePt);
int simRev(System* sys, int startTimePt, int endTimePt, double* time, double** state, Distribution** dists, int* speciesDistKey, int boundHandlingMethod);
void writeStateData(int dataSavePtId, System* sys, FileInterface* fi, std::string varName, int trial, double** state, int startTimePt, int numTimePts, omp_lock_t& fileLock);       
void averageStateData(int numTrials, int numTimePts, int numSpecies, double*** stateData, double** avStateData, double** lastAvStateData);
double calcDist(double** avStateData1, double** avStateData2);

#endif