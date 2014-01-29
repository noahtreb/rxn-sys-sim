#ifndef MAIN_H
#define MAIN_H

#include <string>

class System;
class FileInterface;
class Distribution;

void simFwd(System* sys, int numTimePts, double* time, double** state, Distribution** dists, int* speciesDistKey);
void simRev(System* sys, int numTimePts, double* time, double** state, Distribution** dists, int* speciesDistKey);
void writeStateData(int dataSavePtId, System* sys, FileInterface* fi, std::string varName, int trial, double** state, int numTimePts, omp_lock_t& lock);       
void averageStateData(int numTrials, int numTimePts, int numSpecies, double*** stateData, double** avStateData, double** lastAvStateData);
double calcDist(double** avStateData1, double** avStateData2);

#endif