#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <float.h>
#include <random>
#include <string>
#include "main.h"
#include "Distribution.h"
#include "FileInterface.h"
#include "PriorityQueue.h"
#include "System.h"
#include "Reaction.h"
#include "Species.h"

using namespace std;

int main(int argc, const char* argv[]) {
    bool skipInitFwd = false;
    bool skipInitRev = false;
    bool skipRefine = false;
    int seed;
    double stoppingTol = -1;
    string fileName;
    
    if (argc >= 4) {
        for (int i = 0; i < argc; i++) {
            int temp = 0;
            if (argv[i][0] == '-') {
                switch(argv[i][1]) {
                    case 'r':
                        seed = atoi(argv[i+1]);
                        i++;
                        break;
                    case 's':
                        temp = atoi(argv[i+1]);
                        if (temp == 1) {
                            skipInitFwd = true;
                        }
                        i++;
                        break;
                    case 't':
                        temp = atoi(argv[i+1]);
                        if (temp == 1) {
                            skipInitRev = true;
                        }
                        i++;
                        break;
                    case 'u':
                        temp = atoi(argv[i+1]);
                        if (temp == 1) {
                            skipRefine = true;
                        }
                        i++;
                        break;
                    case 'e':
                        stoppingTol = atof(argv[i+1]);
                        i++;
                        break;
                }
            } else {
                fileName = argv[i];
            }
        }
    } else {
        fprintf(stderr, "Usage: ./RxnSysSim <fileName> -r <seed> [-s <skipInitFwd>] [-t <skipInitRev>] [-u <skipRefine>]\n");
        fprintf(stderr, "                   [-e <stoppingTol>] \n");
        fprintf(stderr, "            <fileName> Name of the netCDF file to read from and write to.\n");
        fprintf(stderr, "    -r          <seed> Integer used to seed the random number generator.\n");
        fprintf(stderr, "Optional arguments:\n");
        fprintf(stderr, "    -s   <skipInitFwd> Use this flag to indicate that preexisting data from the 'initFwdData' variable\n");
        fprintf(stderr, "                       should be used in the forward refinement algorithm.\n");
        fprintf(stderr, "    -t   <skipInitRev> Use this flag to indicate that preexisting data from the 'initRevData' variable\n");
        fprintf(stderr, "                       should be used in the reverse refinement algorithm.\n");
        fprintf(stderr, "    -u    <skipRefine> Use this flag to bypass the forward and reverse refinement algorithms.\n");
        fprintf(stderr, "    -e   <stoppingTol> Maximum allowable L2 distance between successive probability distributions.\n\n");
        abort();
    }
    
    FileInterface* fi = new FileInterface(fileName);
    std::mt19937 rng(seed);   
        
    if (stoppingTol > 0) {
        fi->overwriteStoppingTol(stoppingTol);
    }
    
    int numTrials;
    double startTime;
    double endTime;
    int numTimePts;
    double timeStep;
    int numDataSavePts;
    int* dataSavePts;
    int numBoundedSpeciesStates;
    
    System* masterSys = fi->readFileData(numTrials, startTime, endTime, numTimePts, timeStep, stoppingTol, numDataSavePts, dataSavePts, numBoundedSpeciesStates);
            
    double* time = new double[numTimePts];  
    
    int numThreads = omp_get_max_threads();
    
    double*** state = new double**[numThreads];
    for (int i = 0; i < numThreads; i++) {
        state[i] = new double*[masterSys->numSpecies];
        for (int j = 0; j < masterSys->numSpecies; j++) {
            state[i][j] = new double[numTimePts];
        }
    }      
    
    for (int i = 0; i < numTimePts; i++) {
        time[i] = startTime + timeStep * i;
    }    
    fi->writeTimeData(time, numTimePts);
            
    string varName;
    double** lastFwdStatePt = new double*[numTrials];
    double** lastRevStatePt = new double*[numTrials];
        
    for (int i = 0; i < numTrials; i++) {
        lastFwdStatePt[i] = new double[masterSys->numSpecies];
        lastRevStatePt[i] = new double[masterSys->numSpecies];
    }
    
    omp_lock_t lock;
    omp_init_lock(&lock);
        
    System* sys;
    
    if (skipRefine) {
        numDataSavePts = 0;
    }
    
    double fwdSetupStart[numDataSavePts + 1];
    double fwdSetupEnd[numDataSavePts + 1];
    double fwdTrialStart[numDataSavePts + 1][numTrials];
    double fwdWriteStart[numDataSavePts + 1][numTrials];
    double fwdTrialEnd[numDataSavePts + 1][numTrials];
    double revSetupStart[numDataSavePts + 1];
    double revSetupEnd[numDataSavePts + 1];
    double revTrialStart[numDataSavePts + 1][numTrials];
    double revWriteStart[numDataSavePts + 1][numTrials];
    double revTrialEnd[numDataSavePts + 1][numTrials];
    
    int j = 0;
    int* distSpeciesKey = new int[numBoundedSpeciesStates];
    int* speciesDistKey = new int[masterSys->numSpecies];
    for (int i = 0; i < masterSys->numSpecies; i++) {
        if (masterSys->species[i]->stateBounded) {
            distSpeciesKey[j] = i;
            speciesDistKey[i] = j;
            j++;
        } else {
            speciesDistKey[i] = -1;
        }
    }
    
    Distribution** revDistsPrev = new Distribution*[numBoundedSpeciesStates];
    for (int i = 0; i < numBoundedSpeciesStates; i++) {
        revDistsPrev[i] = new Distribution(masterSys->species[distSpeciesKey[i]], numTrials, rng());
    }
    
    Distribution** revDists = new Distribution*[numBoundedSpeciesStates];
    for (int i = 0; i < numBoundedSpeciesStates; i++) {
        revDists[i] = new Distribution(masterSys->species[distSpeciesKey[i]], numTrials, rng());
    }
    
    Distribution** fwdDistsPrev = new Distribution*[numBoundedSpeciesStates];
    for (int i = 0; i < numBoundedSpeciesStates; i++) {
        fwdDistsPrev[i] = new Distribution(masterSys->species[distSpeciesKey[i]], numTrials, rng());
    }
    
    Distribution** fwdDists = new Distribution*[numBoundedSpeciesStates];
    for (int i = 0; i < numBoundedSpeciesStates; i++) {
        fwdDists[i] = new Distribution(masterSys->species[distSpeciesKey[i]], numTrials, rng());
    }
       
    double progStart = omp_get_wtime();
    if (!skipInitFwd) {        
        fwdSetupStart[0] = omp_get_wtime();
        for (int i = 0; i < numBoundedSpeciesStates; i++) {
            revDists[i]->addNode(revDists[i]->species->state, numTrials);
        }
        fwdSetupEnd[0] = omp_get_wtime();
        
        #pragma omp parallel for default(shared) private(sys)
        for (int i = 0; i < numTrials; i++) {
            fwdTrialStart[0][i] = omp_get_wtime();
            int threadId = omp_get_thread_num();            
                    
            sys = new System(*masterSys);
            sys->seed(rng());
            sys->initFwd();
            
            varName = "initFwdData";
            simFwd(sys, numTimePts, time, state[threadId], revDists, speciesDistKey);
            
            fwdWriteStart[0][i] = omp_get_wtime();
            writeStateData(-1, sys, fi, varName, i, state[threadId], numTimePts, lock);
            fwdTrialEnd[0][i] = omp_get_wtime();
        }
    } else {
        fwdSetupStart[0] = 0;
        fwdSetupEnd[0] = 0;
        
        for (int i = 0; i < numTrials; i++) {
            fwdTrialStart[0][i] = 0;
            fwdWriteStart[0][i] = 0;
            fwdTrialEnd[0][i] = 0;
        }
    }
    
    if (!skipInitRev) {
        revSetupStart[0] = omp_get_wtime();
        
        masterSys->rxnPq->minHeap = false;
        masterSys->time = time[numTimePts - 1];
        fi->readInitDataPt("initFwdData", numTrials, masterSys->numSpecies, numTimePts - 1, lastFwdStatePt);
        
        for (int i = 0; i < numBoundedSpeciesStates; i++) {             
            fwdDists[i]->update(lastFwdStatePt, numTrials);
        }
                
        revSetupEnd[0] = omp_get_wtime();
        
        #pragma omp parallel for default(shared) private(sys)
        for (int i = 0; i < numTrials; i++) {    
            revTrialStart[0][i] = omp_get_wtime();
            int threadId = omp_get_thread_num(); 
            
            sys = new System(*masterSys);
            sys->seed(rng());            
            
            for (int j = 0; j < sys->numSpecies; j++) {
                sys->species[j]->state = lastFwdStatePt[i][j];
            }

            for (int j = 0; j < sys->numRxns; j++) {
                sys->rxns[j]->updateProp(sys->volRatio);
            }
            sys->initRev();

            varName = "initRevData";
            simRev(sys, numTimePts, time, state[threadId], fwdDists, speciesDistKey);
            
            revWriteStart[0][i] = omp_get_wtime();
            writeStateData(-1, sys, fi, varName, i, state[threadId], numTimePts, lock);
            revTrialEnd[0][i] = omp_get_wtime();
        }
    } else {
        revSetupStart[0] = 0;
        revSetupEnd[0] = 0;
        
        for (int i = 0; i < numTrials; i++) {
            revTrialStart[0][i] = 0;
            revWriteStart[0][i] = 0;
            revTrialEnd[0][i] = 0;
        }
    }  
    
    if (!skipRefine) {        
        for (int i = 0; i < numDataSavePts; i++) {
            fwdSetupStart[i + 1] = omp_get_wtime();
            
            masterSys->rxnPq->minHeap = true;
            masterSys->time = time[0];
            if (i == 0) {
                fi->readInitDataPt("initRevData", numTrials, masterSys->numSpecies, 0, lastRevStatePt);
            } else {
                fi->readDataPt("revData", i - 1, numTrials, masterSys->numSpecies, 0, lastRevStatePt);
            }

            for (int i = 0; i < numBoundedSpeciesStates; i++) {             
                revDists[i]->update(lastRevStatePt, numTrials);
            }
            
            fwdSetupEnd[i + 1] = omp_get_wtime();
            
            #pragma omp parallel for default(shared) private(sys)
            for (int j = 0; j < numTrials; j++) {                
                fwdTrialStart[i + 1][j] = omp_get_wtime();
                
                int threadId = omp_get_thread_num(); 

                sys = new System(*masterSys);
                sys->seed(rng());            

                for (int k = 0; k < sys->numSpecies; k++) {
                    sys->species[k]->state = lastRevStatePt[j][k];
                }

                for (int k = 0; k < sys->numRxns; k++) {
                    sys->rxns[k]->updateProp(sys->volRatio);
                }
                sys->initFwd();

                varName = "fwdData";
                simFwd(sys, numTimePts, time, state[threadId], revDists, speciesDistKey);
                
                fwdWriteStart[i + 1][j] = omp_get_wtime();
                writeStateData(i, sys, fi, varName, j, state[threadId], numTimePts, lock);
                fwdTrialEnd[i + 1][j] = omp_get_wtime();
            }
            
            revSetupStart[i + 1] = omp_get_wtime();
            
            masterSys->rxnPq->minHeap = false;
            masterSys->time = time[numTimePts - 1];
            fi->readDataPt("fwdData", i, numTrials, masterSys->numSpecies, numTimePts - 1, lastFwdStatePt);

            for (int i = 0; i < numBoundedSpeciesStates; i++) {             
                fwdDists[i]->update(lastFwdStatePt, numTrials);
            }
            
            revSetupEnd[i + 1] = omp_get_wtime();
            
            #pragma omp parallel for default(shared) private(sys)
            for (int j = 0; j < numTrials; j++) { 
                revTrialStart[i + 1][j] = omp_get_wtime();                
                int threadId = omp_get_thread_num(); 

                sys = new System(*masterSys);
                sys->seed(rng());            

                for (int k = 0; k < sys->numSpecies; k++) {
                    sys->species[k]->state = lastFwdStatePt[j][k];
                }

                for (int k = 0; k < sys->numRxns; k++) {
                    sys->rxns[k]->updateProp(sys->volRatio);
                }
                sys->initRev();

                varName = "revData";
                simRev(sys, numTimePts, time, state[threadId], revDists, speciesDistKey);
                
                revWriteStart[i + 1][j] = omp_get_wtime();
                writeStateData(i, sys, fi, varName, j, state[threadId], numTimePts, lock);
                revTrialEnd[i + 1][j] = omp_get_wtime();
            }
        }
    }
    
    double progEnd = omp_get_wtime();
        
    double avFwdSetupTime = 0;
    double avFwdTrialTime = 0;
    double avFwdRunTime = 0;
    double avFwdWriteTime = 0;
    
    double avRevSetupTime = 0;
    double avRevTrialTime = 0;
    double avRevRunTime = 0;
    double avRevWriteTime = 0;
    
    for (int i = 0; i < numDataSavePts; i++) {
        avFwdSetupTime += fwdSetupEnd[i] - fwdSetupStart[i];
        avRevSetupTime += revSetupEnd[i] - revSetupStart[i];
            
        for (int j = 0; j < numTrials; j++) {
            avFwdTrialTime += fwdTrialEnd[i][j] - fwdTrialStart[i][j];
            avFwdRunTime += fwdWriteStart[i][j] - fwdTrialStart[i][j];
            avFwdWriteTime += fwdTrialEnd[i][j] - fwdWriteStart[i][j];
            
            avRevTrialTime += revTrialEnd[i][j] - revTrialStart[i][j];
            avRevRunTime += revWriteStart[i][j] - revTrialStart[i][j];
            avRevWriteTime += revTrialEnd[i][j] - revWriteStart[i][j];
        }
    }
    
    avFwdSetupTime /= numDataSavePts;
    avFwdTrialTime /= (numDataSavePts * numTrials);
    avFwdRunTime /= (numDataSavePts * numTrials);
    avFwdWriteTime /= (numDataSavePts * numTrials);
    
    avRevSetupTime /= numDataSavePts;
    avRevTrialTime /= (numDataSavePts * numTrials);
    avRevRunTime /= (numDataSavePts * numTrials);
    avRevWriteTime /= (numDataSavePts * numTrials);
        
    double avSetupTime = (avFwdSetupTime + avRevSetupTime) / 2;
    double avTrialTime = (avFwdTrialTime + avRevTrialTime) / 2;
    double avRunTime = (avFwdRunTime + avRevRunTime) / 2;
    double avWriteTime = (avFwdWriteTime + avRevWriteTime) / 2;
    
    fprintf(stdout, "\nForward simulations:\n");
    fprintf(stdout, "    Average setup time:           %e s\n", avFwdSetupTime);
    fprintf(stdout, "    Average simulation run time:  %e s\n", avFwdRunTime);
    fprintf(stdout, "    Average file write time:      %e s\n", avFwdWriteTime);
    fprintf(stdout, "    Average trial execution time: %e s\n\n", avFwdTrialTime);
    
    fprintf(stdout, "Reverse simulations:\n");
    fprintf(stdout, "    Average setup time:           %e s\n", avRevSetupTime);
    fprintf(stdout, "    Average simulation run time:  %e s\n", avRevRunTime);
    fprintf(stdout, "    Average file write time:      %e s\n", avRevWriteTime);
    fprintf(stdout, "    Average trial execution time: %e s\n\n", avRevTrialTime);
    
    fprintf(stdout, "Forward and reverse simulations:\n");
    fprintf(stdout, "    Average setup time:           %e s\n", avSetupTime);
    fprintf(stdout, "    Average simulation run time:  %e s\n", avRunTime);
    fprintf(stdout, "    Average file write time:      %e s\n", avWriteTime);
    fprintf(stdout, "    Average trial execution time: %e s\n\n", avTrialTime);
    
    fprintf(stdout, "Total program execution time:          %e s\n\n", progEnd - progStart);
        
    delete[] dataSavePts;
    delete[] time;
    
    for (int i = 0; i < numTrials; i++) {
        delete[] lastFwdStatePt[i];
        delete[] lastRevStatePt[i];
    }
    delete[] lastFwdStatePt;
    delete[] lastRevStatePt;
    
    for (int i = 0; i < numThreads; i++) {
        for (int j = 0; j < masterSys->numSpecies; j++) {
            delete[] state[i][j];
        }        
        delete[] state[i];
    }
    delete[] state;
    
    delete[] distSpeciesKey;
    delete[] speciesDistKey;
    
    for (int i = 0; i < numBoundedSpeciesStates; i++) {
        delete revDists[i];
        delete revDistsPrev[i];
        delete fwdDists[i];
        delete fwdDistsPrev[i];
    }
    delete[] revDists;
    delete[] revDistsPrev;
    delete[] fwdDists;
    delete[] fwdDistsPrev;
    
    delete masterSys;    
    delete fi;
    
    return 0;
}

void simFwd(System* sys, int numTimePts, double* time, double** state, Distribution** dists, int* speciesDistKey) {
    int boundBreachSpeciesId;
    
    for (int i = 0; i < numTimePts; i++) {
        while(sys->rxnPq->getNextTime() <= time[i]) {
            boundBreachSpeciesId = sys->execRxn(true);
            
            if (boundBreachSpeciesId >= 0) {
                i = 0;
                sys->updateTime(time[i]);
                sys->species[boundBreachSpeciesId]->state = dists[speciesDistKey[boundBreachSpeciesId]]->sample();
                sys->initFwd();
            }
        }
        
        sys->updateTime(time[i]);

        for (int j = 0; j < sys->numSpecies; j++) {
            state[j][i] = sys->species[j]->state;
        }
    }
}

void simRev(System* sys, int numTimePts, double* time, double** state, Distribution** dists, int* speciesDistKey) {
    int boundBreachSpeciesId;
    
    for (int i = numTimePts - 1; i >= 0; i--) {
        double nextTime = sys->rxnPq->getNextTime();
        
        while(nextTime >= time[i] && nextTime != -DBL_MAX) {
            boundBreachSpeciesId = sys->execRxn(false);
            nextTime = sys->rxnPq->getNextTime();
            
            if (boundBreachSpeciesId >= 0) {
                i = numTimePts - 1;
                sys->updateTime(time[i]);
                sys->species[boundBreachSpeciesId]->state = dists[speciesDistKey[boundBreachSpeciesId]]->sample();
                sys->initRev();
            }
        }
        
        sys->updateTime(time[i]);
        
        for (int j = 0; j < sys->numSpecies; j++) {
            state[j][i] = sys->species[j]->state;
        }
    }
}

void writeStateData(int dataSavePtId, System* sys, FileInterface* fi, string varName, int trial, double** state, int numTimePts, omp_lock_t& lock) {
    omp_set_lock(&lock);
    if (dataSavePtId < 0) { 
        fi->writeInitStateData(varName, trial, state, sys->numSpecies, numTimePts);
    } else {        
        fi->writeStateData(varName, dataSavePtId, trial, state, sys->numSpecies, numTimePts);
    }
    omp_unset_lock(&lock);
}

/*
void averageStateData(int numTrials, int numTimePts, int numSpecies, double*** stateData, double** avStateData, double** lastAvStateData) {
    for (int i = 0; i < numTimePts; i++) {
        for (int j = 0; j < numSpecies; j++) {
            lastAvStateData[i][j] = avStateData[i][j];
            avStateData[i][j] = 0;
        }
    }
    
    for (int i = 0; i < numTrials; i++) {
        for (int j = 0; j < numTimePts; j++) {
            for (int k = 0; k < numSpecies; k++) {
                avStateData[j][k] += stateData[i][j][k];
            }
        }
    }
    
    for (int i = 0; i < numTimePts; i++) {
        for (int j = 0; j < numSpecies; j) {
            avStateData[i][j] /= numTrials;
        }
    }
}*/