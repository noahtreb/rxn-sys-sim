#include <assert.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include "main.h"

#include "Distribution.h"
#include "FileInterface.h"
#include "PriorityQueue.h"
#include "Reaction.h"
#include "Species.h"
#include "System.h"

using namespace std;

int main(int argc, const char* argv[]) {
    bool skipInitFwd = false;
    bool skipInitRev = false;
    bool skipRefine = false;
    int seed;
    double stoppingTol = -1;
    string fileName;
    int boundHandlingMethod = 0;
    
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
                    case 'b':
                        boundHandlingMethod = atoi(argv[i+1]);
                        i++;
                        break;
                }
            } else {
                fileName = argv[i];
            }
        }
    } else {
        fprintf(stderr, "Usage: ./RxnSysSim <fileName> -r <seed> [-s <skipInitFwd>] [-t <skipInitRev>] [-u <skipRefine>]\n");
        fprintf(stderr, "                   [-e <stoppingTol>] [-b <boundHandlingMethod>]\n");
        fprintf(stderr, "            <fileName> Name of the netCDF file to read from and write to.\n");
        fprintf(stderr, "    -r          <seed> Integer used to seed the random number generator.\n");
        fprintf(stderr, "Optional arguments:\n");
        fprintf(stderr, "    -s         <skipInitFwd> Use this flag to indicate that preexisting data from the 'initFwdData' variable\n");
        fprintf(stderr, "                             should be used in the forward refinement algorithm.\n");
        fprintf(stderr, "    -t         <skipInitRev> Use this flag to indicate that preexisting data from the 'initRevData' variable\n");
        fprintf(stderr, "                             should be used in the reverse refinement algorithm.\n");
        fprintf(stderr, "    -u          <skipRefine> Use this flag to bypass the forward and reverse refinement algorithms.\n");
        fprintf(stderr, "    -e         <stoppingTol> Maximum allowable L2 distance between successive probability distributions.\n");
        fprintf(stderr, "    -b <boundHandlingMethod> 1 for deletion. 2 for current.\n\n");        
        abort();
    }
    
    FileInterface* fi = new FileInterface(fileName);
    std::mt19937 rng(seed);   
        
    if (stoppingTol > 0) {
        fi->overwriteStoppingTol(stoppingTol);
    }
    
    if (boundHandlingMethod > 0) {
        fi->overwriteBoundHandlingMethod(boundHandlingMethod);
    }
    
    int numTrials;
    double startTime;
    double endTime;
    int numTimePts;
    double timeStep;
    int numDataSavePts;
    int* dataSavePts;
    int numBoundedSpeciesStates;
    
    System* masterSys = fi->readFileData(numTrials, startTime, endTime, numTimePts, timeStep, stoppingTol, numDataSavePts, dataSavePts, boundHandlingMethod, numBoundedSpeciesStates);
    
    if (skipRefine) {
        numDataSavePts = 0;
    }
    
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
    
    double*** boundStatePts = new double**[numThreads];
    for (int i = 0; i < numThreads; i++) {
        boundStatePts[i] = new double*[numTrials];        
        
        for (int j = 0; j < numTrials; j++) {
            boundStatePts[i][j] = new double[masterSys->numSpecies];
        }
    }
    
    omp_lock_t fileLock, boundLock;
    omp_init_lock(&fileLock);
    omp_init_lock(&boundLock);
        
    System* sys;
    
    int boundBreachSpeciesId;
    
    vector< pair<int, int> >* runTrials = new vector< pair<int, int> >();
    vector< pair<int, int> >* rerunTrials = new vector< pair<int, int> >();
    
    if (skipRefine) {
        numDataSavePts = 0;
    }
    
    double fwdSetupStart[numDataSavePts + 1];
    double fwdSetupEnd[numDataSavePts + 1];
    vector<double> fwdTrialStart[numDataSavePts + 1][numTrials];
    vector<double> fwdWriteStart[numDataSavePts + 1][numTrials];
    vector<double> fwdTrialEnd[numDataSavePts + 1][numTrials];
    double revSetupStart[numDataSavePts + 1];
    double revSetupEnd[numDataSavePts + 1];
    vector<double> revTrialStart[numDataSavePts + 1][numTrials];
    vector<double> revWriteStart[numDataSavePts + 1][numTrials];
    vector<double> revTrialEnd[numDataSavePts + 1][numTrials];
    
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
            fwdTrialStart[0][i].push_back(omp_get_wtime());
            int threadId = omp_get_thread_num();            

            sys = new System(*masterSys);
            sys->seed(rng());
            sys->initFwd();

            int speciesId, endTimePt;
            varName = "initFwdData";
            bool boundBreach = simFwd(sys, 0, numTimePts, time, state[threadId], revDists, speciesDistKey, boundHandlingMethod, speciesId, endTimePt);
            
            fwdWriteStart[0][i].push_back(omp_get_wtime());
            writeStateData(-1, sys, fi, varName, i, state[threadId], 0, numTimePts, fileLock);
            fwdTrialEnd[0][i].push_back(omp_get_wtime());

            if (boundBreach) {
                omp_set_lock(&boundLock);
                runTrials->push_back(make_pair(i, endTimePt));
                boundBreachSpeciesId = speciesId;
                omp_unset_lock(&boundLock);
            }
        }
        
        int numRerunTrials = runTrials->size();
        int numPrevTrials[numRerunTrials];
        int numExecPrevTrials[numRerunTrials];
        
        Distribution** threadDists = new Distribution*[numThreads];
        for (int i = 0; i < numThreads; i++) {
            threadDists[i] = new Distribution(masterSys->species[boundBreachSpeciesId], numTrials, rng());
        }
        
        while (numRerunTrials > 0) {            
            std::sort(runTrials->begin(), runTrials->end(), pairSort);
            
            int lastCount = 0;
            int lastTimePt = runTrials->at(0).second;
            numPrevTrials[0] = 0;
            numExecPrevTrials[0] = 0;
            
            for (int i = 1; i < numRerunTrials; i++) {                
                if (runTrials->at(i).second != lastTimePt) {
                    lastCount = i;
                    lastTimePt = runTrials->at(i).second;
                }
                
                numPrevTrials[i] = lastCount;
                numExecPrevTrials[i] = 0;
            }
                        
            #pragma omp parallel for default(shared) private(sys) schedule(static, 1)
            for (int i = 0; i < numRerunTrials; i++) {
                int threadId = omp_get_thread_num();
                int trialId = runTrials->at(i).first;
                int startTimePt = runTrials->at(i).second;
                
                fwdTrialStart[0][trialId].push_back(omp_get_wtime());
                                                
                while(numExecPrevTrials[i] != numPrevTrials[i]) {
                    usleep(1000);
                }
                
                omp_set_lock(&fileLock);
                fi->readInitDataPt(varName, numTrials, masterSys->numSpecies, startTimePt, boundStatePts[threadId]);            
                omp_unset_lock(&fileLock);
                
                threadDists[threadId]->update(boundStatePts[threadId], numTrials, true);
                
                sys = new System(*masterSys);
                sys->seed(rng());                
                sys->updateTime(time[startTimePt]);
                sys->species[boundBreachSpeciesId]->state = threadDists[threadId]->sample();
                sys->initFwd();

                int nextTimePt = 0;
                int nextTimePtIndex = i+1;
                
                nextTimePt = numTimePts; 
                for (int j = nextTimePtIndex; j < numRerunTrials; j++) {
                    if (runTrials->at(j).second != startTimePt) {
                        nextTimePt = runTrials->at(j).second;
                        nextTimePtIndex = j;
                        break;
                    } else {
                        nextTimePtIndex = numRerunTrials;
                    }
                }
                
                while (true) {
                    int endTimePt;
                    varName = "initFwdData";
                    bool boundBreached = simFwd(sys, startTimePt, nextTimePt, time, state[threadId], revDists, speciesDistKey, boundHandlingMethod, boundBreachSpeciesId, endTimePt);

                    fwdWriteStart[0][trialId].push_back(omp_get_wtime());
                    writeStateData(-1, sys, fi, varName, trialId, state[threadId], startTimePt, nextTimePt - startTimePt, fileLock);
                    fwdTrialEnd[0][trialId].push_back(omp_get_wtime());

                    if (boundBreached) {
                        omp_set_lock(&boundLock);
                        rerunTrials->push_back(make_pair(trialId, endTimePt));
                        omp_unset_lock(&boundLock);
                        
                        break;
                    }
                    
                    omp_set_lock(&boundLock);
                    int j = 0;
                    for (j = nextTimePtIndex; j < numRerunTrials; j++) {
                        if (runTrials->at(j).second == nextTimePt) {
                            numExecPrevTrials[j]++;
                        } else {
                            nextTimePtIndex = j;
                            nextTimePt = runTrials->at(j).second;
                            break;
                        }
                    }
                    omp_unset_lock(&boundLock);
                    
                    if (j == numRerunTrials) {
                        break;
                    }
                }
            }
            
            vector< pair<int, int> >* temp = runTrials;
            runTrials = rerunTrials;
            rerunTrials = temp;
            rerunTrials->clear();
            
            numRerunTrials = runTrials->size();
        }
    } else {
        fwdSetupStart[0] = 0;
        fwdSetupEnd[0] = 0;
        
        for (int i = 0; i < numTrials; i++) {
            fwdTrialStart[0][i].push_back(0);
            fwdWriteStart[0][i].push_back(0);
            fwdTrialEnd[0][i].push_back(0);
        }
    }
    
    if (!skipInitRev) {
        revSetupStart[0] = omp_get_wtime();
        
        masterSys->rxnPq->minHeap = false;
        masterSys->time = time[numTimePts - 1];
        fi->readInitDataPt("initFwdData", numTrials, masterSys->numSpecies, numTimePts - 1, lastFwdStatePt);
        
        for (int i = 0; i < numBoundedSpeciesStates; i++) {             
            fwdDists[i]->update(lastFwdStatePt, numTrials, false);
        }
                
        revSetupEnd[0] = omp_get_wtime();
        
        #pragma omp parallel for default(shared) private(sys)
        for (int i = 0; i < numTrials; i++) {    
            revTrialStart[0][i].push_back(omp_get_wtime());
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
            simRev(sys, 0, numTimePts, time, state[threadId], fwdDists, speciesDistKey, boundHandlingMethod);
            
            revWriteStart[0][i].push_back(omp_get_wtime());
            writeStateData(-1, sys, fi, varName, i, state[threadId], 0, numTimePts, fileLock);
            revTrialEnd[0][i].push_back(omp_get_wtime());
        }
    } else {
        revSetupStart[0] = 0;
        revSetupEnd[0] = 0;
        
        for (int i = 0; i < numTrials; i++) {
            revTrialStart[0][i].push_back(0);
            revWriteStart[0][i].push_back(0);
            revTrialEnd[0][i].push_back(0);
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
                revDists[i]->update(lastRevStatePt, numTrials, false);
            }
            
            fwdSetupEnd[i + 1] = omp_get_wtime();
            
            #pragma omp parallel for default(shared) private(sys)
            for (int j = 0; j < numTrials; j++) {                
                fwdTrialStart[i + 1][j].push_back(omp_get_wtime());
                
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

                int temp1, temp2;
                varName = "fwdData";
                simFwd(sys, 0, numTimePts, time, state[threadId], revDists, speciesDistKey, boundHandlingMethod, temp1, temp2);
                
                fwdWriteStart[i + 1][j].push_back(omp_get_wtime());
                writeStateData(i, sys, fi, varName, j, state[threadId], 0, numTimePts, fileLock);
                fwdTrialEnd[i + 1][j].push_back(omp_get_wtime());
            }
            
            revSetupStart[i + 1] = omp_get_wtime();
            
            masterSys->rxnPq->minHeap = false;
            masterSys->time = time[numTimePts - 1];
            fi->readDataPt("fwdData", i, numTrials, masterSys->numSpecies, numTimePts - 1, lastFwdStatePt);

            for (int i = 0; i < numBoundedSpeciesStates; i++) {             
                fwdDists[i]->update(lastFwdStatePt, numTrials, false);
            }
            
            revSetupEnd[i + 1] = omp_get_wtime();
            
            #pragma omp parallel for default(shared) private(sys)
            for (int j = 0; j < numTrials; j++) { 
                revTrialStart[i + 1][j].push_back(omp_get_wtime());                
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
                simRev(sys, 0, numTimePts, time, state[threadId], revDists, speciesDistKey, boundHandlingMethod);
                
                revWriteStart[i + 1][j].push_back(omp_get_wtime());
                writeStateData(i, sys, fi, varName, j, state[threadId], 0, numTimePts, fileLock);
                revTrialEnd[i + 1][j].push_back(omp_get_wtime());
            }
        }
    }
    
    double progEnd = omp_get_wtime();
        
    numDataSavePts++;
    
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
            for (int k = 0; k < fwdTrialEnd[i][j].size(); k++) {
                avFwdTrialTime += fwdTrialEnd[i][j][k] - fwdTrialStart[i][j][k];
                avFwdRunTime += fwdWriteStart[i][j][k] - fwdTrialStart[i][j][k];
                avFwdWriteTime += fwdTrialEnd[i][j][k] - fwdWriteStart[i][j][k];
            }
            
            for (int k = 0; k < revTrialEnd[i][j].size(); k++) {
                avRevTrialTime += revTrialEnd[i][j][k] - revTrialStart[i][j][k];
                avRevRunTime += revWriteStart[i][j][k] - revTrialStart[i][j][k];
                avRevWriteTime += revTrialEnd[i][j][k] - revWriteStart[i][j][k];
            }
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

int simFwd(System* sys, int startTimePt, int endTimePt, double* time, double** state, Distribution** dists, int* speciesDistKey, int boundHandlingMethod, int& boundBreachSpeciesId, int& breachTimePt) {
    bool boundStop = false;
    
    for (int i = startTimePt; i < endTimePt; i++) {
        fprintf(stderr, "%i", i);
        while(sys->rxnPq->getNextTime() <= time[i] && !boundStop) {
            boundBreachSpeciesId = sys->execRxn(true);
            
            if (boundBreachSpeciesId >= 0) {
                switch (boundHandlingMethod) {
                    case 1:                        
                        i = 0;
                        sys->updateTime(time[i]);
                        sys->species[boundBreachSpeciesId]->state = dists[speciesDistKey[boundBreachSpeciesId]]->sample();
                        sys->initFwd();
                        break;
                    case 2:
                        breachTimePt = i;
                        boundStop = true;
                        break;
                }
            }
        }
        
        sys->updateTime(time[i]);

        for (int j = 0; j < sys->numSpecies; j++) {
            state[j][i] = sys->species[j]->state;
        }
        
        if (boundStop) {
            break;
        }        
    }
    
    if (boundStop) {
        for (int i = breachTimePt+1; i < endTimePt; i++) {
            state[boundBreachSpeciesId][i] = -1;
        }
    }
    
    return boundStop;
}

int simRev(System* sys, int startTimePt, int endTimePt, double* time, double** state, Distribution** dists, int* speciesDistKey, int boundHandlingMethod) {
    int boundBreachSpeciesId;
    
    int numTimePts = endTimePt;
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

void writeStateData(int dataSavePtId, System* sys, FileInterface* fi, string varName, int trial, double** state, int startTimePt, int numTimePts, omp_lock_t& fileLock) {
    omp_set_lock(&fileLock);
    if (dataSavePtId < 0) { 
        fi->writeInitStateData(varName, trial, state, sys->numSpecies, startTimePt, numTimePts);
    } else {        
        fi->writeStateData(varName, dataSavePtId, trial, state, sys->numSpecies, startTimePt, numTimePts);
    }
    omp_unset_lock(&fileLock);
}

bool pairSort(pair<int, int> left, pair<int, int> right) {
    return left.second < right.second;
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