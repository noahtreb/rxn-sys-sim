#include <assert.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "Distribution.h"
#include "FileInterface.h"
#include "PriorityQueue.h"
#include "Reaction.h"
#include "Species.h"
#include "System.h"

#include "main.h"

using namespace std;

RerunTrial::RerunTrial(int id, int trialId, int startTimePt, int endTimePt, int boundBreachSpeciesId) {
    this->id = id;
    this->trialId = trialId;
    this->startTimePt = startTimePt;
    this->endTimePt = endTimePt;
    this->boundBreachSpeciesId = boundBreachSpeciesId;
    
    this->numPrevTrials = -1;
    this->threadAssigned = false;
    this->threadId = -1;
}

void rerunTrialFwd(FileInterface* fi, System* sys, Distribution** revDists, Distribution* dist, vector<RerunTrial*>* rerunTrials,
        mt19937* rng, mutex* vecLock, mutex* fileLock, mutex* absCurrLock, string varName, double** boundStatePts, double** state,
        double* time, int** absCurr, int* speciesDistKey, int boundHandlingMethod, int threadId, int dataSavePtId, int numTrials, 
        int numTimePts) {
    while (true) {
        vecLock->lock();
        
        int startTimePt;
        RerunTrial* rt;
        bool rtAssigned = false;
        
        for (int i = 0; i < numTimePts; i++) {
            for (int j = 0; j < rerunTrials[i].size(); j++) {
                if (!rerunTrials[i][j]->threadAssigned) {
                    startTimePt = i;
                    rt = rerunTrials[i][j];
                    rt->threadAssigned = true;
                    rt->threadId = threadId;
                    rtAssigned = true;
                    break;
                }                
            }
            
            if (rtAssigned) {
                break;
            }
        }
                
        vecLock->unlock();
        
        if (!rtAssigned) {
            return;
        }
        
        while (rt->numPrevTrials > 0) {
            this_thread::sleep_for(chrono::seconds(1));
        }
        
        fileLock->lock();
        if (dataSavePtId < 0) {
            fi->readInitDataPt(varName, numTrials, sys->numSpecies, startTimePt, boundStatePts);  
        } else {
            fi->readDataPt(varName, dataSavePtId, numTrials, sys->numSpecies, startTimePt, boundStatePts);            
        }
        fileLock->unlock();
        
        dist->setSpecies(sys->species[rt->boundBreachSpeciesId]);
        dist->update(boundStatePts, numTrials, true);
              
        sys->updateTime(time[startTimePt]);
        sys->species[rt->boundBreachSpeciesId]->state = dist->sample();
        sys->updateProps();
        sys->seed((*rng)());
        sys->initFwd();
        
        int origStartTimePt = startTimePt;
        
        while (true) {            
            int endTimePt = numTimePts;
            
            for (int i = startTimePt + 1; i < numTimePts; i++) {
                if (rerunTrials[i].size() > 0) {
                    endTimePt = i;
                    break;
                }
            }
            
            int boundBreachSpeciesId;
            bool boundBreached = simFwd(sys, startTimePt, endTimePt, time, state, revDists, speciesDistKey, boundHandlingMethod, boundBreachSpeciesId, endTimePt);
            
            writeStateData(dataSavePtId, sys, fi, varName, rt->trialId, state, startTimePt, endTimePt - startTimePt, fileLock);
               
            int stopTimePt;
            if (endTimePt == numTimePts) {
                stopTimePt = endTimePt - 1;
            } else {
                stopTimePt = endTimePt;
            }
            
            vecLock->lock();               
            
            for (int i = startTimePt + 1; i <= stopTimePt; i++) {                
                for (int j = 0; j < rerunTrials[i].size(); j++) {
                    rerunTrials[i][j]->numPrevTrials--;
                }
            } 
            
            if (endTimePt == numTimePts) {
                int rtId = rt->id;
                for (int i = rtId + 1; i < rerunTrials[origStartTimePt].size(); i++) {
                    rerunTrials[origStartTimePt][i]->id--;
                }
                
                delete rt;
                rerunTrials[origStartTimePt].erase(rerunTrials[origStartTimePt].begin() + rtId);
                vecLock->unlock();
                
                break;
            }               
            
            if (boundBreached) {
                absCurrLock->lock();
                absCurr[boundBreachSpeciesId][endTimePt]++;
                absCurrLock->unlock();                
                
                int rtId = rt->id;
                for (int i = rtId + 1; i < rerunTrials[origStartTimePt].size(); i++) {
                    rerunTrials[origStartTimePt][i]->id--;
                }
                
                rerunTrials[origStartTimePt].erase(rerunTrials[origStartTimePt].begin() + rtId);
                rerunTrials[endTimePt].push_back(rt);

                rt->threadAssigned = false;
                rt->threadId = -1;
                rt->numPrevTrials = 0;
                rt->id = rerunTrials[origStartTimePt].size() - 1;
                
                for (int i = 0; i < endTimePt; i++) {
                    for (int j = 0; j < rerunTrials[i].size(); j++) {
                        rt->numPrevTrials++;
                    }
                }
                vecLock->unlock();
                
                break;
            }
                        
            startTimePt = endTimePt;
            vecLock->unlock();
        }
    }
    
    delete sys;
}

void rerunTrialRev(FileInterface* fi, System* sys, Distribution** fwdDists, Distribution* dist, vector<RerunTrial*>* rerunTrials,
        mt19937* rng, mutex* vecLock, mutex* fileLock, mutex* absCurrLock, string varName, double** boundStatePts, double** state,
        double* time, int** absCurr, int* speciesDistKey, int boundHandlingMethod, int threadId, int dataSavePtId, int numTrials, 
        int numTimePts) {
    while (true) {
        vecLock->lock();
        
        int endTimePt;
        RerunTrial* rt;
        bool rtAssigned = false;
        
        for (int i = numTimePts - 1; i >= 0; i--) {
            for (int j = 0; j < rerunTrials[i].size(); j++) {
                if (!rerunTrials[i][j]->threadAssigned) {
                    endTimePt = i;
                    rt = rerunTrials[i][j];
                    rt->threadAssigned = true;
                    rt->threadId = threadId;
                    rtAssigned = true;
                    break;
                }                
            }
            
            if (rtAssigned) {
                break;
            }
        }
                
        vecLock->unlock();
        
        if (!rtAssigned) {
            return;
        }
        
        while (rt->numPrevTrials > 0) {
            this_thread::sleep_for(chrono::seconds(1));
        }
        
        fileLock->lock();
        if (dataSavePtId < 0) {
            fi->readInitDataPt(varName, numTrials, sys->numSpecies, endTimePt, boundStatePts);  
        } else {
            fi->readDataPt(varName, dataSavePtId, numTrials, sys->numSpecies, endTimePt, boundStatePts);            
        }           
        fileLock->unlock();
        
        dist->setSpecies(sys->species[rt->boundBreachSpeciesId]);
        dist->update(boundStatePts, numTrials, true);
              
        sys->updateTime(time[endTimePt]);
        sys->species[rt->boundBreachSpeciesId]->state = dist->sample();
        sys->updateProps();
        sys->seed((*rng)());
        sys->initRev();
        
        int origEndTimePt = endTimePt;        
        
        while (true) {            
            int startTimePt = 0;
            
            for (int i = endTimePt - 1; i >= 0; i--) {
                if (rerunTrials[i].size() > 0) {
                    startTimePt = i;
                    break;
                }
            }
            
            int boundBreachSpeciesId;
            bool boundBreached = simRev(sys, startTimePt, endTimePt, time, state, fwdDists, speciesDistKey, boundHandlingMethod, boundBreachSpeciesId, startTimePt);
            
            writeStateData(dataSavePtId, sys, fi, varName, rt->trialId, state, startTimePt, endTimePt - startTimePt, fileLock);
                        
            vecLock->lock();               
            
            for (int i = startTimePt; i < endTimePt; i++) {                
                for (int j = 0; j < rerunTrials[i].size(); j++) {
                    rerunTrials[i][j]->numPrevTrials--;
                }
            } 
            
            if (startTimePt == 0) {
                int rtId = rt->id;
                for (int i = rtId + 1; i < rerunTrials[origEndTimePt].size(); i++) {
                    rerunTrials[origEndTimePt][i]->id--;
                }
                
                delete rt;
                rerunTrials[origEndTimePt].erase(rerunTrials[origEndTimePt].begin() + rtId);
                vecLock->unlock();
                
                break;
            }               
            
            if (boundBreached) {
                absCurrLock->lock();
                absCurr[boundBreachSpeciesId][startTimePt]++;
                absCurrLock->unlock();
                
                int rtId = rt->id;
                for (int i = rtId + 1; i < rerunTrials[origEndTimePt].size(); i++) {
                    rerunTrials[origEndTimePt][i]->id--;
                }
                
                rerunTrials[origEndTimePt].erase(rerunTrials[origEndTimePt].begin() + rtId);
                rerunTrials[startTimePt].push_back(rt);

                rt->threadAssigned = false;
                rt->threadId = -1;
                rt->numPrevTrials = 0;
                rt->id = rerunTrials[startTimePt].size() - 1;
                
                for (int i = startTimePt + 1; i < numTimePts; i++) {
                    for (int j = 0; j < rerunTrials[i].size(); j++) {
                        rt->numPrevTrials++;
                    }
                }
                vecLock->unlock();
                
                break;
            }
                        
            endTimePt = startTimePt;
            vecLock->unlock();
        }
    }
    
    delete sys;
}

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
        fprintf(stderr, "                   [-e <stoppingTol>] [-b <boundMethod>]\n");
        fprintf(stderr, "            <fileName> Name of the netCDF file to read from and write to.\n");
        fprintf(stderr, "    -r          <seed> Integer used to seed the random number generator.\n");
        fprintf(stderr, "Optional arguments:\n");
        fprintf(stderr, "    -s <skipInitFwd> Use this flag to indicate that preexisting data from the 'initFwdData' variable\n");
        fprintf(stderr, "                     should be used in the forward refinement algorithm.\n");
        fprintf(stderr, "    -t <skipInitRev> Use this flag to indicate that preexisting data from the 'initRevData' variable\n");
        fprintf(stderr, "                     should be used in the reverse refinement algorithm.\n");
        fprintf(stderr, "    -u  <skipRefine> Use this flag to bypass the forward and reverse refinement algorithms.\n");
        fprintf(stderr, "    -e <stoppingTol> Maximum allowable L2 distance between successive probability distributions.\n");
        fprintf(stderr, "    -b <boundMethod> 1 for deletion. 2 for current.\n\n");        
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
    
    int** absCurr = new int*[masterSys->numSpecies];
    for (int i = 0; i < masterSys->numSpecies; i++) {
        absCurr[i] = new int[numTimePts];
        for (int j = 0; j < numTimePts; j++) {
            absCurr[i][j] = 0;
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
    
    omp_lock_t fileLock, boundLock, absCurrLock;
    omp_init_lock(&fileLock);
    omp_init_lock(&boundLock);
    omp_init_lock(&absCurrLock);
        
    System* sys;
    
    vector<RerunTrial*>* rerunTrials = new vector<RerunTrial*>[numTimePts];
    
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
    
    Distribution** threadDists = new Distribution*[numThreads];
    for (int i = 0; i < numThreads; i++) {
        threadDists[i] = new Distribution(masterSys->species[0], numTrials, rng());
    }
        
    thread* threads = new thread[numThreads];
       
    double progStart = omp_get_wtime();
    if (!skipInitFwd) {        
        fwdSetupStart[0] = omp_get_wtime();
        for (int i = 0; i < numBoundedSpeciesStates; i++) {
            revDists[i]->addNode(revDists[i]->species->state, numTrials);
        }
        fwdSetupEnd[0] = omp_get_wtime();
        
        int numRerunTrials = 0;
        #pragma omp parallel for default(shared) private(sys)
        for (int i = 0; i < numTrials; i++) {
            fwdTrialStart[0][i].push_back(omp_get_wtime());
            int threadId = omp_get_thread_num();            

            sys = new System(*masterSys);
            sys->seed(rng());
            sys->initFwd();

            int breachSpeciesId, breachTimePt;
            varName = "initFwdData";
            bool boundBreach = simFwd(sys, 0, numTimePts, time, state[threadId], revDists, speciesDistKey, boundHandlingMethod, breachSpeciesId, breachTimePt);
            
            fwdWriteStart[0][i].push_back(omp_get_wtime());
            writeStateData(-1, sys, fi, varName, i, state[threadId], 0, numTimePts, fileLock);
            fwdTrialEnd[0][i].push_back(omp_get_wtime());

            if (boundBreach) {
                omp_set_lock(&absCurrLock);
                absCurr[breachSpeciesId][breachTimePt]++;
                omp_unset_lock(&absCurrLock);
                
                omp_set_lock(&boundLock);
                int rtId = rerunTrials[breachTimePt].size();
                rerunTrials[breachTimePt].push_back(new RerunTrial(rtId, i, breachTimePt, numTimePts, breachSpeciesId));
                numRerunTrials++;
                omp_unset_lock(&boundLock);
            }
        }
              
        if (numRerunTrials > 0) {
            int cumCount = 0;
            for (int i = 0; i < numTimePts; i++) {
                for (int j = 0; j < rerunTrials[i].size(); j++) {
                    rerunTrials[i][j]->numPrevTrials = cumCount;
                }
                
                cumCount += rerunTrials[i].size();
            }
            
            mutex vecLock, fileLock, absCurrLock;
            
            for (int i = 0; i < numThreads; i++) {
                threads[i] = thread(rerunTrialFwd, fi, new System(*masterSys), revDists, threadDists[i], rerunTrials, &rng, &vecLock, &fileLock, &absCurrLock, varName, boundStatePts[i], state[i], time, absCurr, speciesDistKey, boundHandlingMethod, i, -1, numTrials, numTimePts);
            }
            
            for (int i = 0; i < numThreads; i++) {
                threads[i].join();
            }
        }
        
        fi->writeInitAbsCurrData("initFwdAbsCurr", absCurr, masterSys->numSpecies, 0, numTimePts);
        
        #pragma omp parallel for default(shared)
        for (int i = 0; i < masterSys->numSpecies; i++) {
            for (int j = 0; j < numTimePts; j++) {
                absCurr[i][j] = 0;
            }
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
        
        int numRerunTrials = 0;
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

            int breachSpeciesId, breachTimePt;
            varName = "initRevData";
            bool boundBreach = simRev(sys, 0, numTimePts, time, state[threadId], fwdDists, speciesDistKey, boundHandlingMethod, breachSpeciesId, breachTimePt);
            
            revWriteStart[0][i].push_back(omp_get_wtime());
            writeStateData(-1, sys, fi, varName, i, state[threadId], 0, numTimePts, fileLock);
            revTrialEnd[0][i].push_back(omp_get_wtime());
            
            if (boundBreach) {
                omp_set_lock(&absCurrLock);
                absCurr[breachSpeciesId][breachTimePt]++;
                omp_unset_lock(&absCurrLock);
                
                omp_set_lock(&boundLock);
                int rtId = rerunTrials[breachTimePt].size();
                rerunTrials[breachTimePt].push_back(new RerunTrial(rtId, i, 0, breachTimePt, breachSpeciesId));
                numRerunTrials++;
                omp_unset_lock(&boundLock);
            }
        }
        
        if (numRerunTrials > 0) {
            int cumCount = 0;
            for (int i = numTimePts - 1; i >= 0; i--) {
                for (int j = 0; j < rerunTrials[i].size(); j++) {
                    rerunTrials[i][j]->numPrevTrials = cumCount;
                }
                
                cumCount += rerunTrials[i].size();
            }
                        
            mutex vecLock, fileLock, absCurrLock;
            
            for (int i = 0; i < numThreads; i++) {
                threads[i] = thread(rerunTrialRev, fi, new System(*masterSys), fwdDists, threadDists[i], rerunTrials, &rng, &vecLock, &fileLock, &absCurrLock, varName, boundStatePts[i], state[i], time, absCurr, speciesDistKey, boundHandlingMethod, i, -1, numTrials, numTimePts);
            }
            
            for (int i = 0; i < numThreads; i++) {
                threads[i].join();
            }
        }
        
        fi->writeInitAbsCurrData("initRevAbsCurr", absCurr, masterSys->numSpecies, 0, numTimePts);
        
        #pragma omp parallel for default(shared)
        for (int i = 0; i < masterSys->numSpecies; i++) {
            for (int j = 0; j < numTimePts; j++) {
                absCurr[i][j] = 0;
            }
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

            for (int j = 0; j < numBoundedSpeciesStates; j++) {             
                revDists[j]->update(lastRevStatePt, numTrials, false);
            }
            
            fwdSetupEnd[i + 1] = omp_get_wtime();
            
            int numRerunTrials = 0;
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

                int breachSpeciesId, breachTimePt;
                varName = "fwdData";
                bool boundBreach = simFwd(sys, 0, numTimePts, time, state[threadId], revDists, speciesDistKey, boundHandlingMethod, breachSpeciesId, breachTimePt);
                
                fwdWriteStart[i + 1][j].push_back(omp_get_wtime());
                writeStateData(i, sys, fi, varName, j, state[threadId], 0, numTimePts, fileLock);
                fwdTrialEnd[i + 1][j].push_back(omp_get_wtime());
                
                if (boundBreach) {
                    omp_set_lock(&absCurrLock);
                    absCurr[breachSpeciesId][breachTimePt]++;
                    omp_unset_lock(&absCurrLock);

                    omp_set_lock(&boundLock);
                    int rtId = rerunTrials[breachTimePt].size();
                    rerunTrials[breachTimePt].push_back(new RerunTrial(rtId, j, 0, breachTimePt, breachSpeciesId));
                    numRerunTrials++;
                    omp_unset_lock(&boundLock);
                }
            }
            
            if (numRerunTrials > 0) {
                int cumCount = 0;
                for (int j = 0; j < numTimePts; j++) {
                    for (int k = 0; k < rerunTrials[j].size(); k++) {
                        rerunTrials[j][k]->numPrevTrials = cumCount;
                    }

                    cumCount += rerunTrials[j].size();
                }

                mutex vecLock, fileLock, absCurrLock;

                for (int j = 0; j < numThreads; j++) {
                    threads[j] = thread(rerunTrialFwd, fi, new System(*masterSys), revDists, threadDists[j], rerunTrials, &rng, &vecLock, &fileLock, &absCurrLock, varName, boundStatePts[j], state[j], time, absCurr, speciesDistKey, boundHandlingMethod, j, i, numTrials, numTimePts);
                }

                for (int j = 0; j < numThreads; j++) {
                    threads[j].join();
                }
            }

            fi->writeAbsCurrData("fwdAbsCurr", i, absCurr, masterSys->numSpecies, 0, numTimePts);

            #pragma omp parallel for default(shared)
            for (int j = 0; j < masterSys->numSpecies; j++) {
                for (int k = 0; k < numTimePts; k++) {
                    absCurr[j][k] = 0;
                }
            }
            
            revSetupStart[i + 1] = omp_get_wtime();
            
            masterSys->rxnPq->minHeap = false;
            masterSys->time = time[numTimePts - 1];
            fi->readDataPt("fwdData", i, numTrials, masterSys->numSpecies, numTimePts - 1, lastFwdStatePt);

            for (int i = 0; i < numBoundedSpeciesStates; i++) {             
                fwdDists[i]->update(lastFwdStatePt, numTrials, false);
            }
            
            revSetupEnd[i + 1] = omp_get_wtime();
            
            numRerunTrials = 0;
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
                
                int breachSpeciesId, breachTimePt;
                varName = "revData";
                bool boundBreach = simRev(sys, 0, numTimePts, time, state[threadId], revDists, speciesDistKey, boundHandlingMethod, breachSpeciesId, breachTimePt);
                
                revWriteStart[i + 1][j].push_back(omp_get_wtime());
                writeStateData(i, sys, fi, varName, j, state[threadId], 0, numTimePts, fileLock);
                revTrialEnd[i + 1][j].push_back(omp_get_wtime());
                
                if (boundBreach) {
                    omp_set_lock(&absCurrLock);
                    absCurr[breachSpeciesId][breachTimePt]++;
                    omp_unset_lock(&absCurrLock);

                    omp_set_lock(&boundLock);
                    int rtId = rerunTrials[breachTimePt].size();
                    rerunTrials[breachTimePt].push_back(new RerunTrial(rtId, j, 0, breachTimePt, breachSpeciesId));
                    numRerunTrials++;
                    omp_unset_lock(&boundLock);
                }
            }
            
            if (numRerunTrials > 0) {
                int cumCount = 0;
                for (int j = numTimePts - 1; j >= 0; j--) {
                    for (int k = 0; k < rerunTrials[j].size(); k++) {
                        rerunTrials[j][k]->numPrevTrials = cumCount;
                    }

                    cumCount += rerunTrials[j].size();
                }

                mutex vecLock, fileLock, absCurrLock;

                for (int j = 0; j < numThreads; j++) {
                    threads[j] = thread(rerunTrialRev, fi, new System(*masterSys), fwdDists, threadDists[j], rerunTrials, &rng, &vecLock, &fileLock, &absCurrLock, varName, boundStatePts[j], state[j], time, absCurr, speciesDistKey, boundHandlingMethod, j, i, numTrials, numTimePts);
                }

                for (int j = 0; j < numThreads; j++) {
                    threads[j].join();
                }
            }

            fi->writeAbsCurrData("revAbsCurr", i, absCurr, masterSys->numSpecies, 0, numTimePts);

            #pragma omp parallel for default(shared)
            for (int j = 0; j < masterSys->numSpecies; j++) {
                for (int k = 0; k < numTimePts; k++) {
                    absCurr[j][k] = 0;
                }
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

bool simFwd(System* sys, int startTimePt, int endTimePt, double* time, double** state, Distribution** dists, int* speciesDistKey, int boundHandlingMethod, int& boundBreachSpeciesId, int& breachTimePt) {
    bool boundStop = false;
    
    for (int i = startTimePt; i < endTimePt; i++) {
        while(sys->rxnPq->getNextTime() <= time[i] && !boundStop) {
            boundBreachSpeciesId = sys->execRxn(true);
            
            if (boundBreachSpeciesId >= 0) {
                switch (boundHandlingMethod) {
                    case 1:                        
                        i = 0;
                        sys->updateTime(time[i]);
                        sys->species[boundBreachSpeciesId]->state = dists[speciesDistKey[boundBreachSpeciesId]]->sample();
                        sys->updateProps();
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

bool simRev(System* sys, int startTimePt, int endTimePt, double* time, double** state, Distribution** dists, int* speciesDistKey, int boundHandlingMethod, int& boundBreachSpeciesId, int& breachTimePt) {
    bool boundStop = false;
    
    int numTimePts = endTimePt;
    for (int i = numTimePts - 1; i >= startTimePt; i--) {
        double nextTime = sys->rxnPq->getNextTime();
        
        while(nextTime >= time[i] && nextTime != -DBL_MAX && !boundStop) {
            boundBreachSpeciesId = sys->execRxn(false);
            nextTime = sys->rxnPq->getNextTime();
            
            if (boundBreachSpeciesId >= 0) {
                switch(boundHandlingMethod) {
                    case 1:
                        i = numTimePts - 1;
                        sys->updateTime(time[i]);
                        sys->species[boundBreachSpeciesId]->state = dists[speciesDistKey[boundBreachSpeciesId]]->sample();
                        sys->updateProps();
                        sys->initRev();
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
        for (int i = startTimePt; i < breachTimePt; i++) {
            state[boundBreachSpeciesId][i] = -1;
        }
    }
    
    return boundStop;
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

void writeStateData(int dataSavePtId, System* sys, FileInterface* fi, string varName, int trial, double** state, int startTimePt, int numTimePts, mutex* fileLock) {
    fileLock->lock();
    if (dataSavePtId < 0) { 
        fi->writeInitStateData(varName, trial, state, sys->numSpecies, startTimePt, numTimePts);
    } else {        
        fi->writeStateData(varName, dataSavePtId, trial, state, sys->numSpecies, startTimePt, numTimePts);
    }
    fileLock->unlock();
}

bool rerunTrialSort(RerunTrial left, RerunTrial right) {
    return left.endTimePt < right.endTimePt;
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