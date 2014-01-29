#include "Reaction.h"
#include "Species.h"
#include <stdio.h>

using namespace std;

Reaction::Reaction(const int id,
        Species** species,
        const int numStoichSpecies, 
        int* stoichSpeciesIds, 
        const int* const stoichCoeffs, 
        const int rateLaw, 
        const int numRateConsts, 
        double* rateConsts, 
        const int numRateSpecies, 
        int* rateSpeciesIds, 
        const int numDeps,
        const int* const deps,
        double vol) : 
        id(id),
        numStoichSpecies(numStoichSpecies),
        stoichCoeffs(stoichCoeffs), 
        rateLaw(rateLaw), 
        numRateConsts(numRateConsts),
        rateConsts(rateConsts),
        numRateSpecies(numRateSpecies),
        numDeps(numDeps),
        deps(deps) {
    this->prop = 0;   
    this->oldProp = 0;
    
    this->stoichSpecies = new Species*[numStoichSpecies];    
    for (int i = 0; i < numStoichSpecies; i++) {
        this->stoichSpecies[i] = species[stoichSpeciesIds[i]];
    }
    
    this->rateSpecies = new Species*[numRateSpecies];
    for (int i = 0; i < numRateSpecies; i++) {
        this->rateSpecies[i] = species[rateSpeciesIds[i]];
    }
    
    switch(rateLaw) { // Convert to mesoscopic units.
        case 3: // Second-order (bimolecular) mass action            
            //this->rateConsts[0] = this->rateConsts[0] / vol / 6.022e23;
            break;
    }
}

Reaction::Reaction(const Reaction& other) : 
        id(other.id),
        numStoichSpecies(other.numStoichSpecies), 
        stoichCoeffs(other.stoichCoeffs), 
        rateLaw(other.rateLaw), 
        numRateConsts(other.numRateConsts),
        rateConsts(other.rateConsts),
        numRateSpecies(other.numRateSpecies),
        numDeps(other.numDeps),
        deps(other.deps) {
    this->prop = other.prop;
    this->oldProp = other.oldProp;
    
    this->stoichSpecies = new Species*[this->numStoichSpecies];    
    for (int i = 0; i < this->numStoichSpecies; i++) {
        this->stoichSpecies[i] = other.stoichSpecies[i];
    }
    
    this->rateSpecies = new Species*[this->numRateSpecies];
    for (int i = 0; i < this->numRateSpecies; i++) {
        this->rateSpecies[i] = other.rateSpecies[i];
    }
}

Reaction::~Reaction() {
    for (int i = 0; i < this->numStoichSpecies; i++) {
        this->stoichSpecies[i] = NULL;
    }
    
    for (int i = 0; i < this->numRateSpecies; i++) {
        this->rateSpecies[i] = NULL;
    }
    
    delete[] this->stoichSpecies;
    delete[] this->stoichCoeffs;
    delete[] this->rateConsts;
    delete[] this->rateSpecies;
}

void Reaction::updateProp(double volRatio) {
    this->oldProp = this->prop;    
    this->prop = this->rateConsts[0];
    
    for (int i = 0; i < this->numStoichSpecies; i++) {
        if (this->stoichCoeffs[i] < 0) {
            for (int j = 0; j < -this->stoichCoeffs[i]; j++) {
                this->prop *= 1.0 * (this->stoichSpecies[i]->state - j);// / (j + 1);
            }
        }
    }
    
    /*switch(this->rateLaw) {
        case 2: // First-order mass action
            this->prop = this->rateConsts[0] * speciesState[this->rateSpeciesIds[0]];
            break;
        case 3: // Second-order (bimolecular) mass action
            this->prop = this->rateConsts[0] * speciesState[this->rateSpeciesIds[0]] * speciesState[this->rateSpeciesIds[1]] / volRatio;
            break;
    }*/
}

void Reaction::print() const {
    int stoichCoeff;
    int stoichSpeciesId;
    
    for (int i = 0; i < this->numStoichSpecies; i++) {
        stoichCoeff = this->stoichCoeffs[i];
        stoichSpeciesId = this->stoichSpecies[i]->id;
        
        if (i == 0 && stoichCoeff > 0) {            
            fprintf(stdout, " -> ");
        }
        
        if (stoichCoeff < 0) {
            fprintf(stdout, "%i(%i)", -stoichCoeff, stoichSpeciesId);
        } else {            
            fprintf(stdout, "%i(%i)", stoichCoeff, stoichSpeciesId);
        }
        
        if (i != this->numStoichSpecies - 1) {
            if (stoichCoeff * this->stoichCoeffs[i + 1] > 0) {
                fprintf(stdout, " + ");
            } else {
                fprintf(stdout, " -> ");
            }
        }
        
        if (i == this->numStoichSpecies - 1 && stoichCoeff < 0) {
            fprintf(stdout, " -> ");
        }
    }
    
    fprintf(stdout, "\n    ID: %i, Rate Constant: %e, Propensity: %e, Old Propensity: %e\n", this->id, this->rateConsts[0], this->prop, this->oldProp);
}