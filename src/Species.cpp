#include "Species.h"
#include <stdlib.h>

using namespace std;

Species::Species(const int id,
        double state, 
        const bool stateChanges,
        const double stateUpperBound, 
        const double stateLowerBound, 
        const bool stateBounded) : 
        id(id),
        stateChanges(stateChanges),
        stateUpperBound(stateUpperBound),
        stateLowerBound(stateLowerBound),
        stateBounded(stateBounded) {
    this->state = state;
}

Species::Species(const Species& other) :
        id(other.id),
        stateChanges(other.stateChanges),
        stateUpperBound(other.stateUpperBound),
        stateLowerBound(other.stateLowerBound),
        stateBounded(other.stateBounded) {
    this->state = other.state;
}

Species::~Species() {
    
}