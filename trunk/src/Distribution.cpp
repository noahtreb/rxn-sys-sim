#include "Distribution.h"
#include "Species.h"
#include <algorithm>
#include <math.h>
#include <random>
#include <vector>

using namespace std;

DistNode::DistNode(double state, int count, int cumCount) {
    this->state = state;
    this->count = count;
    this->cumCount = cumCount;
}

DistNode::~DistNode() {
    
}

void DistNode::print() const {
    fprintf(stdout, "%e, %i, %i\n", this->state, this->count, this->cumCount);
}

Distribution::Distribution(Species* species, int maxNodes, int seed) {    
    this->nodes = new DistNode*[maxNodes];
    this->species = species;
    
    for (int i = 0; i < maxNodes; i++) {
        this->nodes[i] = new DistNode(-1, 0, 0);
    }
    
    this->maxNodes = maxNodes;
    this->numNodes = 0;
    
    this->rng = new std::mt19937(seed);
}

Distribution::~Distribution() {
    for (int i = 0; i < this->maxNodes; i++) {
        delete this->nodes[i];
    }
    
    this->species = NULL;
    delete[] this->nodes;
    delete this->rng;
}

void Distribution::addNode(double state, int count) {
    this->nodes[this->numNodes]->state = state;
    this->nodes[this->numNodes]->count = count;
    
    if (this->numNodes == 0) {
        this->nodes[this->numNodes]->cumCount = count;
    } else {
        this->nodes[this->numNodes]->cumCount = count + this->nodes[this->numNodes - 1]->cumCount;
    }
    
    this->numNodes++;
}

void Distribution::clear() {
    this->numNodes = 0;
}

void Distribution::update(double** statePt, int length, bool bounded) {
    this->numNodes = 0;
    
    std::vector<double> states;
    for (int i = 0; i < length; i++) {
        states.push_back(statePt[i][this->species->id]);
    }
    
    std::sort(states.begin(), states.end());
    
    int startIndex = 0;
    if (bounded) {
        for (int i = 0; i < length; i++) {
            if (states[i] > this->species->stateLowerBound) {
                startIndex = i;
                break;
            }
        }
    }
    
    int count = 0;
    double lastState = states[startIndex];
    for (int i = startIndex; i < length; i++) {
        if (bounded && states[i] >= this->species->stateUpperBound) {
            break;
        } else if (states[i] == lastState) {
            count++;
        } else {
            this->addNode(lastState, count);
            
            lastState = states[i];
            count = 1;
        }
    }
    
    this->addNode(lastState, count);
}

int Distribution::getNumNodes() const {
    return this->numNodes;
}

double Distribution::getState(int nodeId) const {
    return this->nodes[nodeId]->state;
}

int Distribution::getCount(int nodeId) const {
    return this->nodes[nodeId]->count;
}

double Distribution::calcDistance(Distribution* other) const {
    double distance = 0;
    
    int i = 0;
    int j = 0;
    
    int numNodes = this->numNodes;
    int otherNumNodes = other->getNumNodes();
    
    double state = this->nodes[i]->state;
    double otherState = other->getState(j);
    
    int count = this->nodes[i]->count;
    int otherCount = other->getCount(j);
    
    while (i != numNodes || j != otherNumNodes) {        
        if (state < otherState) {
            distance += count * count;
            i++;            
            
            state = this->nodes[i]->state;
            count = this->nodes[i]->count;
        } else if (state > otherState) {
            distance += otherCount * otherCount;
            j++;
            
            otherState = other->getState(j);
            otherCount = other->getCount(j);
        } else {
            distance += (count - otherCount) * (count - otherCount);
            i++;
            j++;
            
            state = this->nodes[i]->state;
            count = this->nodes[i]->count;
            
            otherState = other->getState(j);
            otherCount = other->getCount(j);
        }
    }
    
    for (int k = i; k < numNodes; k++) {
        count = this->nodes[k]->count;
        distance += count * count;
    }
    
    for (int k = j; k < otherNumNodes; k++) {        
        otherCount = other->getCount(k);
        distance += otherCount * otherCount;
    }
    
    return sqrt((double) distance);
}

double Distribution::sample() const {
    double count = 1.0 * (*this->rng)() / this->rng->max() * this->maxNodes;
    
    int i;
    for (i = 0; i < this->numNodes; i++) {
        if (count < this->nodes[i]->cumCount) {
            break;
        }
    }
    
    return this->nodes[i]->state;
}

void Distribution::print() const {
    for (int i = 0; i < this->numNodes; i++) {
        this->nodes[i]->print();
    }
}