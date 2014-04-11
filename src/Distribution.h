#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <random>

class Species;

class DistNode {
public:
    double state;
    int count;
    int cumCount;
    
    DistNode(double state, int count, int cumCount);
    virtual ~DistNode();
    
    void print() const;
private:
};

class Distribution {
public:
    Species* species;
    
    Distribution(Species* species, int maxNodes, int seed);
    virtual ~Distribution();   
    
    void addNode(double state, int count);    
    void clear();
    
    void update(double** statePt, int length, bool bounded);
    
    int getNumNodes() const;
    double getState(int nodeId) const;
    int getCount(int nodeId) const;
    
    double calcDistance(Distribution* other) const;
    double sample() const;
    
    void print() const;
private:
    int numNodes;
    int maxNodes;
    DistNode** nodes;
    
    std::mt19937* rng;
};
#endif