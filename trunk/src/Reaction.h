#ifndef REACTION_H
#define REACTION_H

class Species;

class Reaction{
public:
    const int id;
    
    const int numStoichSpecies;
    Species** stoichSpecies;
    const int* const stoichCoeffs;
    
    const int rateLaw;
    const int numRateConsts;
    
    const int numRateSpecies;
    Species** rateSpecies;
    
    const int numDeps;
    const int* const deps;
    
    double prop;
    double oldProp;
    
    Reaction(const int id, Species** species, const int numStoichSpecies, int* stoichSpeciesIds, const int* const stoichCoeffs, const int rateLaw, const int numRateConsts, double* rateConsts, const int numRateSpecies, int* rateSpeciesIds, const int numDeps, const int* const deps, double vol);
    Reaction(const Reaction& other);
    virtual ~Reaction();
    
    void updateProp(double volRatio);
    
    void print() const;
private:    
    double* rateConsts;
};

#endif