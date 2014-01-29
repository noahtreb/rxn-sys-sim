#ifndef SPECIES_H
#define SPECIES_H

class Species {
public:
    const int id;
    double state;
    const bool stateChanges;
    const double stateUpperBound;
    const double stateLowerBound;
    const bool stateBounded;
    
    Species(const int id, double state, const bool stateChanges, const double stateUpperBound, const double stateLowerBound, const bool stateBounded);
    Species(const Species& other);
    virtual ~Species();
private:
};

#endif