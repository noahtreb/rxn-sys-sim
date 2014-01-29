#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H

class Reaction;

class Node{
public:
    Reaction* rxn;
    double time;
    
    Node(Reaction* rxn, double time);
    Node(const Node& other);
    virtual ~Node();
private:
};

class PriorityQueue{
public:        
    bool minHeap;
    
    PriorityQueue(int maxNodes, bool minHeap);
    PriorityQueue(const PriorityQueue& other);
    virtual ~PriorityQueue();
    
    void addNode(Reaction* rxn, double time);
    void clearQueue();
    
    void heapSort();
    int search(int rxnId) const;
    
    int getNextRxnId() const;
    Reaction* getNextRxn() const;
    double getNextTime() const;
    
    void updatePqByRxnId(int rxnId, double newTime);
    void updatePqByNodeId(int nodeId, double newTime);
    
    void updateNodeByNodeId(int nodeId, double newTime);
    void updateRxnByNodeId(int nodeId, Reaction* newRxn);
    
    double getTimeByNodeId(int nodeId) const;
    int getRxnIdByNodeId(int nodeId) const;
    Reaction* getRxnByNodeId(int nodeId) const;
    
    int getNumNodes() const;
    
    void print() const;
private:        
    int numNodes;
    int maxNodes;
    Node** nodes;
    
    void heapify(int nodeId);
    void swap(int nodeId1, int nodeId2);
};

#endif