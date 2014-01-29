#include "PriorityQueue.h"
#include "Reaction.h"
#include <float.h>
#include <math.h>
#include <stdio.h>

using namespace std;

Node::Node(Reaction* rxn, double time) {
    this->rxn = rxn;
    this->time = time;
}

Node::Node(const Node& other) {
    this->rxn = other.rxn;
    this->time = other.time;
}

Node::~Node() {
    this->rxn = NULL;
}

PriorityQueue::PriorityQueue(int maxNodes, bool minHeap) {
    this->minHeap = minHeap;
    
    this->nodes = new Node*[maxNodes];
    
    for (int i = 0; i < maxNodes; i++) {
        this->nodes[i] = new Node(0, DBL_MAX);
    }
    
    this->maxNodes = maxNodes;
    this->numNodes = 0;
}

PriorityQueue::PriorityQueue(const PriorityQueue& other) {
    this->minHeap = other.minHeap;
    
    this->nodes = new Node*[other.maxNodes];
       
    this->maxNodes = other.maxNodes;
    this->numNodes = other.numNodes;
    
    for (int i = 0; i < this->maxNodes; i++) {
        this->nodes[i] = new Node(*(other.nodes[i]));
    }
}

PriorityQueue::~PriorityQueue() {
    for (int i = 0; i < maxNodes; i++) {
        delete this->nodes[i];
    }
    
    delete[] this->nodes;
}

void PriorityQueue::addNode(Reaction* rxn, double time) {
    this->nodes[this->numNodes]->rxn = rxn;
    this->nodes[this->numNodes]->time = time;
    
    this->numNodes++;
    
    this->heapify(this->numNodes - 1);
}

void PriorityQueue::clearQueue() {
    this->numNodes = 0;
}

void PriorityQueue::heapSort() {
    for (int i = this->numNodes - 1; i >= 0; i--) {
        this->heapify(i);
    }
}

void PriorityQueue::heapify(int nodeId) {
    int parentId = (int) floor((nodeId + 1.0) / 2.0) - 1;
    if (parentId < 0) {
        parentId = 0;
    }
    
    int leftId = (nodeId + 1) * 2 - 1;
    int rightId;
        
    if (leftId >= this->numNodes) {
        leftId = nodeId;
        rightId = nodeId;
    } else {
        rightId = leftId + 1;
        if (rightId >= this->numNodes) {
            rightId = nodeId;
        } 
    }
    
    double time = this->nodes[nodeId]->time;
    double parentTime = this->nodes[parentId]->time;
    double leftTime = this->nodes[leftId]->time;
    double rightTime = this->nodes[rightId]->time;
    
    if (this->minHeap) {
        if (time < parentTime) {
            this->swap(nodeId, parentId);
            this->heapify(parentId);
        } else if (time > leftTime) {
            this->swap(nodeId, leftId);
            this->heapify(leftId);
        } else if (time > rightTime) {
            this->swap(nodeId, rightId);
            this->heapify(rightId);
        }
    } else {
        if (time > parentTime) {
            this->swap(nodeId, parentId);
            this->heapify(parentId);
        } else if (time < leftTime) {
            this->swap(nodeId, leftId);
            this->heapify(leftId);
        } else if (time < rightTime) {
            this->swap(nodeId, rightId);
            this->heapify(rightId);
        }
    }
}

void PriorityQueue::swap(int nodeId1, int nodeId2) {
    Node* tempNodePtr = this->nodes[nodeId1];
    this->nodes[nodeId1] = this->nodes[nodeId2];
    this->nodes[nodeId2] = tempNodePtr;
}

int PriorityQueue::getNextRxnId() const {
    return this->nodes[0]->rxn->id;
}

Reaction* PriorityQueue::getNextRxn() const {
    return this->nodes[0]->rxn;
}

double PriorityQueue::getNextTime() const {
    return this->nodes[0]->time;
}

int PriorityQueue::search(int rxnId) const{
    int nodeId = -1;
    
    for (int i = 0; i < this->numNodes; i++) {
        if (this->nodes[i]->rxn->id == rxnId) {
            nodeId = i;
            break;
        }
    }
    
    return nodeId;
}

void PriorityQueue::updatePqByRxnId(int rxnId, double newTime) {
    int nodeId = this->search(rxnId);    
    this->updatePqByNodeId(nodeId, newTime);
}

void PriorityQueue::updatePqByNodeId(int nodeId, double newTime) {
    this->updateNodeByNodeId(nodeId, newTime);
    this->heapify(nodeId);
}

void PriorityQueue::updateNodeByNodeId(int nodeId, double newTime) {
    this->nodes[nodeId]->time = newTime;
}

void PriorityQueue::updateRxnByNodeId(int nodeId, Reaction* newRxn) {
    this->nodes[nodeId]->rxn = newRxn;
}

double PriorityQueue::getTimeByNodeId(int nodeId) const {
    return this->nodes[nodeId]->time;
}

int PriorityQueue::getRxnIdByNodeId(int nodeId) const {
    return this->nodes[nodeId]->rxn->id;
}

Reaction* PriorityQueue::getRxnByNodeId(int nodeId) const {
    return this->nodes[nodeId]->rxn;
}

int PriorityQueue::getNumNodes() const {
    return this->numNodes;
}

void PriorityQueue::print() const {
    for (int i = 0; i < this->numNodes; i++) {
        fprintf(stdout, "Reaction ID: %i, Time: %e\n", getRxnIdByNodeId(i), getTimeByNodeId(i));
    }
}