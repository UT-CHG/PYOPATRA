//
// Created by Georgia Stuart on 2/3/21.
//

#ifndef PYOPATRA_PARTICLELIST_H
#define PYOPATRA_PARTICLELIST_H

#include "particle.h"

class Node {
public:
    Node *next, *previous;
    Particle particle;

    Node(double latitude, double longitude, double diameter, double density, double depth);
};

class ParticleList {
public:
    Node *head, *tail;


    ParticleList();
    ParticleList(Node *head);
    Node* find_tail(Node *start);
    void insert(double latitude, double longitude, double diameter, double density, double depth);
    void insert(int index, double latitude, double longitude, double diameter, double density, double depth);
    void remove(int index);
    Node& operator[](int index) const;
};

#endif //PYOPATRA_PARTICLELIST_H
