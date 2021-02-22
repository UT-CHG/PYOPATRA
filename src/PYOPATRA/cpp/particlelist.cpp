//
// Created by Georgia Stuart on 2/3/21.
//

#include <iostream>
#include "particlelist.h"

ParticleList::ParticleList()
    : head(nullptr)
    , tail(nullptr)
{}

ParticleList::ParticleList(Node *head)
    : head(head)
    , tail(find_tail(head))
{}

Node * ParticleList::find_tail(Node *start) {
    Node *cur = start;
    while (cur->next != nullptr) {
        cur = cur->next;
    }

    return cur;
}

Node& ParticleList::operator[](int index) const {
    Node *cur = head;
    for (int i = 0; i < index; i++) {
        cur = cur->next;

        if (cur == nullptr) {
            throw std::out_of_range("Index out of range.");
        }
    }

    return *cur;
}

void ParticleList::insert(double latitude, double longitude, double diameter, double density, double depth) {
    tail->next = new Node(latitude, longitude, diameter, density, depth);
    tail->next->previous = tail;
    tail = tail->next;
}

void ParticleList::insert(int index, double latitude, double longitude, double diameter, double density, double depth) {
    Node *new_node = new Node(latitude, longitude, diameter, density, depth);
    Node *before = &(*this)[index];
    Node *after = before->next;

    new_node->previous = before;
    new_node->next = after;

    before->next = new_node;
    after->previous = new_node;
}