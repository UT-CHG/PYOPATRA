//
// Created by Georgia Stuart on 2/23/21.
//

#ifndef PYOPATRA_ILLNODE_H
#define PYOPATRA_ILLNODE_H

template <class T>
class ILLNode {
public:
    ILLNode *next, *prev;
    T *owner;

    ILLNode() : next(this), prev(this), owner(nullptr) {};
    explicit ILLNode(T *owner) : next(this), prev(this), owner(owner) {};

    void remove();
    void insert_after(ILLNode<T> &node);


};

#include "illnode.inl"

#endif //PYOPATRA_ILLNODE_H
