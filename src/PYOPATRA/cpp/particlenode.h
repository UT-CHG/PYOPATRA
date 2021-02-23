//
// Created by georgia on 2/23/21.
//

#ifndef PYOPATRA_PARTICLENODE_H
#define PYOPATRA_PARTICLENODE_H

// Forward Declaration of Particle
class Particle;

class ParticleNode {
private:
    Particle *next, *prev;
};

#endif //PYOPATRA_PARTICLENODE_H
