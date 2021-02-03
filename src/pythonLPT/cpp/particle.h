//
// Created by Georgia Stuart on 2/3/21.
//

#ifndef LPTCPP_PARTICLE_H
#define LPTCPP_PARTICLE_H

class Particle {
public:
    double latitude, longitude, diameter, density, depth;

    Particle();
    Particle(double latitute, double longitude, double diameter, double density, double depth);
};

#endif //LPTCPP_PARTICLE_H
