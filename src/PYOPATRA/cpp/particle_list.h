//
// Created by Georgia Stuart on 5/11/21.
//

#ifndef PYOPATRA_PARTICLE_LIST_H
#define PYOPATRA_PARTICLE_LIST_H

#include "list.h"
#include "particle.h"

template <int dim>
class ParticleList{
public:
    using ParticleN = ParticleBase<dim>;
    using Vector = typename ParticleN::Vector;

    void create_particle(const Vector& location) {
        ParticleN* temp = new ParticleN();
        temp->set_location(location);
        list.push(temp->get_node());
    }
    int get_num_particles() { return list.length; }
    Eigen::MatrixXd get_all_particle_locations() const {
        Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(list.length, dim);
        auto current = list.get_head();
        size_t index = 0;
        while (current) {
            temp.row(index++) = current->owner->get_location();
            current = current->next;
        }
        return temp;
    }

private:
    List<ParticleN> list;
};

using ParticleList2D = ParticleList<2>;

#endif //PYOPATRA_PARTICLE_LIST_H
