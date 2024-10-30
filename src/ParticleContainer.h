//
// Created by jerry_park on 10/30/24.
//

#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

#include "Particle.h"
#include <vector>
#include <array>

#pragma once


class ParticleContainer {

  public:

   ParticleContainer(std::vector<Particle> particles, double start_time, double end_time, double delta_t,
    std::string outputformat);

   ParticleContainer();

   virtual ~ParticleContainer();


    void addParticle(Particle particle);
    /**** forward declaration of tㅡhe calculation functions ****/
    /**
     * calculate the force for all particles
     */
    void calculateF();

    /**
     * calculate the position for all particles
     */
    void calculateX();

    /**
     * calculate the position for all particles
     */
    void calculateV();

    void calculate();

    void plotParticles(int iterations);

  private:
    double start_time = 0;
    double end_time = 1000;
    double delta_t;
    std::string outputformat;
    std::vector<Particle> particles;
    std::vector<std::array<double,3>> positions;

};



#endif //PARTICLECONTAINER_H
