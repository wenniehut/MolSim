
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

  //Destructor  
   ~ParticleContainer();

  
    void addParticle(Particle particle);

    //calculate the force for all particles
    void calculateF();

     //calculate the position for all particles
    void calculateX();

    //calculate the velocity for all particles
    void calculateV();

    //calculate the position, force and velocity for all particles
    void calculate();

    //Plots the particles in the container
    void plotParticles(int iterations);

  
  private:
    //Member variables 
    std::vector<Particle> particles;
    std::vector<std::array<double,3>> positions;
    double start_time = 0;
    double end_time = 1000;
    double delta_t;
    std::string outputformat;


};

#endif //PARTICLECONTAINER_H
