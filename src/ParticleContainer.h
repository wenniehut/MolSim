
#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

#include "Particle.h"
#include <vector>
#include <array>

#pragma once


/**
 * @brief ParticleContainer class
 */
class ParticleContainer {

  public:

   ParticleContainer(std::vector<Particle> particles, double start_time, double end_time, double delta_t,
    std::string outputformat);

  /**
   * @brief Destructor
   */
   ~ParticleContainer();

  
    /**
     * @brief Add a particle to the container
     */
    void addParticle(Particle particle);

    /**
     * @brief Calculate the force between all particles version1
     */
    void calculateF_v1();

    /**
     * @brief Calculate the force between all particles version2 , optimized version by memory access
     */
    void calculateF_v2();

     /**
      * @brief Calculate the position for all particles
      */
    void calculateX();

    /**
     * @brief Calculate the velocity for all particles
     */
    void calculateV();

    /**
     * @brief Calculate the position, force and velocity for all particles
     * @param version 1 or 2
     */
    void calculate(int version);

    /** 
     * @brief Plot the particles in the container
     */
    void plotParticles(int iterations);

  
  private:
    /**
     * @brief Member variables
     */

    std::vector<Particle> particles;
    std::vector<std::array<double,3>> positions;
    double start_time = 0;
    double end_time = 1000;
    double delta_t;
    std::string outputformat;


};

#endif //PARTICLECONTAINER_H
