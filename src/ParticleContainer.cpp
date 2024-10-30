//
// Created by jerry_park
//

#include "ParticleContainer.h"
#include "Particle.h"
#include <iostream>
#include "utils/ArrayUtils.h"
#include <vector>
#include <array>
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"

ParticleContainer::ParticleContainer(std::vector<Particle> particles,double start_time,
                                     double end_time , double delta_t,std::string outputformat){
  this->particles = particles;
  this->delta_t = delta_t;
  this->start_time = start_time;
  this->end_time = end_time;
  this->outputformat = outputformat;
  std::cout <<"ParticleContainer generated!\n";
 }
//ParticleContainer::ParticleContainer(){}
//void ParticleContainer:: setParticles(std::vector<Particle> particles){
//  this->particles = particles;
//}

ParticleContainer::~ParticleContainer(){
  std::cout << "ParticleContainer destructed!\n";
}

void ParticleContainer::addParticle(Particle particle){
  this->particles.push_back(particle);
}

void ParticleContainer :: calculateF() {
  std::vector<Particle>::iterator iterator;
  iterator = particles.begin();

  for (auto &p1 : particles) {
    p1.setOldF(p1.getF()); // change Old force and new Force

    std::array<double, 3> cal_f = {0,0,0};

    for (auto &p2 : particles) {
      if(&p1 != &p2) {
        std::array<double, 3> vec = p2.getX() - p1.getX();

        double tmp = 0.0;

        for(int i = 0; i < 3; i++) {
          tmp += vec[i] * vec[i];
        }
        tmp = pow(sqrt(tmp), 3);

        if (tmp == 0 ) continue;// divider must not be 0

        double val = p1.getM() * p2.getM();
        val /= tmp;

        for(int i=0 ; i< 3; i++) {
          cal_f[i] += (vec[i]* val);
        }
        // std :: array<double, 3> x_vector = p1.x;
        // @TODO: insert calculation of forces here!
      }
      p1.setF(cal_f);
    }

  }
   // std::cout << "force: " << p1 << "\n";
}

void ParticleContainer:: calculateX(){
   for (auto &p : particles) {
    std::array<double, 3> vec = p.getX();

    for(int i =0 ; i< 3; i++) {
      vec[i] += delta_t * p.getV()[i];
      vec[i] += (delta_t * delta_t * p.getOldF()[i]) / (2 * p.getM());
    }
    p.setX(vec);
    //std::cout << "location: " << p << "\n";
    // @TODO: insert calculation of position updates here!
  }
}


void ParticleContainer:: calculateV(){
  for (auto &p : particles) {
    std::array<double, 3> vec = p.getV();
    for(int i =0 ; i< 3; i++) {
      vec[i] += delta_t * (p.getF()[i] + p.getOldF()[i]) / (2 * p.getM());
    }
    p.setV(vec);
    //std::cout << "v : "<< p << "\n";
    // @TODO: insert calculation of veclocity updates here!
  }
}

void ParticleContainer:: calculate(){
   int iteration  = 0;
   double current_time = start_time;
   while (current_time < end_time) {
    // calculateX();
    //
    // calculateF();
    // calculateV();

    // calculate new x
    calculateX();
    // calculate new f
    calculateF();
    // calculate new v
    calculateV();

    iteration++;
    if (iteration % 10 == 0) {
      plotParticles(iteration); // temporary
    }
    std::cout << "Iteration " << iteration << " finished." << "\n";
    current_time += delta_t;
    }
    std::cout << "output written. Terminating..." << "\n";

}

void ParticleContainer :: plotParticles(int iteration) {
  std::cout << "Plotting Particles..." << "\n";

  std::string out_name("MD_vtk");
  if(outputformat.compare(".vtu") == 0) {
    outputWriter::VTKWriter writer2;
    writer2.writeFile(out_name, iteration,particles);
  }
  else if(outputformat.compare(".xyz") == 0) {
    outputWriter::XYZWriter writer;
    writer.plotParticles(particles, out_name, iteration);
  }

}
