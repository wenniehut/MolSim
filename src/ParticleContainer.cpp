#include <iostream>
#include <vector>
#include <array>
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"
#include "Particle.h"
#include "ParticleContainer.h"


//Constructor, which initializes the particle container with a vector of particles , start time, end time, delta_t and outputformat
//outputformat can be either .vtu or .xyz
ParticleContainer::ParticleContainer(std::vector<Particle> particles,double start_time,
                                     double end_time , double delta_t,std::string outputformat){

  this->particles = particles;
  this->delta_t = delta_t;
  this->start_time = start_time;
  this->end_time = end_time;
  this->outputformat = outputformat;

  positions.reserve(particles.size());
  for(int i = 0 ; i< particles.size(); i++) {
    positions[i] = particles[i].getX();
  }

  std::cout << "ParticleContainer generated!\n";
 }


//Destructer 
ParticleContainer::~ParticleContainer(){
  std::cout << "ParticleContainer destructed!\n";
}

//In case the user wants to add a particle to the container
void ParticleContainer::addParticle(Particle particle){
  this->particles.push_back(particle);
  this->positions.push_back(particle.getX());
}

void ParticleContainer :: calculateF_v1() {
  std::vector<Particle>::iterator iterator;
  iterator = particles.begin();
  for (auto &p1 : particles) {
    p1.setOldF(p1.getF()); // update old force
    std::array<double, 3> cal_f = {0,0,0};
    for (auto &p2 : particles) {
      if(&p1 != &p2) {
        std::array<double, 3> vec = p2.getX() - p1.getX();
        double tmp = 0.0;
        for(int i = 0; i < 3; i++) {
          tmp += vec[i] * vec[i];
        }
        tmp = pow(sqrt(tmp), 3);
        if (tmp == 0 ) continue;// avoid division by zero
        double val = p1.getM() * p2.getM();
        val /= tmp;
        for(int i=0 ; i< 3; i++) {
          cal_f[i] += (vec[i]* val);
        }
        // @TODO: insert calculation of forces here!
      }
      p1.setF(cal_f);
    }
  }
}


//Calculates the force for all particles
void ParticleContainer :: calculateF_v2() {

  for (int i = 0; i < particles.size(); i++) {
    particles[i].setOldF(particles[i].getF()); // update old force

    std::array<double, 3> cal_f = {0,0,0};

    for (int j =0 ; j< particles.size(); j++) {

      if(i != j) { // avoid self interaction as in the formula
        std::array<double, 3> vec;

        double tmp = 0.0;

        for(int k = 0; k < 3; k++) {
          vec[k] = positions[j][k] - positions[i][k];
          tmp += vec[k] * vec[k];
        }
        tmp = pow(sqrt(tmp), 3);

        if (tmp == 0 ) continue;// avoid division by zero

        double val = particles[i].getM() * particles[j].getM() / tmp;

        for(int k=0 ; k< 3; k++) {
          cal_f[k] += (vec[k]* val);
        }
        // @TODO: insert calculation of forces here!
      }
    }
    particles[i].setF(cal_f);
  }
}

//Calculates the position for all particles
void ParticleContainer:: calculateX(){
   for (int i = 0; i < particles.size(); i++) {
    std::array<double, 3> vec = particles[i].getX();

    for(int j =0 ; j< 3; j++) {
      vec[j] += delta_t * particles[i].getV()[j];
      vec[j] += (delta_t * delta_t * particles[i].getOldF()[j]) / (2 * particles[i].getM());
    }
    particles[i].setX(vec);
    positions[i] = vec;
  }
  
}

//Calculates the velocity for all particles
void ParticleContainer:: calculateV(){
  for (auto &p : particles) {
    std::array<double, 3> vec = p.getV();
    for(int i =0 ; i< 3; i++) {
      vec[i] += delta_t * (p.getF()[i] + p.getOldF()[i]) / (2 * p.getM());
    }
    p.setV(vec);
  }
  
}

//Calculates the position, force and velocity for all particles
void ParticleContainer:: calculate(int version) {
   int iteration  = 0;
   double current_time = start_time;

   while (current_time < end_time) {

    // calculate new x
    calculateX();
    // calculate new f
    if(version == 1) calculateF_v1();
    else calculateF_v2();
    // calculate new v
    calculateV();

    iteration++;
    if (iteration % 10 == 0) {
      plotParticles(iteration); // plot every 10th iteration
    }

    std::cout << "Iteration " << iteration << " finished." << "\n";
    current_time += delta_t;
    }
    std::cout << "output written. Terminating..." << "\n";
    return;
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
