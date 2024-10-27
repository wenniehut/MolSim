
#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"

#include <iostream>
#include <list>


#include <cmath> // here
#include "outputWriter/VTKWriter.h" // here

/**** forward declaration of the calculation functions ****/

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

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

constexpr double start_time = 0;
constexpr double end_time = 1000;
constexpr double delta_t = 0.014;

// TODO: what data structure to pick?
std::list<Particle> particles;

int main(int argc, char *argsv[]) {

  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc != 2) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename" << std::endl;
  }

  FileReader fileReader;
  fileReader.readFile(particles, argsv[1]);

  double current_time = start_time;

  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    // calculate new x
    calculateX();
    // calculate new f
    calculateF();
    // calculate new v
    calculateV();

    iteration++;
    if (iteration % 10 == 0) {
      plotParticles(iteration);
    }
    std::cout << "Iteration " << iteration << " finished." << std::endl;

    current_time += delta_t;
  }

  std::cout << "output written. Terminating..." << std::endl;
  return 0;
}


void calculateF() {
  std::list<Particle>::iterator iterator;
  iterator = particles.begin();

  for (auto &p1 : particles) {
    p1.setOldF(p1.getF()); // change Old force and Force
    std::array<double, 3> cal_f = {0,0,0};
    for (auto &p2 : particles) {

      if(&p1 == &p2) {continue;}

      std::array<double, 3> vec = p2.getX() - p1.getX();

      double val = p1.getM() * p2.getM();
      double tmp = 0.0;
      for(int i = 0; i < 3; i++) {
        tmp += vec[i] * vec[i];
      }
      tmp = pow(tmp, 3/2);
      val /= tmp;

      for(int i =0 ; i< 3; i++) {
        vec[i] *= val;
      }
      for(int i=0 ; i< 3; i++) {
        cal_f[i] += vec[i];
      }
      // std :: array<double, 3> x_vector = p1.x;
      // @TODO: insert calculation of forces here!
    }
    p1.setF(cal_f);
  }
}
void calculateX() {
  for (auto &p : particles) {
    std::array<double, 3> vec = p.getX();
    for(int i =0 ; i< 3; i++) {
      vec[i] += delta_t * p.getV()[i];
      vec[i] += (delta_t * delta_t * p.getF()[i]);
      vec[i] /= 2 * p.getM();
    }
    p.setX(vec);
    std::cout << p << std::endl;
    // @TODO: insert calculation of position updates here!
  }
}

void calculateV() {
  for (auto &p : particles) {
    std::array<double, 3> vec = p.getV();
    for(int i =0 ; i< 3; i++) {
      vec[i] += delta_t * (p.getF()[i] + p.getOldF()[i]) / (2 * p.getM());
    }
    p.setV(vec);
    // @TODO: insert calculation of veclocity updates here!
  }
}

void plotParticles(int iteration) {

  std::string out_name("MD_vtk");

  outputWriter::XYZWriter writer;
  writer.plotParticles(particles, out_name, iteration);

  outputWriter::VTKWriter writer2;
  writer2.writeFile(out_name, iteration,particles);


}
