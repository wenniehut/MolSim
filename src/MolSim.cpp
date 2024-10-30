
#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"

#include <iostream>
#include <list>

///////////////////////////////////
#include <cmath>
#include "outputWriter/VTKWriter.h"
#include <getopt.h>
#include <cstdlib>
#include <float.h>
#include <unistd.h>
//////////////////////////////////

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

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

void printHelp();
// constexpr double start_time = 0;
// constexpr double end_time = 1000;
// constexpr double delta_t = 0.014;
double start_time = 0;
double end_time = 1000;
double delta_t = 0.014;
std::string outputformat = ".vtu";

// TODO: what data structure to pick?

//std::list<Particle> particles;
std::vector<Particle> particles;


int main(int argc, char *argsv[]) {

  std::cout << "Hello from MolSim for PSE!" << "\n"
  << "To run this program, please provide the file name as an argument, like this: `./Molsim abc.txt " << "\n"
  << "To see more options, type ./Molsim help or ./Molsim --help" << "\n\n";

  for (int i = 0; i < argc; i++) {
    if (std::string(argsv[i]) == "--help" || std::string(argsv[i]) == "help") {
      printHelp();
      return 0;
    }
  }

  if (argc < 2) {
    std::cout << "Error: Filename is required." << "\n";
    printHelp();
    return 1;
  }

  FileReader fileReader;
  fileReader.readFile(particles, argsv[1]);
  
  if (argc >= 3) delta_t = std::atof(argsv[2]);
  if (argc >= 4) end_time = std::atof(argsv[3]);
  if (argc >= 5) {
    outputformat = argsv[4];
    if (outputformat != ".xyz" && outputformat != ".vtu") {
      std::cout << "Invalid output format! Choose either 'xyz' or 'vtu'." << "\n\n";
      return 1;
    }
  }

  std::cout << "Testfilename: " << argsv[1] << "\n";
  std::cout << "Start Time: " << start_time << "\n";
  std::cout << "Time End: " << end_time << "\n" ;
  std::cout << "Delta Time: " << delta_t << "\n";
  std::cout << "Output format: " << outputformat << "\n\n";

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
      plotParticles(iteration); // temporary
    }
    std::cout << "Iteration " << iteration << " finished." << "\n";

    current_time += delta_t;

  }

  std::cout << "output written. Terminating..." << "\n";
  return 0;
}


void calculateF() {
  std::vector<Particle>::iterator iterator;
  iterator = particles.begin();

  for (auto &p1 : particles) {
    p1.setOldF(p1.getF()); // change Old force and Force

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
void calculateX() {
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

void calculateV() {
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

void plotParticles(int iteration) {
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

void printHelp() {
  std::cout << "MolSim for PSE - Usage:" << "\n";
  std::cout << "./MolSim filename [delta_t] [end_time] [output_format]" << "\n";
  std::cout << "Options:" << "\n";
  std::cout << "  1.  Specify the input file name (required)" << "\n";
  std::cout << "  2.  Set the delta time (optional, default is 0.014)" << "\n";
  std::cout << "  3.  Set the end time (optional, default is 1000)" << "\n";
  std::cout << "  4.  Specify output format: either \".xyz\" or \".vtu\" (optional, default is \".vtu \")" << "\n";
  std::cout << "  --help, help message" << "\n\n";
}