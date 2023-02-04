#include <iostream>
#include <vector>
#include "Mesh.h"
#include "Device.h"
#include "MaterialModel.h"
#include "Problem.h"
#include <armadillo>


namespace setup {
    const unsigned int NEL {50}; 
    const unsigned int NED {50};
    const double DOMAIN_LENGTH {2}; 
    const double DOMAIN_DEPTH {2};

    const double INIT_CONTROL {0.5};
    double LEVELS [] {0.0, 0.5, 1.0, 1.5};
    unsigned int N_LEVELS {4};
    double VELS [] {0.5, 1.0, 3.5, 4.5};

    const double PULSE_INTENSITY {1.00};
    const double PULSE_CENTRAL_FREQUENCY {2.0};
    const double OBSERVATION_PERIOD {2.6};
    const double TIME_STEP {0.002};
}


int main(){
    
    //Creating problem:
    Problem real_problem {
                          setup::NEL, 
                          setup::NED,
                          setup::DOMAIN_LENGTH,
                          setup::DOMAIN_DEPTH, 
                          setup::PULSE_INTENSITY, 
                          setup::PULSE_CENTRAL_FREQUENCY, 
                          setup::OBSERVATION_PERIOD, 
                          setup::TIME_STEP,
                          setup::INIT_CONTROL
                          }; 
                          
    real_problem.add_device("receiver",1.0,1.0,0);
    real_problem.add_device("source",  1.5,1.5,0);
//    real_problem.add_device("source",  1.0,1.5,0);
//    real_problem.add_device("source",  1.5,1.0,0);
    real_problem.set_control(setup::N_LEVELS, setup::LEVELS, setup::VELS);
    real_problem.build();
    real_problem.solve("forward_only");
    real_problem.write_output();
    
    



    std::cout << " " << std::endl;
    std::cout << "This program ended succesfully. " << std::endl; 
    return 0;
}