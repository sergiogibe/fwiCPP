#include <iostream>
#include <fstream>
#include <vector>
#include "Mesh.h"
#include "Device.h"
#include "MaterialModel.h"
#include "Problem.h"
#include <armadillo>
#include <json.hpp>


int main(){
    
    //Reading data:
    std::ifstream file("settings.json");
    nlohmann::json settings;
    file >> settings;
    
    //Creating problem:
    Problem real_problem {
                            settings["mesh"]["nel"],
                            settings["mesh"]["ned"],
                            settings["domain"]["length"],
                            settings["domain"]["depth"],
                            settings["time"]["pulse_int"],
                            settings["time"]["pulse_freq"],
                            settings["time"]["period"],
                            settings["time"]["dt"],
                            settings["levelset"]["init"]
                         };
                          
    real_problem.add_device("receiver", 1.0, 1.0, 0);
    real_problem.add_device("source",   1.5, 1.5, 0);
    real_problem.add_device("source",   1.5, 1.5, 1);
    real_problem.add_device("source",   1.5, 1.5, 2);
//    real_problem.add_device("source",  1.0,1.5,0);
//    real_problem.add_device("source",  1.5,1.0,0);

    real_problem.set_control(settings["levelset"]["n_lvls"], 
                             settings["levelset"]["levels"], 
                             settings["levelset"]["vels"]);
    
    real_problem.build();
    real_problem.solve("forward_only");
    real_problem.write_output(settings["output"]["save_solution"],
                              settings["output"]["shot_id"], 
                              settings["output"]["sample_size"]);


    std::cout << " " << std::endl;
    std::cout << "This program ended succesfully. " << std::endl; 
    return 0;
}