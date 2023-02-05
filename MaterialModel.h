#ifndef _MATMOD_H_
#define _MATMOD_H_
#include "Mesh.h"
#include <armadillo>
#include <vector>

class MaterialModel {
//Attributes
private:
    unsigned int nn; //init list
    
    unsigned int n_lvls;
    arma::colvec* lvls;
    arma::colvec* velocities;
    arma::colvec* control_func;

//Methods
public:
    MaterialModel(unsigned int nnodes, double init_value);
    ~MaterialModel();
    
    void set(unsigned int n_levels, std::vector<double> levels, std::vector<double> velocities_arr);
    
    arma::colvec* get_levels();
    arma::colvec* get_velocities();
    arma::colvec* get_control_function();
    unsigned int get_nlvls();
    unsigned int get_nvels();
    
    void log_control_func();
    void log_levels();
    void log_velocities();

};

#endif // _MATMOD_H_
