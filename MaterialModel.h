#ifndef _MATMOD_H_
#define _MATMOD_H_
#include "Mesh.h"

class MaterialModel {
//Attributes
private:
    unsigned int nn; //init
    
    unsigned int n_lvls;
    double* lvls;
    double* velocities;
    double* control_func;

//Methods
public:
    MaterialModel(unsigned int nnodes, double init_value);
    ~MaterialModel();
    
    void set_levels(unsigned int n_levels, double* levels);
    void set_velocities(double* velocities_arr);
    
    double* get_levels();
    double* get_velocities();
    double* get_control_function();
    unsigned int get_nlvls();
    unsigned int get_nvels();
    
    void log_control_func();
    void log_levels();
    void log_velocities();

};

#endif // _MATMOD_H_
