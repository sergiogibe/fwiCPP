#include "MaterialModel.h"
#include <iostream>

MaterialModel::MaterialModel(unsigned int nnodes, double init_value) : nn(nnodes) {
    control_func = new double[nn];
    for (unsigned int i {0}; i < nn; i++) {
        control_func[i] = init_value;
    }
}

MaterialModel::~MaterialModel() {
    delete [] control_func;
    delete [] lvls;
    delete [] velocities;
}

void MaterialModel::set_levels(unsigned int n_levels, double* levels) {
    n_lvls = n_levels;
    lvls = new double[n_lvls];
    lvls = levels;
}

void MaterialModel::set_velocities(double* velocities_arr) {
    velocities = new double[n_lvls];
    velocities = velocities_arr;
}

double* MaterialModel::get_levels() {
    return lvls;
}

double* MaterialModel::get_velocities() {
    return velocities;
}

double* MaterialModel::get_control_function() {
    return control_func;
}

unsigned int MaterialModel::get_nlvls() {
    return n_lvls;
}

unsigned int MaterialModel::get_nvels() {
    return n_lvls;
}

void MaterialModel::log_control_func() {
    if (nn < 150) {
        std::cout << "------------------" << std::endl;
        for (unsigned int i {0}; i < nn; i++) {
            std::cout << control_func[i] << std::endl;
        }
        std::cout << "------------------" << std::endl;
    }
}

void MaterialModel::log_levels() {
    std::cout << "--------------------------------------------------------------------" << std::endl;
    for (unsigned int i {0}; i < n_lvls; i++) {
        std::cout << "    Level " << i+1 << ": " << lvls[i];
    }
    std::cout << "\n--------------------------------------------------------------------" << std::endl;
}

void MaterialModel::log_velocities() {
    std::cout << "--------------------------------------------------------------------" << std::endl;
    for (unsigned int i {0}; i < n_lvls; i++) {
        std::cout << "    Vel " << i+1 << ": " << velocities[i];
    }
    std::cout << "\n--------------------------------------------------------------------" << std::endl;
}
