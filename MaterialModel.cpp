#include "MaterialModel.h"
#include <iostream>
#include <armadillo>
#include <vector>

MaterialModel::MaterialModel(unsigned int nnodes, double init_value) : nn(nnodes) {
    control_func = new arma::colvec(nn);
    for (unsigned int i {0}; i < nn; i++) {
        (*control_func)(i) = init_value;
    }
    
    lvls = nullptr;
    velocities = nullptr;
}

MaterialModel::~MaterialModel() {
    delete control_func;
    delete lvls;
    delete velocities;
}

void MaterialModel::set(unsigned int n_levels, std::vector<double> levels, std::vector<double> velocities_arr) {
    n_lvls = n_levels;
    
    delete lvls;
    lvls = new arma::colvec(n_lvls);
    for (arma::uword i {0}; i < n_lvls; i++) {
        (*lvls)(i) = levels[i];
    }
    
    delete velocities;
    velocities = new arma::colvec(n_lvls);
    for (arma::uword i {0}; i < n_lvls; i++) {
        (*velocities)(i) = velocities_arr[i];
    }
}

arma::colvec* MaterialModel::get_levels() {
    return lvls;
}

arma::colvec* MaterialModel::get_velocities() {
    return velocities;
}

arma::colvec* MaterialModel::get_control_function() {
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
        for (arma::uword i {0}; i < nn; i++) {
            std::cout << (*control_func)(i) << std::endl;
        }
        std::cout << "------------------" << std::endl;
    }
}

void MaterialModel::log_levels() {
    std::cout << "--------------------------------------------------------------------" << std::endl;
    for (arma::uword i {0}; i < n_lvls; i++) {
        std::cout << "    Level " << i+1 << ": " << (*lvls)(i);
    }
    std::cout << "\n--------------------------------------------------------------------" << std::endl;
}

void MaterialModel::log_velocities() {
    std::cout << "--------------------------------------------------------------------" << std::endl;
    for (arma::uword i {0}; i < n_lvls; i++) {
        std::cout << "    Vel " << i+1 << ": " << (*velocities)(i);
    }
    std::cout << "\n--------------------------------------------------------------------" << std::endl;
}
