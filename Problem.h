#ifndef _PROBLEM_H_
#define _PROBLEM_H_
#include "Mesh.h"
#include "MaterialModel.h"
#include "Device.h"
#include <vector>
#include <string>
#include <armadillo>


struct integration_output {
    arma::mat N;
    arma::mat B;
};


class Problem {

/*   ===============   Attributes    ===============   */
private:

    //Init list attributes
    const unsigned int number_elements_length;
    const unsigned int number_elements_depth;
    const double domain_length;
    const double domain_depth;
    const double pulse_intensity;
    const double pulse_frequency;
    const double delta_time;
    const unsigned int n_steps;
    const double control_init_value;
    
    //Key instance attributes
    SquareLinearMesh* problem_mesh;
    MaterialModel* control;
    arma::mat* ricker_pulse; 
    arma::sp_mat* global_stiffness_consistent;
    arma::sp_mat* global_mass_consistent;
    arma::mat* global_mass;
    std::vector<Device> receivers;
    std::vector<Device> sources;
    arma::Mat<unsigned int> receiver_nodes;
    arma::Mat<unsigned int> source_nodes;
    arma::cube* solution;
    arma::mat pml_thickness;
    unsigned int pml_max_nlayers;
    double pml_saturation;
    std::map<unsigned int, double> map_pml;
    
    
/*    ===============    Public Methods    ===============    */    
public:

    //Constructors
    Problem(unsigned int nel,
            unsigned int ned,
            double d_length,
            double d_depth,
            double pulse_int,
            double pulse_freq,
            double obs_period,
            double deltaT,
            double control_init);
    ~Problem();
    
    //Solutions
    void solve(std::string mode);
    void build();
    void pml_build();
    
    //Writers
    void print_levels();
    void print_velocities();
    void print_receivers();
    void print_sources();
    void print_pulse();
    void write_output(bool save_valid, bool save_solution);
    
    //Getters
    
    //Setters
    void set_control(unsigned int n_levels, std::vector<double> levels, std::vector<double> velocities_arr);
    void add_device(std::string mode, double posx, double posy, unsigned int id);
    void pml_set(arma::mat thickness, double saturation);
    


/*    ===============    Private Methods    ===============    */   
private:

    void generate_ricker_pulse();
    void calculate_stiffness_matrix();
    void calculate_mass_matrix();
    void lumping_mass(std::string mode);
    void update();
    
    arma::mat load_gauss_points();
    integration_output integrate(double r, double s, arma::mat element_coord, double& det_jacob);
    void find_element_coordinates(arma::mat& element_coord, unsigned int element);
    void find_element_global_nodes(arma::mat& aux_matrix, unsigned int element);
    void find_element_stiffness(arma::mat& element_stiff, arma::mat& element_coord);
    void find_element_mass(arma::mat& element_mass, arma::mat& element_coord);
    double find_element_slowness(unsigned int element); 
    
    double pml_saturation_curve(unsigned int layer);

};

#endif // _PROBLEM_H_
