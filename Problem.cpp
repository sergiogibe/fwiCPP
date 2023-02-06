#include "Problem.h"
#include "Device.h"
#include "Mesh.h"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <armadillo>
#include <filesystem>
#include "omp.h"


Problem::Problem(unsigned int nel,
                 unsigned int ned,
                 double d_length,
                 double d_depth,
                 double pulse_int = 1.0,
                 double pulse_freq = 2.0,
                 double obs_period = 2.6,
                 double deltaT = 0.002,
                 double control_init = -0.01) 
                 :     
                 number_elements_length(nel),
                 number_elements_depth(ned),
                 domain_length(d_length),
                 domain_depth(d_depth),
                 pulse_intensity(pulse_int),
                 pulse_frequency(pulse_freq),
                 delta_time(deltaT),
                 n_steps(static_cast<int>(floor(obs_period/deltaT)) + 1),
                 control_init_value(control_init) {
                   
    //Generating the problem's mesh
    problem_mesh = new SquareLinearMesh(number_elements_length, number_elements_depth, domain_length, domain_depth);
    
    //Generating control function
    control = new MaterialModel(problem_mesh->get_number_of_nodes(), control_init_value);
    
    //Initializing other members as nullptrs
    ricker_pulse = nullptr;
    global_stiffness_consistent = nullptr;
    global_mass_consistent = nullptr;
    global_mass = nullptr;
    solution = nullptr;
}

Problem::~Problem() {
    delete problem_mesh;
    delete control;
    delete ricker_pulse;
    delete global_stiffness_consistent;
    delete global_mass_consistent;
    delete global_mass;
    delete solution;
}


/*   ===============  Solutions  ===============  */

void Problem::solve(std::string mode) {
    
    //Always update the problem first (update the mass)
    update();
    
    //Number of shots
    unsigned int n_shots = sources.size();
    
    //Simplifying variables' names for cleaning syntax
    unsigned int nn = problem_mesh->get_number_of_nodes();
    double dt = delta_time;
    
    //Solution allocation
    solution = new arma::cube;
    solution->zeros(nn,n_steps,n_shots);
    
    arma::colvec p = arma::colvec(nn);
    arma::colvec dotp = arma::colvec(nn);
    arma::colvec ddotp = arma::colvec(nn);
    arma::colvec paux = arma::colvec(nn);
    arma::colvec bar_stiffness;
    double force {0.0};

    // Code below is temporary while we can't validate if the sparse matrix yields the exact same results
    // as the dense matrix. We can remove it after we add automated testing to see that the results are
    // indeed the same.
#define USE_SPARSE_MATRIX
#ifdef USE_SPARSE_MATRIX
    auto stiff_mat = arma::sp_mat(*global_stiffness_consistent);
#else
    auto stiff_mat = arma::mat(*global_stiffness_consistent);
#endif
#undef USE_SPARSE_MATRIX

    //Shots loop (PARALLEL OPENMP DIRECTIVE)
    //#pragma omp parallel for
    for (arma::uword s = 0; s < n_shots; s++) {
        std::cout << "Solving shot: " << s+1 << "   ...   " << std::flush;
        
        //Auxiliary p, dotp, ddotp vectors
        p.zeros(nn);
        dotp.zeros(nn);
        ddotp.zeros(nn);
        paux.zeros(nn);
        
        //Time loop
        for (arma::uword t = 1; t < n_steps; t++) {
                
            //Calculating term related to the constant stiffness
            bar_stiffness = stiff_mat*(p + dt * dotp + 0.5 * dt * dt * ddotp);

            //Node loop (PARALLEL OPENMP DIRECTIVE)
            //#pragma omp parallel for
            for (arma::uword n = 0; n < nn; n++) {
                    
                //Here we look whether the node "n" has a force or not (if it is source node or not)
                force = 0.0;
                if (n+1 == source_nodes(s)) {
                    force = (*ricker_pulse)(t);
                }
                    
                //Calculating ddotp
                ddotp(n) = (force - bar_stiffness(n)) / (*global_mass)(n);
                
                //Updating p(solution), dotp, and auxiliary vector "paux"
                p(n) += dt*dotp(n) + dt*dt*0.5*paux(n);
                dotp(n) += dt*0.5*(paux(n)+ddotp(n));
                paux(n) = ddotp(n);
                
                (*solution)(n,t,s) = p(n);
        
            }
        }
        std::cout << "Done. " << std::endl;
    }
    
}

void Problem::update() {
    
    //Update the problem calculating the mass matrix (and lumped one)
    calculate_mass_matrix();
    lumping_mass("diagonal_scale");
}

void Problem::build() {
    
    //Time setup
    generate_ricker_pulse();
    
    //Stiffness frame calculation
    calculate_stiffness_matrix();
    arma::sp_mat aux_sparse(*global_stiffness_consistent);
    global_stiffness_sparse = aux_sparse;
    
    //Setting up consistent and lumped global mass matrix
    global_mass_consistent = new arma::mat;
    global_mass_consistent->zeros(problem_mesh->get_number_of_nodes(),problem_mesh->get_number_of_nodes());
    global_mass = new arma::mat;
    global_mass->zeros(problem_mesh->get_number_of_nodes());
    
    //Source and receiver nodes (arrays that go into the solver)
    receiver_nodes.zeros(receivers.size());
    source_nodes.zeros(sources.size());
    int i {0};
    for (auto device: receivers) {
        receiver_nodes(i) = device.get_global_node();
        i++;
    }
    i = 0;
    for (auto device: sources) {
        source_nodes(i) = device.get_global_node();
        i++;
    }
}


/*   ===============  FEM methods   ===============  */

void Problem::lumping_mass(std::string mode = "row_sum") {
    
    //Lumping mass (Diagonal Scaling Method)
    if (mode == "diagonal_scale") {
        //Getting the sum of all elements of mass consistent and divide by the trace
        double total = arma::accu(*global_mass_consistent)/arma::trace(*global_mass_consistent);
        for (unsigned int i = 0; i < problem_mesh->get_number_of_nodes(); i++) {
            (*global_mass)(i) = total * ((*global_mass_consistent)(i,i));
        }
    }
    //Lumping mass (Row Sum Method)
    else if (mode == "row_sum") {
        (*global_mass) = arma::sum((*global_mass_consistent),1);
    }
}

void Problem::calculate_mass_matrix() {
    
    //Setting element mass matrix
    arma::mat element_mass;
    element_mass.zeros(4,4);
    
    //Summing all four gauss points contribution to element mass matrix
    //
    //Obs.: As the coordinates of nodes does not change, all element mass matrix
    //are equal. Therefore, we calculate it only considering the first element as a 
    //constant frame.
    //
    arma::mat element_coord;
    element_coord.zeros(4,2);
    unsigned int which_element {1}; //First element starts on 1.
    find_element_coordinates(element_coord, which_element);
    find_element_mass(element_mass, element_coord);
    
    //Looping all finite elements (assembling global matrix)
    arma::mat aux_matrix;
    aux_matrix.zeros(4);
    for (unsigned int e = 1; e <= problem_mesh->get_number_of_elements(); e++) {
        
        //Four nodes per element
        find_element_global_nodes(aux_matrix, e);
        
        //Finding material slowness property
        double mu = find_element_slowness(e);
        
        //Matrix assembling
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                (*global_mass_consistent)(aux_matrix(i)-1, aux_matrix(j)-1) += mu*element_mass(i,j);
            }
        }
    }
}

void Problem::calculate_stiffness_matrix() {
    
    //Setting element stiffness matrix
    arma::mat element_stiffness;
    element_stiffness.zeros(4,4);
    
    //Summing all four gauss points contribution to element stiffness matrix
    //
    //Obs.: As the coordinates of nodes does not change, all element stiffness matrix
    //are equal. Therefore, we calculate it only considering the first element as a 
    //constant frame.
    //
    arma::mat element_coord;
    element_coord.zeros(4,2);
    unsigned int which_element {1}; //First element starts on 1.
    find_element_coordinates(element_coord, which_element);
    find_element_stiffness(element_stiffness, element_coord);
    
    //Setting this->global stiffness matrix (arma matrix instance pointer to heap)
    global_stiffness_consistent = new arma::mat;
    global_stiffness_consistent->zeros(problem_mesh->get_number_of_nodes(),problem_mesh->get_number_of_nodes());
    arma::mat aux_matrix;
    aux_matrix.zeros(4);
    
    //Looping all finite elements (assembling global matrix)
    for (unsigned int e = 1; e <= problem_mesh->get_number_of_elements(); e++) {
        
        //Four nodes per element
        find_element_global_nodes(aux_matrix, e);
        
        //Matrix assembling
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                (*global_stiffness_consistent)(aux_matrix(i)-1, aux_matrix(j)-1) += element_stiffness(i,j);
            }
        }
    }
}

integration_output Problem::integrate(double r, double s, arma::mat element_coord, double& det_jacob) {
    
    //Finding shape functions nodal value (not always used, unless finding MASS)
    arma::mat shape_mat = 
    {
        0.25 * (1 - r) * (1 - s),
        0.25 * (1 + r) * (1 - s),
        0.25 * (1 + r) * (1 + s),
        0.25 * (1 - r) * (1 + s)
    };
    
    //Finding derivatives in respect to r and s
    arma::mat derivatives = 
    {
        {-0.25 * (1 - s), 0.25 * (1 - s), 0.25 * (1 + s), -0.25 * (1 + s)},
        {-0.25 * (1 - r), -0.25 * (1 + r), 0.25 * (1 + r), 0.25 * (1 - r)}
    };
    
    //Finding the jacobian matrix for the element
    arma::mat jacob = derivatives*element_coord;
    det_jacob = arma::det(jacob);
    
    //Calculating the deformation matrix for the element
    arma::mat deformation_mat = arma::solve(jacob, derivatives);
    
    //Return
    integration_output out{shape_mat, deformation_mat};

    return out;
}

arma::mat Problem::load_gauss_points() {
    
    //Gauss quadrature constants
    const double value {0.57735026905};
    const double weight {1.00000};
    
    //Generating auxiliary gauss points matrix
    arma::mat pts =
    {
        { value,   value,  weight},
        {-value,   value,  weight},
        { value,  -value,  weight},
        {-value,  -value,  weight}
    };
    
    return pts;
}

void Problem::find_element_coordinates(arma::mat& element_coord, unsigned int element) {
    
    //Element X-coordinates
    element_coord(0,0) = problem_mesh->get_coordinates()[problem_mesh->get_connectivity()[element-1].global_node1-1].x;
    element_coord(1,0) = problem_mesh->get_coordinates()[problem_mesh->get_connectivity()[element-1].global_node2-1].x;
    element_coord(2,0) = problem_mesh->get_coordinates()[problem_mesh->get_connectivity()[element-1].global_node3-1].x;
    element_coord(3,0) = problem_mesh->get_coordinates()[problem_mesh->get_connectivity()[element-1].global_node4-1].x;
    
    //Element Y-coordinates
    element_coord(0,1) = problem_mesh->get_coordinates()[problem_mesh->get_connectivity()[element-1].global_node1-1].y;
    element_coord(1,1) = problem_mesh->get_coordinates()[problem_mesh->get_connectivity()[element-1].global_node2-1].y;
    element_coord(2,1) = problem_mesh->get_coordinates()[problem_mesh->get_connectivity()[element-1].global_node3-1].y;
    element_coord(3,1) = problem_mesh->get_coordinates()[problem_mesh->get_connectivity()[element-1].global_node4-1].y;
}

void Problem::find_element_global_nodes(arma::mat& aux_matrix, unsigned int element) {
    
    //Four nodes per element
    aux_matrix(0) = problem_mesh->get_connectivity()[element-1].global_node1;
    aux_matrix(1) = problem_mesh->get_connectivity()[element-1].global_node2;
    aux_matrix(2) = problem_mesh->get_connectivity()[element-1].global_node3;
    aux_matrix(3) = problem_mesh->get_connectivity()[element-1].global_node4;
}
    
void Problem::find_element_stiffness(arma::mat& element_stiff, arma::mat& element_coord) {
    
    //Loading matrix containing all gauss points for spacial integration
    arma::mat pts = load_gauss_points();
    
    //Setting reference to the determinant of Jacobian
    double det_jacob {0.0};
    
    //Looping all 4 Gauss points 
    for (int i = 0; i < 4; i++) {
        double r = pts(i,0);
        double s = pts(i,1);
        double weight = pts(i,2);
        
        //Finding deformation matrix "B"
        arma::mat def_matrix = integrate(r, s, element_coord, det_jacob).B;
        
        //Transpose ("BTB term")
        arma::mat deft_matrix = def_matrix.t();
        
        //Element stiffness composition
        element_stiff += (deft_matrix * def_matrix) * det_jacob * weight;
    }
}

void Problem::find_element_mass(arma::mat& element_mass, arma::mat& element_coord) {
    
    //Loading matrix containing all gauss points for spacial integration
    arma::mat pts = load_gauss_points();
    
    //Setting reference to the determinant of Jacobian
    double det_jacob {0.0};
    
    //Looping all 4 Gauss points 
    for (int i = 0; i < 4; i++) {
        double r = pts(i,0);
        double s = pts(i,1);
        double weight = pts(i,2);
        
        //Finding shape functions vector "N"
        arma::mat shape_matrix = integrate(r, s, element_coord, det_jacob).N;
        
        //Transpose ("NTN term")
        arma::mat shapet_matrix = shape_matrix.t();
        
        //Element stiffness composition
        element_mass += (shapet_matrix * shape_matrix) * det_jacob * weight;
    }
}

double Problem::find_element_slowness(unsigned int element) {
    
    //Finding the control value at each node for element "element"
    arma::mat aux_element_control;
    aux_element_control.zeros(4);
    aux_element_control(0) = (*control->get_control_function())(problem_mesh->get_connectivity()[element-1].global_node1-1);
    aux_element_control(1) = (*control->get_control_function())(problem_mesh->get_connectivity()[element-1].global_node2-1);
    aux_element_control(2) = (*control->get_control_function())(problem_mesh->get_connectivity()[element-1].global_node3-1);
    aux_element_control(3) = (*control->get_control_function())(problem_mesh->get_connectivity()[element-1].global_node4-1);
    
    //Finding the slowness for each velocity
    unsigned int n_vels = control->get_nvels();
    arma::mat slowness_values;
    slowness_values.zeros(n_vels);
    for (unsigned int i = 0; i < n_vels; i++) {
        slowness_values(i) = 1/(pow((*control->get_velocities())(i),2));
    }
    
    //Using the mean to calculate the control value "cv" for this element
    double cv {0.0};
    for (int i = 0; i < 4; i++) {
        cv += 0.25*aux_element_control(i);
    }

    //Finding in which region the element belongs and calculate the resultant "mu"
    //
    //Obs.: Using hat-shape function: 
    //mu_e = mu(i+1)*(x-l(i)/l(i+1)-l(i)) + mu(i)*(1- (x-l(i)/l(i+1)-l(i)))
    //
    
    //First we need to identify who l(i+1) and l(i) are. Upper and lower bounds, respectively.
    //
    arma::colvec levels = *control->get_levels();
    unsigned int n_lvls = control->get_nlvls();
    double lower_b {0.0};
    unsigned int lower_b_index {0};
    double upper_b {levels(1)};
    for (unsigned int i = 0; i < n_lvls-1; i++) {
        if (cv-levels(i+1) < 0) {
            upper_b = levels(i+1);
            lower_b = levels(i);
            lower_b_index = i;
            break;
        }
        else if (cv >= levels(n_lvls-1)) {
            return slowness_values(n_vels-1);
        }
    }
    
    //Now we can use the hat-shape function
    double aux_term = (cv-lower_b)/(upper_b-lower_b);
    double element_mu = slowness_values(lower_b_index+1)*(aux_term) + slowness_values(lower_b_index)*(1-aux_term);
    
    return element_mu;
}

/*   ===============  Control methods  ===============  */

void Problem::set_control(unsigned int n_levels, std::vector<double> levels, std::vector<double> velocities_arr) {
    control->set(n_levels, levels, velocities_arr);
}

void Problem::print_levels() {
    control->log_levels();
}

void Problem::print_velocities() {
    control->log_velocities();
}


/*  ===============   Ricker's Pulse methods  ===============  */

void Problem::generate_ricker_pulse() {
    
    //Dinamically allocated local scope auxiliary time array (needs to be released after)
    double* time_arr = new double[n_steps];
    time_arr[0] = 0.0;
    for (unsigned int i {0}; i < n_steps - 1; i++) {
        time_arr[i+1] = time_arr[i] + delta_time;
    }
    
    //Math constants
    const double pi {3.14159265359};
    const double e  {2.71828182845};
    
    //Dinamically allocated arma matrix pointer attribute "this->ricker_pulse" 
    ricker_pulse = new arma::mat;
    ricker_pulse->zeros(n_steps);
    for (unsigned int i {0}; i < n_steps; i++) {
        time_arr[i] -= 1 / pulse_frequency;
        (*ricker_pulse)(i) = pulse_intensity * ((1 - (2 * pow(pi, 2)) * (pow(pulse_frequency, 2)) * (pow(time_arr[i], 2))) * pow(e,((-1) * (pow(pi, 2)) * (pow(pulse_frequency, 2)) * (pow(time_arr[i], 2)))));
    }
    
    //Releasing space in memory
    delete [] time_arr;
}

void Problem::print_pulse() {
    ricker_pulse->print();
}


/*  ===============   Devices methods  ===============  */

void Problem::add_device(std::string mode, double posx, double posy, unsigned int id) {
    
    //Adding devices accordingly with mode: "source" or "receiver"
    Device device {*problem_mesh, posx, posy, id};
    if (mode == "receiver") {
        this->receivers.push_back(device);
    }
    else if (mode == "source") {
        this->sources.push_back(device);
    }
    else {
        std::cout << "LOG: Incorrect device mode while creating the problem, please fix it." << std::endl;
        exit(1);
    }
    
}

void Problem::print_receivers() {
    for (auto device: receivers) {
        device.log();
    }
}

void Problem::print_sources() {
    for (auto device: sources) {
        device.log();
    }
}


/*  ===============   Data Output  ===============  */

void Problem::write_output() {
    
    std::cout << "Writing output to './Output/data.csv'    ...   " << std::flush;
    std::filesystem::create_directory("./Output/");
    //Open the CSV file for writing
    std::ofstream file("./Output/data.csv");
    if (!file.is_open()) {
        std::cerr << "Error opening file to write output." << std::endl;
    }
    
    //Useful shortcuts and settings
    unsigned int nn = problem_mesh->get_number_of_nodes();
    unsigned int ne = problem_mesh->get_number_of_elements();
    int decimal_places = 14;
    float factor = pow(10, decimal_places);
    unsigned int i_max {0};
    
    //Write the configuration as the first line (ROW 0)
    file 
    << nn                     << ","
    << number_elements_length << "," 
    << number_elements_depth  << ","
    << domain_length          << ","
    << domain_depth           << ","
    << pulse_intensity        << ","
    << pulse_frequency        << ","
    << delta_time             << ","
    << n_steps                << ","
    << receivers.size()       << ","
    << sources.size()         << ","
    << control->get_nlvls()   << "\n"; 
    
    //Write the Ricker pulse to the file (ROW 1)
    for (unsigned int i = 0; i < n_steps; ++i) {
        file << (*ricker_pulse)(i);
        if (i < n_steps - 1) {
            file << ",";
        }
    }
    file << "\n";
    
    //Write the control function to the file (ROW 2)
    for (unsigned int i = 0; i < nn; ++i) {
        file << (*control->get_control_function())(i);
        if (i < nn - 1) {
            file << ",";
        }
    }
    file << "\n";
    
    //Write control function levels to the file (ROW 3)
    for (unsigned int i = 0; i < control->get_nlvls(); ++i) {
        file << (*control->get_levels())(i);
        if (i < control->get_nlvls() - 1) {
            file << ",";
        }
    }
    file << "\n";
    
    //Write the receivers x coordinates to the file (ROW 4)
    for (unsigned int i = 0; i < receivers.size(); ++i) {
        file << receivers[i].get_xcoord();
        if (i < receivers.size() - 1) {
            file << ",";
        }
    }
    file << "\n";
    
    //Write the receivers y coordinates to the file (ROW 5)
    for (unsigned int i = 0; i < receivers.size(); ++i) {
        file << receivers[i].get_ycoord();
        if (i < receivers.size() - 1) {
            file << ",";
        }
    }
    file << "\n";
    
    //Write the sources x coordinates to the file (ROW 6)
    for (unsigned int i = 0; i < sources.size(); ++i) {
        file << sources[i].get_xcoord();
        if (i < sources.size() - 1) {
            file << ",";
        }
    }
    file << "\n";
    
    //Write the sources y coordinates to the file (ROW 7)
    for (unsigned int i = 0; i < sources.size(); ++i) {
        file << sources[i].get_ycoord();
        if (i < sources.size() - 1) {
            file << ",";
        }
    }
    file << "\n";
    
    //Write the connectivity to the file (ROW 8)
    for (unsigned int e = 0; e < ne; ++e) {
       file << problem_mesh->get_connectivity()[e].global_node1 << ",";
        file << problem_mesh->get_connectivity()[e].global_node2 << ",";
        file << problem_mesh->get_connectivity()[e].global_node3 << ",";
        file << problem_mesh->get_connectivity()[e].global_node4;
        i_max++;
        if (i_max != problem_mesh->get_number_of_elements()*4) {
            file << ",";
        }
    }
    file << "\n";
    
    //Write the solution to the file (ROW 9)
    i_max = 0;
    for (unsigned int s = 0; s < sources.size(); ++s) {
        for (unsigned int t = 0; t < n_steps; ++t) {
            for (unsigned int n = 0; n < nn; ++n) {
                file << round((*solution)(n,t,s) * factor) / factor;
                i_max++;
                if (i_max != sources.size()*n_steps*nn) {
                    file << ",";
                }
            }
        }
    }
    file << "\n";
    
    //Write the stiffness matrix (constant matrix) to the file (ROW 10)
    i_max = 0;
    for (arma::uword j = 0; j < nn; ++j) {
        for (arma::uword i = 0; i < nn; ++i) {
            file << round((*global_stiffness_consistent)(i,j) * factor) / factor;
            i_max++;
            if (i_max != nn*nn) {
                file << ",";
            }
        }
    }
    file << "\n";
    
    //Write the mass matrix to the file (ROW 10)
    i_max = 0;
    for (arma::uword j = 0; j < nn; ++j) {
        for (arma::uword i = 0; i < nn; ++i) {
            file << round((*global_mass_consistent)(i,j) * factor) / factor;
            i_max++;
            if (i_max != nn*nn) {
                file << ",";
            }
        }
    }
    file << "\n";

    //Close the file
    file.close();
    std::cout << "Done. " << std::endl;
}