#include "Mesh.h"
#include <iostream>

    
SquareLinearMesh::SquareLinearMesh(unsigned int number_elements_length, 
                                   unsigned int number_elements_depth,
                                   double length, 
                                   double depth) 
                                   :
                                   nel(number_elements_length), 
                                   ned(number_elements_depth), 
                                   element_length(length/number_elements_length), 
                                   element_depth(depth/number_elements_depth),
                                   domain_length(length), 
                                   domain_depth(depth),
                                   nn((number_elements_length+1)*(number_elements_depth+1)), 
                                   ne(number_elements_length*number_elements_depth){
    gen_coordinates();
    gen_connectivity();
}

SquareLinearMesh::~SquareLinearMesh(){
        delete coordinates;
        delete connectivity;
    }

void SquareLinearMesh::gen_coordinates() {
    coordinates = new Node[nn];
    for (unsigned int i {0}; i < nel+1; i++){  //X-coordinates
        for (unsigned int j {0}; j < ned+1; j++){
            coordinates[i+(nel+1)*j].x = element_length*i;
        }
    }
    for (unsigned int i {0}; i < ned+1; i++){  //Y-coordinates
        for (unsigned int j {0}; j < nel+1; j++){
            coordinates[j+(nel+1)*i].y = element_depth*i;
        }
    }
}     

void SquareLinearMesh::gen_connectivity() {
    connectivity = new Element[ne];    
    for (unsigned int n {0}; n < ned; n++){
        for (unsigned int i {0}; i < nel; i++){
            connectivity[i+nel*n].global_node1 = 1 + i + nel*n + n;
            connectivity[i+nel*n].global_node2 = 2 + i + nel*n + n;
            connectivity[i+nel*n].global_node3 = connectivity[i+nel*n].global_node2 + nel + 1;
            connectivity[i+nel*n].global_node4 = connectivity[i+nel*n].global_node1 + nel + 1;
        }
    }
}    
    
Node * SquareLinearMesh::get_coordinates() const{
    return coordinates;
}

Element * SquareLinearMesh::get_connectivity() const {
    return connectivity;
}

unsigned int SquareLinearMesh::get_nel() const {
    return nel;
}

unsigned int SquareLinearMesh::get_ned() const {
    return ned;
}

double SquareLinearMesh::get_domain_length() const {
    return domain_length;
}

double SquareLinearMesh::get_domain_depth() const {
    return domain_depth;
}

unsigned int SquareLinearMesh::get_number_of_nodes() const {
    return nn;
}

unsigned int SquareLinearMesh::get_number_of_elements() const {
    return ne;
}
    
void SquareLinearMesh::log_coordinates() const {
    if (nn < 100) {
        std::cout << "---------------------------------------------" << std::endl;
        for (unsigned int i {0} ; i < nn; i++){
            std::cout << "Node: " << i+1 << "    x:  " << coordinates[i].x << "    y:  " << coordinates[i].y << std::endl;
        }
        std::cout << "---------------------------------------------" << std::endl;
    }
}
    
void SquareLinearMesh::log_connectivity() const {
    if (ne < 100) {
        std::cout << "---------------------------------------------" << std::endl;
        for (unsigned int i {0}; i < ne; i++){
            std::cout << "Element: " << i+1 << "    node 1: " << connectivity[i].global_node1 << "    node 2: " << connectivity[i].global_node2 << "    node 3: " << connectivity[i].global_node3 
            << "    node 4: " << connectivity[i].global_node4 << std::endl; 
        }
        std::cout << "---------------------------------------------" << std::endl;
    }
}
    