#ifndef _MESH_H_
#define _MESH_H_

    

struct Node {
    double x;
    double y;
};

struct Element {
    int global_node1;
    int global_node2;
    int global_node3;
    int global_node4;
};
    
    
class SquareLinearMesh {

//Attributes
private: 
    unsigned int nel;
    unsigned int ned;
    double element_length;        
    double element_depth;         
    double domain_length;
    double domain_depth;
    unsigned int nn;                 //number of nodes
    unsigned int ne;                 //number of elements

    Node* coordinates;
    Element* connectivity;


//Methods
public:

    //Constructors
    SquareLinearMesh(unsigned int number_elements_length, 
                     unsigned int number_elements_depth, 
                     double length, 
                     double depth); 
    ~SquareLinearMesh();
    
    
    //Getters
    Node* get_coordinates() const;
    Element* get_connectivity() const;
    unsigned int get_nel() const;
    unsigned int get_ned() const;
    double get_domain_length() const;
    double get_domain_depth() const;
    unsigned int get_number_of_nodes() const;
    unsigned int get_number_of_elements() const;
    double get_element_length() const;
    double get_element_depth() const;
    
    //Loggers
    void log_coordinates() const;
    void log_connectivity() const;
    
private:
    void gen_coordinates();
    void gen_connectivity();
        
};
                                                     




#endif