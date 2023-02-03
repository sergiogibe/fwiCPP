#ifndef _DEVICE_H_
#define _DEVICE_H_
#include "Mesh.h"

class Device {
//Attributes
private:
    double x;
    double y;
    unsigned int id; 
    
    unsigned int global_node;

//Methods
public:

    //Constructors
    Device(const SquareLinearMesh & mesh_obj, double x_pos, double y_pos, unsigned int id);
    ~Device();
    
    //Getters
    unsigned int get_global_node();
    double get_xcoord();
    double get_ycoord();
    
    //Loggers
    void log();
};

#endif // _DEVICE_H_
