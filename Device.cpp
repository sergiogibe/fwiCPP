#include "Device.h"
#include <iostream>
#include <cmath>

Device::Device(const SquareLinearMesh & mesh, 
               double x_pos, 
               double y_pos, 
               unsigned int dev_id) 
               : 
               x(x_pos), 
               y(y_pos), 
               id(dev_id) {
    
    //find X nodal position:
    double side_ratio {x / mesh.get_domain_length()};
    if (side_ratio > 1.0) { 
        side_ratio = 1.0; 
    }
    if (side_ratio < 0.0) { 
        side_ratio = 0.0; 
    }
    
    int nodeX {static_cast<int>(round((mesh.get_nel() + 1.0) * side_ratio))};
    if (nodeX == 0) { 
        nodeX = 1; 
    }
    
    //find Y nodal position:
    side_ratio = y / mesh.get_domain_depth();
    if (side_ratio > 1.0) { 
        side_ratio = 1.0; 
    }
    if (side_ratio < 0.0) { 
        side_ratio = 0.0; 
    }
    
    int nodeY {static_cast<int>(round((mesh.get_ned() + 1.0) * side_ratio))};
    if (nodeY == 0) { 
        nodeY = 1; 
    }
    
    global_node = nodeX + (mesh.get_nel() + 1) * (nodeY - 1);
} 

Device::~Device() {
}

unsigned int Device::get_global_node(){
    return global_node;
}

double Device::get_xcoord(){
    return x;
}

double Device::get_ycoord(){
    return y;
}

void Device::log(){
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Device: " << id << "    x: " << x << "    y: " << y << std::endl;
    std::cout << "Global node: " << global_node << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
}

