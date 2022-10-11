#ifndef TOPOLOGIES_H
#define TOPOLOGIES_H

#include <iostream>
#include <string>

using namespace Eigen;

class Topologies{
    public:
        

        void create_row(int connectivity[], int occupancy_a[], int occupancy_b[], int index);
        int get_n_atoms();
        int * get_feature_vector();
        Topologies(int a = 0) {
            n_atoms = a;
        };
    
    private:
        int n_atoms;
        int bonds[3] = {0,0,0};
        

};

#endif  