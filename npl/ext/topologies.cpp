#include <iostream>
#include "topologies.h"


using namespace std;

int Topologies::get_n_atoms() {
  return n_atoms;
}

void Topologies::create_row(int connectivity[], int occupancy_a[], int occupancy_b[], int index)
{
  cout << index;
  for (int i = 0; i < get_n_atoms(); i++) {
      if (occupancy_a[index] == 1) {
        bonds[2] += connectivity[i] * occupancy_b[i];
        bonds[0] += connectivity[i] * occupancy_a[i];
      } else {
        bonds[2] += connectivity[i] * occupancy_a[i];
        bonds[1] += connectivity[i] * occupancy_b[i];

      }
    }
    
}

int * Topologies::get_feature_vector() {
  return bonds;
}
