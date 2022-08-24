#include <vector>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <string>
#include <pybind11/stl.h>

class SomeClass {
   float multiplier; 


public:
   SomeClass(float multiplier_) : multiplier(multiplier_) {};
   float multiply(float input) {
      return multiplier * input;
   }

   std::vector<float> multiply_items(std::vector<float> items) {
      for (auto i = 0; i < items.size(); i++) {
         items[i] = multiply(items.at(i));
      }
      return items;
   }

   std::vector<std::vector<uint8_t>> make_image() {
      auto out = std::vector<std::vector<uint8_t>>();
      for (auto i = 0; i < 128; i++) {
         out.push_back(std::vector<uint8_t>(64));
      }
      for (auto i = 0; i < 30; i++) {
         for (auto j = 0; j < 30; j++) { out[i][j] = 255; }
      }
    return out;
  }

};