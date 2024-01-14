#ifndef VERTEX_H
#define VERTEX_H
#include "Enums.h"
// Define the Vertex class template
template <typename IT>
class Vertex {
public:
    IT TreeField;
    IT BridgeField;
    IT ShoreField;
    IT AgeField;

    // Constructor
    Vertex(IT age,short int Label)
        : TreeField(-1),
          BridgeField(-1), 
          ShoreField(-1), 
          AgeField(age) {}

    // Copy constructor
    Vertex(const Vertex& other)
        : TreeField(other.TreeField), 
          BridgeField(other.BridgeField),
          ShoreField(other.ShoreField), 
          AgeField(other.AgeField) {}
          

    // Default constructor
    Vertex() : TreeField(-1),
          BridgeField(-1), 
          ShoreField(-1), 
          AgeField(-1) {}
    
    // Method to check if the vertex is reached
    bool IsReached() const {
        return AgeField != -1;
    }

    // Method to check if the vertex is reached
    bool IsEven() const {
        return AgeField % 2 == 0;
    }

    // Method to check if the vertex is reached
    bool IsOdd() const {
        return AgeField % 2;
    }

    // Utility function to print vertex information
    void print() const {
        std::cout 
                  << "TreeField: " << TreeField << ", "
                  << "BridgeField: " << BridgeField << ", "
                  << "ShoreField: " << ShoreField << ", "
                  << "AgeField: " << AgeField << ", "
                  << "IsReached: " << IsReached() << ", "
                  << std::endl;

    }
};
#endif //VERTEX_H
