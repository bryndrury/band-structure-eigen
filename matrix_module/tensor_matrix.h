#ifndef TENSOR_MATRIX_H
#define TENSOR_MATRIX_H

#include <vector>
#include <iostream>

#include "tensor_lattice.h"
#include "value_definitions.h"

class tensor_matrix
{
public:
    std::vector< std::vector< std::vector< double > > > data;
    int size;

    tensor_matrix() {};

    tensor_matrix(tensor_lattice &K);

private:
    double calculate_potential(vec3 k);
};

tensor_matrix::tensor_matrix(tensor_lattice &K)
{
    size = K.size;
    data = std::vector< std::vector< std::vector< double > > >(size, std::vector< std::vector< double > >(size, std::vector< double >(size)));

    for (int i = 0; i < size; i++) 
    {
        for (int j = 0; j < size; j++) 
        {
            for (int k = 0; k < size; k++) 
            {
                data[i][j][k] = calculate_potential(K.data[i][j][k]);
            }
        }
    }
}

double tensor_matrix::calculate_potential(vec3 k)
{
    double Kn = k.norm();

    if (Kn > 0) 
    {
        return U0 * exp(-rc/d) * ( sin(rc*Kn) / (d*Kn * ( pow(d*Kn,2) + 1)) + cos(rc*Kn) / (pow(d*Kn,2) + 1) );
    } 
    else 
    {
        return U0 * exp(-rc/d) * (rc/d+1);
    }
}

#endif