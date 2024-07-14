#ifndef TENSOR_LATTICE_H
#define TENSOR_LATTICE_H

#include <Eigen/Dense>
#include <vector>
#include <iostream>

#include "vec3.h"

class tensor_lattice
{
public:
    std::vector< std::vector< std::vector< vec3 > > > data;
    int size;

    tensor_lattice(int N, vec3 b1, vec3 b2, vec3 b3)
    {
        size = N;
        data = std::vector< std::vector< std::vector< vec3 > > >(N, std::vector< std::vector< vec3 > >(N, std::vector< vec3 >(N)));

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    data[i][j][k] = b1 * check_index(i, N) + b2 * check_index(j, N) + b3 * check_index(k, N);
                }
            }
        }
    }

private:
    double check_index(int index, int N)
    {
        if (index+1 > (N+1)/2) 
        {
            return index - N;
        } 
        else 
        {
            return index;
        }
    }
};

#endif