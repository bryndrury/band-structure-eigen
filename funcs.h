#pragma once

#include <vector>
#include <cmath>

#include "vec3.h"
#include "value_definitions.h"

inline double check_index(int index, int N)
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

std::vector< std::vector< std::vector< vec3 > > > generate_lattice(int N, vec3 b1, vec3 b2, vec3 b3)
{
    std::vector< std::vector< std::vector< vec3 > > > K(N, std::vector< std::vector< vec3 > >(N, std::vector< vec3 >(N)));

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                K[i][j][k] = b1 * check_index(i, N) + b2 * check_index(j, N) + b3 * check_index(k, N);
            }
        }
    }

    return K;
}

inline double calculate_potential(vec3 k)
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

std::vector< std::vector< std::vector< double > > > generate_potential(const std::vector< std::vector< std::vector< vec3 > > > &K)
{
    std::vector< std::vector< std::vector< double > > > U_lattice(K.size(), std::vector< std::vector< double > >(K.size(), std::vector< double >(K.size())));

    for (int i = 0; i < K.size(); i++) 
    {
        for (int j = 0; j < K.size(); j++) 
        {
            for (int k = 0; k < K.size(); k++) 
            {
                U_lattice[i][j][k] = calculate_potential(K[i][j][k]);
            }
        }
    }

    return U_lattice;
}

std::vector< std::vector< double > > generate_H_block_template(const std::vector< std::vector< std::vector< vec3 > > > &K, const std::vector< std::vector< std::vector< double > > > &U_K, vec3 k_vec)
{
    int N = K.size();
    
    // The empty H_block
    std::vector< std::vector< double > > H_block(N*N*N, std::vector< double >(N*N*N));

    // The empty U_K_flattened
    std::vector< double > U_K_flattened(N*N*N);

    for (int i = 0; i < N; i++) 
    {
        for (int j = 0; j < N; j++) 
        {
            for (int k = 0; k < N; k++) 
            {
                // Fill the U_K_flattened
                U_K_flattened[i + j*N + k*N*N] = U_K[i][j][k];

                // Fill the H_block
                int index_xy = i + j*N + k*N*N;
                H_block[index_xy][index_xy] = pow((K[i][j][k] + k_vec).norm(),2) * c;
            }
        }
    }

    for (int i = 0; i < N*N*N; i++) 
    {
        for (int j = 0; j < N*N*N; j++) 
        {
            H_block[i][j] += U_K_flattened[ abs(j - i) ];
        }
    }

    return H_block;
}



void update_H_block(Eigen::MatrixXd &H_block, const std::vector< std::vector< std::vector< vec3 > > > &K, const double U_K000, vec3 k_vec)
{
    int N = K.size();

    for (int i = 0; i < N; i++) 
    {
        for (int j = 0; j < N; j++) 
        {
            for (int k = 0; k < N; k++) 
            {
                H_block(i + j*N + k*N*N, i + j*N + k*N*N) = U_K000 + pow((K[i][j][k] + k_vec).norm(),2) * c;
            }
        }
    }
}

