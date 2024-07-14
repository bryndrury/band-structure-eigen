#pragma once

#include <Eigen/Dense>
#include <vector>
#include <iostream>

#include "vec3.h"

struct tensor_lattice
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

struct tensor_matrix
{
public:
    std::vector< std::vector< std::vector< double > > > data;
    int size;

    tensor_matrix(tensor_lattice &K)
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

private:
    double calculate_potential(vec3 k)
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
};

struct H_matrix
{
    Eigen::MatrixXd eigen_data;
    int size;

    Eigen::VectorXd eigen_values;

    H_matrix(const tensor_lattice &K, const tensor_matrix &U_K, vec3 k_vec)
    {
        size = K.size;
        eigen_data = Eigen::MatrixXd(size*size*size, size*size*size);

        std::vector< double > U_K_flattened(size*size*size);

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                for (int k = 0; k < size; k++)
                {
                    // Fill the U_K_flattened
                    U_K_flattened[i + j*size + k*size*size] = U_K.data[i][j][k];

                    // Fill the H_block data
                    int index_xy = i + j*size + k*size*size;
                    eigen_data(index_xy, index_xy) = pow((K.data[i][j][k] + k_vec).norm(),2) * c;
                }
            }
        }

        for (int i = 0; i < size*size*size; i++)
        {
            for (int j = 0; j < size*size*size; j++)
            {
                eigen_data(i, j) += U_K_flattened[ abs(j - i) ];
            }
        }
    };

    // Update the H_block for new q value (just the diagonal)
    void update_H_block(const tensor_lattice &K, const double U_K_0, vec3 k_vec)
    {
        int N = K.size;

        for (int i = 0; i < N; i ++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    int index = i + j*N + k*N*N;
                    eigen_data(index, index) = U_K_0 + pow((K.data[i][j][k] + k_vec).norm(),2) * c;
                }
            }
        }
    };

    void calculate_eigenvalues()
    {
        int N = size * size * size;

        Eigen::MatrixXd eigen_vectors(N, N);

        Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > es(eigen_data);
        eigen_values = es.eigenvalues();

        // Sort the eigenvalues by magnitude
        std::vector< double > eigenvalues_vector(N);
        for (int i = 0; i < N; i++)
        {
            eigenvalues_vector[i] = eigen_values(i);
        }

        std::sort(eigenvalues_vector.begin(), eigenvalues_vector.end(), 
            [](double a, double b) { 
                return std::abs(a) < std::abs(b); 
            }
        );

        // Return the sorted eigenvalues to the eigen_values
        for (int i = 0; i < N; i++)
        {
            eigen_values(i) = eigenvalues_vector[i];
        }
    };
};