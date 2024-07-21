#ifndef MATRIX_H
#define MATRIX_H

#include <Eigen/Dense>
#include <vector>
#include <iostream>

#include "value_definitions.h"
#include "tensor_matrix.h"

class H_matrix
{
public:
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

#endif