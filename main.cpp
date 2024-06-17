// Create the band structure of a 3D lattice of Aluminium using the tight binding model
// Created by: Bryn Drury
// Github: github.com/bryndrury/band-structure-cpp
// Written in C style C++

#include <iostream>
#include <fstream>

#include <vector>

#include <chrono>

#include <Eigen/Dense>

#include <omp.h>

#include "vec3.h"
#include "value_definitions.h"
#include "funcs.h"

// the problem is with the vector class and I think its in the generate lattice function

int main()
{
    omp_set_num_threads(6);

    // start the timer
    auto start = std::chrono::high_resolution_clock::now();

    // define the lattice vectors
    double a = 4.05; // angstroms
    vec3 a1( a/2, a/2, 0);
    vec3 a2( a/2, 0, a/2);
    vec3 a3( 0, a/2, a/2);

    // define the reciprocal lattice vectors
    // vec3 b1 = (a2.cross(a3) * 2 * pi) / a1.dot(a2.cross(a3));
    // b1 = 2 * np.pi * np.cross(a2, a3) / a1.dot( np.cross(a2, a3) )
    vec3 b1 = (a2.cross(a3) * 2 * pi) / a1.dot(a2.cross(a3));
    vec3 b2 = (a3.cross(a1) * 2 * pi) / a1.dot(a2.cross(a3));
    vec3 b3 = (a1.cross(a2) * 2 * pi) / a1.dot(a2.cross(a3));

    // High symmetry points
    double reciprocal_constant = (2 * pi / a);
    
    vec3 L          (0.5,  0.5,  0.5 );
    vec3 gamma      (0,    0,    0   );
    vec3 X          (0,    1,    0   );
    vec3 pointU     (0.25, 1,    0.25);

    L       *= reciprocal_constant;
    gamma   *= reciprocal_constant;
    X       *= reciprocal_constant;
    pointU  *= reciprocal_constant;

    // vec3 L          (0.5 * reciprocal_constant, 0.5 * reciprocal_constant, 0.5 * reciprocal_constant);
    // vec3 gamma      (0   * reciprocal_constant, 0   * reciprocal_constant, 0   * reciprocal_constant);
    // vec3 X          (0   * reciprocal_constant, 1   * reciprocal_constant, 0   * reciprocal_constant);
    // vec3 pointU     (0.25* reciprocal_constant, 1   * reciprocal_constant, 0.25* reciprocal_constant);

    // The path 
    std::vector< vec3 > paths = { L, gamma, X, pointU, gamma };
    std::vector< std::string > path_labels = { "L", "Gamma", "X", "U", "Gamma" };

    // Define the size of the central equation and the number of steps
    int N = 5;
    int ksteps = 25;
    int bands = 6;

    // label locations
    std::vector< int > label_locations(paths.size());

    // generate the H_block
    std::vector< std::vector< std::vector< vec3 > > > K = generate_lattice(N, b1, b2, b3);
    std::vector< std::vector< std::vector< double > > > U_K = generate_potential(K);
    std::vector< std::vector< double > > H_block_template = generate_H_block_template(K, U_K, vec3());

    // open the file to write the eigenvalues
    std::ofstream file;
    file.open("../energies.data");
    double cum_kstep;

    for (int path_index = 0; path_index < paths.size()-1; path_index++)
    {
        // write the label
        label_locations.push_back(path_index * ksteps);

        // calculate the kstep
        vec3 kstep = ( paths[path_index+1]-paths[path_index] ) / ksteps;

        // 'integrate' along the path
        // #pragma omp parallel for
        for (int step = 0; step < ksteps; step++)
        {
            // Create the Eigen matrices
            Eigen::MatrixXd H_block_eigen(H_block_template.size(), H_block_template[0].size());
            Eigen::MatrixXd eigen_vectors_eigen(H_block_template.size(), H_block_template[0].size());
            Eigen::VectorXd eigen_values_eigen(H_block_template.size());

            // Populate the H_block with the template
            for (int i = 0; i < H_block_template.size(); i++)
            {
                for (int j = 0; j < H_block_template[0].size(); j++)
                {
                    H_block_eigen(i, j) = H_block_template[i][j];
                }
            }

            // Calculate the q value

            vec3 qval = paths[path_index] + (kstep * step);

            // Update the H_block
            update_H_block(H_block_eigen, K, U_K[0][0][0], qval);

            // Solve the eigenvalue problem using Eigen
            Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > es(H_block_eigen);
            eigen_values_eigen = es.eigenvalues();

            // Store the eigenvalues in vector
            std::vector< double > eigenvalues_vector(H_block_template.size());
            for (int i = 0; i < H_block_template.size(); i++)
            {
                eigenvalues_vector[i] = eigen_values_eigen(i);
            }

            // Sort the eigenvalues by magnitude
            std::sort(eigenvalues_vector.begin(), eigenvalues_vector.end(), 
                [](double a, double b) { 
                    return std::abs(a) < std::abs(b); 
                }
            );

            // calculate the cumulative kstep (for the x-axis of the plot)
            cum_kstep += kstep.norm();

            // write the eigenvalues to the file
            #pragma omp critical
            {
                file << cum_kstep;
                for (int i = 0; i < bands; i++)
                {
                    file << " " << eigenvalues_vector[i];
                }
                file << std::endl;
            }

            std::cout << "Step: " << step << " of " << ksteps << " Path: " << path_index+1 << " of " << paths.size()-1 << "    \t\t\r" << std::flush;
        }
    }

    file.close();

    // end the timer
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "\n\n";
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microseconds" << std::endl;
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds" << std::endl;
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;
}