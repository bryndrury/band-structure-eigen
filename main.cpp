#include <iostream>
#include <fstream>

#include <vector>

#include <chrono>

#include <Eigen/Dense>

// #include <omp.h>

#include "vec3.h"
#include "value_definitions.h"
#include "matrix.h"

int main(int argc, char *argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now();

    // Key parameters (take cli if valid)
    int N = (argc > 1) ? std::stoi(argv[1]) : 5;
    int ksteps = (argc > 2) ? std::stoi(argv[2]) : 100;
    int bands;
    if (argc > 3 && std::stoi(argv[3]) < N*N*N) {
        bands = std::stoi(argv[3]);
    } else {
        bands = 6;
    }

    double e_F = 11.6; // eV

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

    // Generate the reciprocal lattice space
    tensor_lattice K( N, b1, b2, b3 );
    tensor_matrix U_K( K );
    H_matrix H( K, U_K, vec3() );


    // Define the High symmetry points
    double reciprocal_constant = (2 * pi / a);
    
    vec3 L          (0.5,  0.5,  0.5 );
    vec3 gamma      (0,    0,    0   );
    vec3 X          (0,    1,    0   );
    vec3 pointU     (0.25, 1,    0.25);

    L       *= reciprocal_constant;
    gamma   *= reciprocal_constant;
    X       *= reciprocal_constant;
    pointU  *= reciprocal_constant;

    // The path 
    std::vector< vec3 > paths = { L, gamma, X, pointU, gamma };
    std::vector< std::string > path_labels = { "L", "Gamma", "X", "U", "Gamma" };

    double U_K_0 = U_K.data[0][0][0];

    // Open the file to write the eigenvalues
    double cumulative_kstep = 0;
    std::ofstream file;
    file.open("../eigenvalues.txt");

    for (int path_index = 0; path_index < paths.size()-1; path_index++)
    {
        vec3 kstep = ( paths[path_index+1]-paths[path_index] ) / ksteps;

        for (int step = 0; step < ksteps; step++)
        {
            vec3 qval = paths[path_index] + (kstep * step);

            H.update_H_block(K, U_K_0, qval);

            H.calculate_eigenvalues();

            cumulative_kstep += kstep.norm();

            // write the eigenvalues to file
            file << cumulative_kstep << " ";
            bool surface = false;

            for (int i = 0; i < bands; i++)
            {
                file << H.eigen_values[i] << " ";
                if (H.eigen_values[i] > ( e_F - (k_B * T) ) && H.eigen_values[i] < ( e_F + (k_B * T) ))
                {
                    surface = true;
                    std::cout << "Surface at " << cumulative_kstep << "    \n";
                }
            }
            file << "\n";

            std::cout << "Step: " << step+1 << " of " << ksteps << " Path: " << path_index+1 << " of " << paths.size()-1 << "    \t\t\r" << std::flush;
        }
    }
}
