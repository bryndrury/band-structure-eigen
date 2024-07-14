#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
// #include <omp.h>

#include <Eigen/Dense>

#include "matrix.h"

void write_eigenvalues_to_file(std::ofstream &file, const H_matrix &H, const int bands, const int cumulative_kstep);

int main(int argc, char *argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now();

    // Key parameters (take cli if valid)
    int N = (argc > 1) ? std::stoi(argv[1]) : 5;
    int ksteps = (argc > 2) ? std::stoi(argv[2]) : 25;
    int bands;
    if (argc > 3 && std::stoi(argv[3]) < N*N*N) {
        bands = std::stoi(argv[3]);
    } else {
        bands = 6;
    }

    // Generate the reciprocal lattice space
    tensor_lattice K( N, b1, b2, b3 );
    tensor_matrix U_K( K );
    H_matrix H( K, U_K, vec3() );

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
            file << cumulative_kstep << " ";

            write_eigenvalues_to_file(file, H, bands, cumulative_kstep);

            std::cout << "Step: " << step+1 << " of " << ksteps << " Path: " << path_index+1 << " of " << paths.size()-1 << "    \t\t\r" << std::flush;
        }
    }
}

void write_eigenvalues_to_file(std::ofstream &file, const H_matrix &H, const int bands, const int cumulative_kstep)
{
    for (int i = 0; i < bands; i++)
    {
        file << H.eigen_values[i] << " ";
    }
    file << "\n";
};