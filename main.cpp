#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <omp.h>

#include <cstdlib>

#include <Eigen/Dense>

#include <matrix.h>

void write_eigenvalues_to_file(std::ofstream &file, const H_matrix &H, const int bands, const double cumulative_kstep);
double calculate_cumulative_kstep(const std::vector< vec3 > &paths, const int path_index, const int ksteps, const int step);

int main(int argc, char *argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now();

    // Key parameters (take cli if valid)
    int N = (argc > 1) ? std::stoi(argv[1]) : 12;
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
    H_matrix H_template( K, U_K, vec3() );

    double U_K_0 = U_K.data[0][0][0];

    // Open the file to write the eigenvalues
    std::ofstream file;
    std::string filename = "../eigenvalues.txt";
    // std::string filename = "../eigenvalues_" + std::to_string(N) + "_" + std::to_string(ksteps) + "_" + std::to_string(bands) + ".txt";
    file.open(filename);

    #pragma omp parallel for collapse(2)
    for (int path_index = 0; path_index < paths.size()-1; path_index++)
    {
        vec3 kstep = ( paths[path_index+1]-paths[path_index] ) / ksteps;

        for (int step = 0; step < ksteps; step++)
        {
            vec3 qval = paths[path_index] + (kstep * step);

            H_matrix H = H_template;
            H.update_H_block(K, U_K_0, qval);
            H.calculate_eigenvalues();

            write_eigenvalues_to_file(file, H, bands, calculate_cumulative_kstep(paths, path_index, ksteps, step));
        }
    }
    file.close();

    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start_time).count() << "ms\n";
}

void write_eigenvalues_to_file(std::ofstream &file, const H_matrix &H, const int bands, const double cumulative_kstep)
{
    #pragma omp critical
    {
        file << cumulative_kstep << " ";
        for (int i = 0; i < bands; i++)
        {
            file << H.eigen_values[i] << " ";
        }
        file << "\n";
    }
};

double calculate_cumulative_kstep(const std::vector< vec3 > &paths, const int path_index, const int ksteps, const int step)
{
    double cumulative_kstep = 0;

    for (int i = 0; i <= path_index; i++)
    {
        vec3 kstep = ( paths[i+1]-paths[i] ) / ksteps;

        if (i != path_index)
        {
            for (int j = 0; j <= ksteps; j++)
            {
                cumulative_kstep += kstep.norm();
            }
        }
        else
        {
            for (int j = 0; j <= step; j++)
            {
                cumulative_kstep += kstep.norm();
            }
        }
    }

    return cumulative_kstep;
}
