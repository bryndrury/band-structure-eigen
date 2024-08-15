#pragma once

#include "vec3.h"

#define pi 3.14159265358979323846
#define eV_to_Ryd 1 / 13.6
#define c 3.80998 // hbar**2 / (2*m) in Angstrom^2 * eV

#define d 0.350 // angstroms
#define rc 0.943 // angstroms
#define U0 -31.3 // angstroms
#define k_B 8.617333262145e-5 // eV / K
#define T 300 // K

double e_F = 11.6; // eV

// define the lattice vectors
double a = 4.05; // angstroms
vec3 a1( a/2, a/2, 0);
vec3 a2( a/2, 0, a/2);
vec3 a3( 0, a/2, a/2);

// define the reciprocal lattice vectors
vec3 b1 = (a2.cross(a3) * 2 * pi) / a1.dot(a2.cross(a3));
vec3 b2 = (a3.cross(a1) * 2 * pi) / a1.dot(a2.cross(a3));
vec3 b3 = (a1.cross(a2) * 2 * pi) / a1.dot(a2.cross(a3));

// Define the High symmetry points
double reciprocal_constant = (2 * pi / a);

// Define the High symmetry points with the multiplication applied directly
vec3 L          (0.5 * reciprocal_constant,  0.5 * reciprocal_constant,  0.5 * reciprocal_constant);
vec3 gamma      (0 * reciprocal_constant,    0 * reciprocal_constant,    0 * reciprocal_constant);
vec3 X          (0 * reciprocal_constant,    1 * reciprocal_constant,    0 * reciprocal_constant);
vec3 pointU     (0.25 * reciprocal_constant, 1 * reciprocal_constant,    0.25 * reciprocal_constant);

// The path 
std::vector< vec3 > paths = { L, gamma, X, pointU, gamma };
std::vector< std::string > path_labels = { "L", "Gamma", "X", "U", "Gamma" };