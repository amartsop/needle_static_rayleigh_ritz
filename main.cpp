#include <iostream>
#include <armadillo>
#include <vector>

#include "modes_magnitude.h"
#include "rayleigh_ritz_beam.h"
#include "post_processing_rr.hpp"

#include "gnuplot-iostream.h"
#include "numerical_integration.hpp"
#include "needle_animation.hpp"

int main(int argc, char *argv[])
{
    // Axial and bending dofs
    int nu = 4, nv = 4, nw = 4;

    // Spatial beam using Rayleigh-Ritz 
    RayleighRitzBeam needle(nu, nv, nw);

    // /********************* Simulation ************************/ 
    // Initial conditions
    arma::dvec q0 = arma::zeros<arma::dvec>(needle.get_model_size());
    arma::dvec q0_dot = arma::zeros<arma::dvec>(needle.get_model_size());

    // State vector initialization
    std::vector<arma::dvec> state_vector;
    arma::dvec state0 = arma::join_vert(q0, q0_dot);
    state_vector.push_back(state0);

    // Timing
    double t_final = 10.0; // Final time (s)
    double fs = 1e3;  // Simulation frequency (Hz)
    double h = 1.0 / fs; // Integration time step (s)
    double t = 0; // Initial time (s) 

    // Time vector initialization 
    std::vector<double> time_vector;
    time_vector.push_back(t);

    // Post processing 
    PostProcessingRR<RayleighRitzBeam> post_processing_rr(&needle);

    // Problem solver 
    NumericalIntegration<RayleighRitzBeam> ni(&needle, h, 4);

    // Iteration counter 
    uint counter = 0;

    // Deformation vector 
    std::vector<arma::dvec> r_vec;

    while (t <= t_final)
    {
        // System solution 
        arma::dvec x = ni.solve(state_vector.at(counter), t);
        state_vector.push_back(x);

        if (!x.is_finite()) {
            std::cout << "Error" << std::endl; 
            break; 
        };

        // Update time vector 
        time_vector.push_back(t);

        // Update time and counter
        t += h; counter += 1;

        std::cout << t << std::endl;
    }

    /******** Animation ********/
    NeedleAnimation needle_animation(&needle, &post_processing_rr);
    double animation_fs = 1e2; // Animation frequency
    int steps = fs / animation_fs;

    // Rigid body position and orientation
    arma::dvec roa_g_g = {0.0, 0.0, 0.0};
    arma::dvec euler_angles = {0.0, 0.0, 0.0};

    for(size_t i = 0; i <= state_vector.size(); i = i + steps)
    {
        // Current state and time 
        arma::dvec x_current = state_vector.at(i);
        double t_current = time_vector.at(i);

        // Elastic coordinates
        arma::dvec qf = x_current.rows(0, needle.get_model_size() - 1);

        // Animation
        needle_animation.animate(roa_g_g, euler_angles, qf);
    }
    
}
