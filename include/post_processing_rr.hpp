#pragma once 

#include <iostream>
#include <armadillo>
#include <vector>

#include "gnuplot-iostream.h"
#include "dynamics_math.h"
#include "euler_rotations.h"

template <class T>
class PostProcessingRR
{
public:
    // Postprocessing for plotting and analysis (Rayleigh-Ritz)
    PostProcessingRR(T *needle);
    
    // Get the beam coordinates in the floating frame of reference
    arma::dmat get_beam_coordinates_ffr(arma::dmat elastic_coordinates);

    // Get the beam coordinates in the inertial frame of reference
    arma::dmat get_beam_coordinates_inertial(arma::dvec roc_g_g, 
        arma::dvec euler_angles, arma::dvec elastic_coordinates) ;

private:

    // Needle handle 
    T* m_needle_ptr;

    // Number of points in beam
    const int npel = 200;     

    // Beam length 
    double m_beam_length;
};

template <class T>
PostProcessingRR<T>::PostProcessingRR(T *needle)
{
    // Needle handle
    m_needle_ptr = needle;

    // Beam length
    m_beam_length = m_needle_ptr->get_beam_length();
}


template <class T>
arma::dmat PostProcessingRR<T>::get_beam_coordinates_ffr(arma::dmat elastic_coordinates)
{
    // Beam points 
    std::vector<double> m_beam_points_x, m_beam_points_y, m_beam_points_z;

    // Beam points armadillo 
    arma::dvec rx, ry, rz;
    
    // Discretize beam
    arma::dvec x_vec = arma::linspace<arma::dvec> (0, m_beam_length, npel);

    for (int i = 0; i < npel; i++) 
    { 
        // Dimensionless length
        double ksi = x_vec(i) / m_beam_length;

        // Position of point p0 wrt to f
        arma::dvec rap0_f_f = {x_vec(i), 0.0, 0.0};

        // Deformation of point p
        arma::dvec rp0p_f_f = m_needle_ptr->get_deflection(ksi, elastic_coordinates);

        // Final position of pj wrt to f origin
        arma::dvec rap_f_f = rap0_f_f + rp0p_f_f;

        m_beam_points_x.push_back(rap_f_f(0));
        m_beam_points_y.push_back(rap_f_f(1));
        m_beam_points_z.push_back(rap_f_f(2));
    }

    // Armadillo matrices of beam elements 
    rx = arma::vec(m_beam_points_x); ry = arma::vec(m_beam_points_y);
    rz = arma::vec(m_beam_points_z);

    return (arma::join_horiz(rx, ry, rz));
}


template <class T>
arma::dmat PostProcessingRR<T>::get_beam_coordinates_inertial(arma::dvec roa_g_g, 
    arma::dvec euler_angles, arma::dvec elastic_coordinates)
{
    // Rigid body coordinates
    arma::dmat rot_f_g = EulerRotations::rotation(euler_angles);
    
    // Beam points 
    std::vector<double> m_beam_points_x, m_beam_points_y, m_beam_points_z;

    // Beam points armadillo 
    arma::dvec rx, ry, rz;

    // Discretize beam
    arma::dvec x_vec = arma::linspace<arma::dvec> (0, m_beam_length, npel);

    for (int i = 0; i < npel; i++) 
    { 
        // Dimensionless length
        double ksi = x_vec(i) / m_beam_length;

        // Position of point p0 wrt to f
        arma::dvec rap0_f_f = {x_vec(i), 0.0, 0.0};

        // Deformation of point p
        arma::dvec rp0p_f_f = m_needle_ptr->get_deflection(ksi, elastic_coordinates);

        // Final position of pj wrt to f origin
        arma::dvec rap_f_f = rap0_f_f + rp0p_f_f;

        // Final position of pj wrt to G origin 
        arma::dvec rop_g_g = roa_g_g + rot_f_g *  rap_f_f;

        // Final position of pj wrt to f origin
        m_beam_points_x.push_back(rop_g_g(0));
        m_beam_points_y.push_back(rop_g_g(1));
        m_beam_points_z.push_back(rop_g_g(2));
    }

    // Armadillo matrices of beam elements 
    rx = arma::vec(m_beam_points_x); ry = arma::vec(m_beam_points_y);
    rz = arma::vec(m_beam_points_z);

    return (arma::join_horiz(rx, ry, rz));
}

