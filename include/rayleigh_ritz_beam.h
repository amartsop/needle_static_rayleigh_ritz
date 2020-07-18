#pragma once 

#include <iostream>
#include <armadillo>
#include "dynamics_math.h"
#include "needle_properties.hpp"
#include "modes_magnitude.h"


class RayleighRitzBeam : public NeedleProperties
{
public:
    RayleighRitzBeam(uint axial_dofs, uint bending_y_dofs, uint bending_z_dofs);

    // Mass matrix getter
    arma::dmat get_mass_matrix(void){ return m_mass; }
    
    // Stiffness matrix getter
    arma::dmat get_stiffness_matrix(void){ return m_stiffness; }

    // Update flexible body matrices and force vector
    void update(double t, arma::dvec q, arma::dvec q_dot);

    // Get system size 
    uint get_system_size(void) { return m_dofs; }

private:

    // Number of axial dofs 
    uint m_axial_dofs;

    // Number of bending dofs y direction
    uint m_bending_y_dofs;

    // Number of bending dofs z direction
    uint m_bending_z_dofs;

    // Number of dofs
    uint m_dofs;

private:
    // Beam lenght (m)
    double m_beam_length;

    // Beam radius (m)
    double m_beam_radius;

    // Beam cross-sectional area (m^2)
    double m_beam_area;

    // Beam mass (kg)    
    double m_beam_mass;
    
    // Beam area moment of inertia iyy (m^4)
    double m_iyy;

    // Beam area moment of inertia izz (m^4)
    double m_izz;

    // Beam young modulus (N / m^2)
    double m_beam_young_modulus;

    // Beam density (kg / m^2)
    double m_beam_density;

private:
    // Flexible body mass matrix
    arma::dmat m_mass;
    
    // Flexible body stiffness matrix
    arma::dmat m_stiffness;

private:
    // Mass matrix of flexible body calculation 
    void mass_matrix_calculation(void); 

    // Stiffness matrix of flexible body calculation 
    void stiffness_matrix_calculation(void); 

    // External force calculation
private:
    // Locator vectors
    arma::ivec m_lu, m_lv, m_lw;

    // Locator matrices
    arma::dmat m_lu_mat, m_lv_mat, m_lw_mat; 

    // Shape function 
    arma::dmat shape_function(void);

private:
    // Natural frequency y dircetion
    arma::dvec m_sy = {1.8751, 4.694, 7.8547, 10.9955};

    // Natural frequency z dircetion
    arma::dvec m_sz = {1.8751, 4.694, 7.8547, 10.9955};

    // Axial frequencies (rad / sec)
    arma::dvec m_u_freq;

    // Bending y frequencies (rad / sec)
    arma::dvec m_v_freq;

    // Bending z frequencies (rad / sec)
    arma::dvec m_w_freq;

private:
    // External force body frame (position l)
    arma::dvec external_force(double t, arma::dvec q, arma::dvec q_dot);

    // External traction body frame
    arma::dvec external_traction(double t, arma::dvec q, arma::dvec q_dot);
};


