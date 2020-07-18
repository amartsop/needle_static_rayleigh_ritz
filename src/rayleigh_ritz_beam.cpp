#include "rayleigh_ritz_beam.h"


RayleighRitzBeam::RayleighRitzBeam(uint axial_dofs, uint bending_y_dofs, uint bending_z_dofs)
{
    // Number of axial dofs 
    m_axial_dofs = axial_dofs;

    // Number of bending dofs y direction
    m_bending_y_dofs = bending_y_dofs;

    // Number of bending dofs z direction
    m_bending_z_dofs = bending_z_dofs;

    // Number of dofs
    m_dofs = m_axial_dofs + m_bending_y_dofs + m_bending_z_dofs;

    // Beam length 
    m_beam_length = m_needle_length;

    // Beam length 
    m_beam_radius = m_needle_radius;

    // Cross-sectional area
    m_beam_area = M_PI * pow(m_beam_radius, 2.0);

    // Beam mass (kg)    
    m_beam_mass = m_beam_density * m_beam_area * m_beam_length;

    // Beam area moment of inertia iyy
    m_iyy = (M_PI / 4) * pow(m_beam_radius, 4.0);

    // Beam area moment of inertia izz
    m_izz = (M_PI / 4) * pow(m_beam_radius, 4.0);

    // Beam young modulus
    m_beam_young_modulus = m_needle_young_modulus;

    // Beam density
    m_beam_density = m_needle_density;

    // Locator vectors
    m_lu = arma::regspace<arma::ivec>(1, 1, m_axial_dofs);
    m_lv = arma::regspace<arma::ivec>(m_axial_dofs + 1, 1,
        m_axial_dofs + m_bending_y_dofs);
    m_lw = arma::regspace<arma::ivec>(m_axial_dofs + m_bending_y_dofs + 1, 1,
        m_axial_dofs + m_bending_y_dofs + m_bending_z_dofs);

    // Locator matrices
    m_lu_mat = dm::locator_matrix(m_lu, m_dofs);
    m_lv_mat = dm::locator_matrix(m_lv, m_dofs);
    m_lw_mat = dm::locator_matrix(m_lw, m_dofs);


    // Initialization of frequencies vectors 
    m_u_freq = arma::zeros(m_axial_dofs);
    m_v_freq = arma::zeros(m_bending_y_dofs);
    m_w_freq = arma::zeros(m_bending_z_dofs);

    // Calculate mass matrix
    mass_matrix_calculation();
    
    // Calculate stiffness matrix 
    stiffness_matrix_calculation();
}


void RayleighRitzBeam::mass_matrix_calculation(void)
{
    arma::dmat mu = arma::eye<arma::dmat>(m_axial_dofs, m_axial_dofs);
    arma::dmat mv = arma::eye<arma::dmat>(m_bending_y_dofs, m_bending_y_dofs);
    arma::dmat mw = arma::eye<arma::dmat>(m_bending_z_dofs, m_bending_z_dofs);

    m_mass = m_lu_mat.t() * mu * m_lu_mat + m_lv_mat.t() * mv * m_lv_mat + 
        m_lw_mat.t() * mw * m_lw_mat;
}


void RayleighRitzBeam::stiffness_matrix_calculation(void)
{
    // Wave coefficient
    double c = sqrt(m_beam_young_modulus / m_beam_density);


    // Axial stiffness matrix 
    arma::dmat ku = arma::zeros<arma::dmat>(m_axial_dofs, m_axial_dofs);
    for (uint i = 0; i < m_axial_dofs; i++)
    {
        uint n = i + 1;
        double un_freq = (2.0 * n - 1.0) * (M_PI * c) / (2.0 * m_beam_length);
        ku(i, i) = pow(un_freq, 2.0);
    }

    // Bending y stiffness matrix 
    arma::dmat kv = arma::eye<arma::dmat>(m_bending_y_dofs, m_bending_y_dofs);
    for (uint i = 0; i < m_bending_y_dofs; i++)
    {
        uint n = i + 1; double sy_n = arma::as_scalar(m_sy(i));

        double bn = sy_n / m_beam_length;

        double vn_freq = sqrt((m_beam_young_modulus * m_izz) 
            / (m_beam_density * m_beam_area)) * pow((sy_n / m_beam_length), 2.0);

        double rv_n = ModesMagnitude::r(m_beam_length, bn);
        double rv_tilde_n = ModesMagnitude::r_tilde(m_beam_length, bn);

        kv(i, i) = (rv_n / rv_tilde_n) * pow(vn_freq, 2.0);
    }

    // Bending z stiffness matrix 
    arma::dmat kw = arma::eye<arma::dmat>(m_bending_z_dofs, m_bending_z_dofs);
    for (uint i = 0; i < m_bending_z_dofs; i++)
    {
        uint n = i + 1; double sz_n = arma::as_scalar(m_sz(i));

        double gn = sz_n / m_beam_length;

        double wn_freq = sqrt((m_beam_young_modulus * m_iyy) 
            / (m_beam_density * m_beam_area)) * pow((sz_n / m_beam_length), 2.0);

        double rw_n = ModesMagnitude::r(m_beam_length, gn);
        double rw_tilde_n = ModesMagnitude::r_tilde(m_beam_length, gn);

        kw(i, i) = (rw_n / rw_tilde_n) * pow(wn_freq, 2.0);
    }

    m_stiffness = m_lu_mat.t() * ku * m_lu_mat + m_lv_mat.t() * kv * m_lv_mat + 
        m_lw_mat.t() * kw * m_lw_mat;
}

arma::dmat RayleighRitzBeam::shape_function(double x)
{

    for (int i = 0)
    // /****** Axial direction ******/
    // // Magnitude
    // double un_hat = sqrt(2.0 / m_beam_mass);

    // // Frequency 

    // // Mode function 
    // double un_tilde = sin(un)

}