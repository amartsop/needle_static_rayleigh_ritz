#pragma once 

#include <iostream>
#include <armadillo>


class ModesMagnitude
{
public:
    ModesMagnitude() {};

    // Magnitude R
    static double r(double l, double bn);

    // Magnitude R_tilde
    static double r_tilde(double l, double bn);

    // Coefficient
    double static jn(double l, double bn);

private:
    // Integrals
    double static a1n(double l, double bn);
    double static a2n(double l, double bn);
    double static a3n(double l, double bn);
    double static a4n(double l, double bn);
    double static a5n(double l, double bn);
    double static a6n(double l, double bn);
    double static a7n(double l, double bn);
    double static a8n(double l, double bn);
    double static a9n(double l, double bn);
    double static a10n(double l, double bn);

};



