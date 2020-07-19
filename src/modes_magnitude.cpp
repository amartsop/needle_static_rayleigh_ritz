#include "modes_magnitude.h"
// Magnitude R
double ModesMagnitude::r(double l, double bn)
{
    double n_a1n = a1n(l, bn); double n_a2n = a2n(l, bn);
    double n_a3n = a3n(l, bn); double n_a4n = a4n(l, bn);
    double n_a5n = a5n(l, bn); double n_a6n = a6n(l, bn);
    double n_a7n = a7n(l, bn); double n_a8n = a8n(l, bn);
    double n_a9n = a9n(l, bn); double n_a10n = a10n(l, bn);

    double n_jn = jn(l, bn);
    
    return n_a1n - 2.0 * n_a2n + n_a3n + pow(n_jn, 2.0) * (n_a4n - 2.0 * n_a5n + 
        n_a6n) + 2.0 * n_jn * (n_a7n - n_a8n - n_a9n + n_a10n);
}

// Magnitude R_tilde
double ModesMagnitude::r_tilde(double l, double bn)
{
    double n_a1n = a1n(l, bn); double n_a2n = a2n(l, bn);
    double n_a3n = a3n(l, bn); double n_a4n = a4n(l, bn);
    double n_a5n = a5n(l, bn); double n_a6n = a6n(l, bn);
    double n_a7n = a7n(l, bn); double n_a8n = a8n(l, bn);
    double n_a9n = a9n(l, bn); double n_a10n = a10n(l, bn);
    double n_jn = jn(l, bn);
    
    return n_a1n + 2.0 * n_a2n + n_a3n + pow(n_jn, 2.0) * (n_a4n + 2.0 * n_a5n + 
        n_a6n) + 2.0 * n_jn * (n_a7n + n_a8n + n_a9n + n_a10n);
}


// Integral A1n
double ModesMagnitude::a1n(double l, double bn)
{
    return (l / 2.0) - sin(2.0 * bn * l) / (4.0 * bn);
}

// Integral A2n
double ModesMagnitude::a2n(double l, double bn)
{
    return (cosh(bn * l) * sin(bn * l) - sinh(bn * l) * cos(bn * l)) / (2.0 * bn);
}

// Integral A3n
double ModesMagnitude::a3n(double l, double bn)
{
    return sinh(2.0 * bn * l) / (4.0 * bn) - (l / 2.0);
}

// Integral A4n
double ModesMagnitude::a4n(double l, double bn)
{
    return (l / 2.0) + sin(2.0 * bn * l) / (4.0 * bn) ;
}

// Integral A5n
double ModesMagnitude::a5n(double l, double bn)
{
    return (sinh(bn * l) * cos(bn * l) + cosh(bn * l) * sin(bn * l)) / (2.0 * bn);
}

// Integral A6n
double ModesMagnitude::a6n(double l, double bn)
{
    return (l / 2.0) + sinh(2.0 * bn * l) / (4.0 * bn) ;
}

// Integral A7n
double ModesMagnitude::a7n(double l, double bn)
{
    return pow(sin(bn * l), 2.0) / (2.0 * bn);
}

// Integral A8n
double ModesMagnitude::a8n(double l, double bn)
{
    return (sinh(bn * l) * sin(bn * l) - cosh(bn * l) * cos(bn * l) + 1.0) / (2.0 * bn);
}

// Integral A9n
double ModesMagnitude::a9n(double l, double bn)
{
    return (cosh(bn * l) * cos(bn * l) + sinh(bn * l) * sin(bn * l) - 1.0) / (2.0 * bn);
}

// Integral A10n
double ModesMagnitude::a10n(double l, double bn)
{
    return pow(sinh(bn * l), 2.0) / (2.0 * bn);
}

// Coefficient 
double ModesMagnitude::jn(double l, double bn)
{
    return (cos(bn * l) + cosh(bn * l)) / (sin(bn * l) - sinh(bn * l));
}