#ifndef MANDEL_PRESSURE_EG_H_
#define MANDEL_PRESSURE_EG_H_
// Mandel pressure solution that has a second component valued 0, in order to be used with EG pressure solution
#include "DealiiHeader.h"
using namespace std;
using namespace dealii;

class MandelPressureEG : public Function<dim>
{
public:
    MandelPressureEG(double _t) : Function<dim>(2) { t = _t; };
    virtual void vector_value(const Point<dim> &p, Vector<double> &values) const;

private:
    double t;
    /*
    // Data from Bin Wang
    double a = 100; // computational domain is axb
    double b = 10;
    double k = 1e-13;
    double mu_f =1e-3;
    double F = 5.94e8;
    double E = 5.94e9;
    double nu = 0.2;
    double alpha = 1;
    // c_0 = inv_M;
    double c_0 = 1 / 1.65e10;
    */

    // Data from Phillips Phillips
    double a = 1; // computational domain is axb
    double b = 1;
    double k = 1e-2;
    double mu_f =1.0;
    double F = 2.0e3;
    double E = 1.0e4;
    double nu = 0.2;
    double alpha = 1;
    // c_0 = inv_M;
    double c_0 = 1e-4;

    double PI = atan(1) * 4;
    double lambda = E * nu / (1 + nu) / (1 - 2 * nu);
    double mu = E / 2 / (1 + nu);
    double K_bulk = lambda + 2. / 3 * mu;
    double K_u = K_bulk + alpha * alpha / c_0;
    // fluid diffusivity coefficient
    double c_f = 1. / c_0 * k/mu_f * (K_bulk + 4. / 3 * mu) / (K_u + 4. / 3 * mu);
    // Skempton's coefficient
    double B = alpha / c_0 / K_u;
    // double B = 5.0 / 6; // Bin Wang Data
    // undrained Poisson's ratio
    double nu_u = (3*nu + alpha* B *(1-2*nu))/(3-alpha*B*(1-2*nu));
    // double nu_u = 0.44; // Bin Wang Data
    double cc = (1 - nu) / (nu_u - nu);
};

#endif