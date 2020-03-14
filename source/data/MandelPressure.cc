#include "MandelPressure.h"
using namespace std;
double MandelPressure::value(const Point<dim> &p, const unsigned int component) const
{
    // Data from Bin Wang
    double PI = atan(1) * 4;
    // computational domain is axb
    double a = 100;
    double b = 10;
    double k = 1e-13;
    double mu_f =1e-3;
    double F = 5.94e8;
    double E = 5.94e9;
    double nu = 0.2;
    double alpha = 1;
    // c_0 = inv_M;
    double c_0 = 1 / 1.65e10;
    double lambda = E * nu / (1 + nu) / (1 - 2 * nu);
    double mu = E / 2 / (1 + nu);
    double K_bulk = lambda + 2. / 3 * mu;
    double K_u = K_bulk + alpha * alpha / c_0;
    // fluid diffusivity coefficient
    double c_f = 1. / c_0 * k/mu_f * (K_bulk + 4. / 3 * mu) / (K_u + 4. / 3 * mu);
    // Skempton's coefficient
    // double B = alpha / c_0 / K_u;
    double B = 5.0 / 6;
    // undrained Poisson's ratio
    // double nu_u = (3*nu + alpha* B *(1-2*nu))/(3-alpha*B*(1-2*nu));
    double nu_u = 0.44;
    double cc = (1 - nu) / (nu_u - nu);
    // solve tan(a_n) = (1-nu)/(nu_u-nu)*a_n numerically
    vector<double> vector_a;
    vector_a.push_back(0);
    double n_terms = 50;
    double left = 0;
    double right = PI / 2;
    double middle = (left + right) / 2;
    double error_mod = 1;
    double x0;
    for (int n = 1; n < n_terms; n++)
    {

        while (error_mod > 1e-8)
        {
            double error = cc * middle - tan(middle);
            if (error > 0.0)
            {
                left = middle;
            }
            else if (error < 0.0)
            {
                right = middle;
            }
            else
            {
                left = middle;
                right = middle;
            }
            middle = (left + right) / 2.0;
            error_mod = abs(tan(middle) - cc * middle);
            x0 = middle;
        }
        vector_a.push_back(x0);
        error_mod = 1.0;
        left = x0 + PI;
        right = (n + 0.5) * PI;
        middle = (left + right) / 2.0;
    }
    double x = p(0);
    double pressure = 0;
    for (int i = 1; i < n_terms; i++)
    {
        double a_n =vector_a[i];
        pressure += 2 * F * B * (1 + nu_u) / 3 / a * sin(a_n) / (a_n - sin(a_n) * cos(a_n)) * (cos(a_n * x / a) - cos(a_n)) * exp(-a_n * a_n * c_f * t / a / a);
    }
    if (t==0){
        pressure =0;
    }
    return pressure;
}