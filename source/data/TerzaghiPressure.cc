#include "TerzaghiPressure.h"
using namespace std;

double TerzaghiPressure::value(const Point<dim> &p, const unsigned int component) const
{
    double k = 1e-3;
    double E = 1e5;
    double nu = 0.2;
    double alpha = 1;
    double inv_M = 0.1;
    double lambda = E * nu / (1 + nu) / (1 - 2 * nu);
    double mu = E / 2 / (1 + nu);
    double K_bulk = lambda + 2. / 3 * mu;
    double K_u = K_bulk + alpha * alpha * inv_M;
    double c_f = inv_M * k * (K_bulk + 4. / 3 * mu) / (K_u + 4. / 3 * mu);

    double F = 1e3; // traction

    double PI = atan(1) * 4;
    double pressure = 0;

    double y = p(1);
    int n_terms = 100000;

    for (int i = 0; i < n_terms; i++)
    {
        double M = PI * (2 * i + 1) / 2;
        pressure += alpha * F * inv_M / (K_u + 4. / 3 * mu) * 2 / M * sin(M * y) * exp(-M * M * c_f * t);
        // pressure +=  2 / M * sin(M * y) * exp(-M * M * c_f * t);
    }
    if (t == 0)
    {
        pressure = 0;
    }
    return pressure;
}