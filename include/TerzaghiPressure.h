#ifndef TERZAGHI_PRESSURE_H_
#define TERZAGHI_PRESSURE_H_
#include "DealiiHeader.h"
using namespace std;
using namespace dealii;
class TerzaghiPressure : public Function<dim>
{
public:
    TerzaghiPressure(double _t) : Function<dim>(1) { t = _t; };

    virtual double value(const Point<dim> &p, const unsigned int component) const;

private:
    double t;
    // Follows Phillips thesis
    double k = 1e-3;
    double E = 1e5;
    double nu = 0.2;
    double alpha = 1;
    // c_0 = inv_M;
    double c_0 = 0.1;
    double lambda = E * nu / (1 + nu) / (1 - 2 * nu);
    double mu = E / 2 / (1 + nu);
    double K_bulk = lambda + 2. / 3 * mu;
    double K_u = K_bulk + alpha * alpha / c_0;
    // fluid diffusivity coefficient
    double c_f = 1./c_0 * k * (K_bulk + 4. / 3 * mu) / (K_u + 4. / 3 * mu);

    double F = 1e3; // traction

    double PI = atan(1) * 4;
};
#endif