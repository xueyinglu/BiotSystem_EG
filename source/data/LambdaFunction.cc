#include "LambdaFunction.h"
using namespace std;

double LambdaFunction::value(const Point<dim> &p, const unsigned int component) const
{
    double E = 5e6;  // Young's modulus
    double nu = 0.2; // Poisson's ratio
    double x = p(0);
    double y = p(1);
    if (x >= 0.1 && x <= 0.9 && y >= 0.5 - 1. / 64 && y <= 0.5 + 1. / 64)
    {
        E = 1e4;
        nu = 0.05;
    }

    if (x >= 0.1 && x <= 0.9 && y >= 0.25 - 1. / 64 && y <= 0.25 + 1. / 64)
    {
        E = 1e4;
        nu = 0.05;
    }
    if (x >= 0.1 && x <= 0.9 && y >= 0.75 - 1. / 64 && y <= 0.75 + 1. / 64)
    {
        E = 1e4;
        nu = 0.05;
    }

    double lame_lambda = E * nu / (1 + nu) / (1 - 2 * nu);
    return lame_lambda;
}
