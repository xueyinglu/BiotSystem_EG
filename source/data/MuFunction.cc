#include "MuFunction.h"
using namespace std;

double MuFunction::value(const Point<dim> &p, const unsigned int component) const
{
    double E = 5e6;  // Young's modulus
    double nu = 0.2; // Poisson's ratio
    double x = p(0);
    double y = p(1);
    
    if (x >= 0.2 && x <= 0.8 && y >= 0.5 - 1. / 64 && y <= 0.5)
    {
        E = 1e4;
        nu = 0.05;
    }
    if (x >= 0.2 && x <= 0.8 && y >= 0.25 - 1. / 64 && y <= 0.25)
    {
        E = 1e4;
        nu = 0.05;
    }
    if (x >= 0.2 && x <= 0.8 && y >= 0.75 - 1. / 64 && y <= 0.75)
    {
        E = 1e4;
        nu = 0.05;
    }
    
    double lame_mu = E / 2 / (1 + nu);
    return lame_mu;
}
