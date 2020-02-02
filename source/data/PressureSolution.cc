#include "PressureSolution.h"
using namespace std;

double PressureSolution::value(const Point<dim> &p, const unsigned int component) const
{
    double k = 0.05;
    double alpha = 0.75;
    double inv_M = 3. / 28;
    double PI = atan(1) * 4;
    double A = 2 * PI * PI * k / (alpha + inv_M);
    return exp(-A * t) * sin(PI * p(0)) * sin(PI * p(1));
}