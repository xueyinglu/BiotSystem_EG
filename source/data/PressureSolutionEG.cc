#include "PressureSolutionEG.h"
using namespace std;

void PressureSolutionEG::vector_value(const Point<dim> &p, Vector<double> & values) const{
    double k = 0.05;
    double alpha = 0.75;
    double inv_M = 3. / 28;
    double PI = atan(1) * 4;
    double A = 2 * PI * PI * k / (alpha + inv_M);
    values[0] = exp(-A * t) * sin(PI * p(0)) * sin(PI * p(1));
    values[1] = 0;
}