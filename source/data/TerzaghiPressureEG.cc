#include "TerzaghiPressureEG.h"
using namespace std;

void TerzaghiPressureEG::vector_value(const Point<dim> &p, Vector<double> &values) const
{
    values[1] = 0;
    double pressure = 0;

    double y = p(1);
    int n_terms = 100000;

    for (int i = 0; i < n_terms; i++)
    {
        double M = PI * (2 * i + 1) / 2;
        pressure += alpha * F /c_0 / (K_u + 4. / 3 * mu) * 2 / M * sin(M * y) * exp(-M * M * c_f * t);
        // pressure +=  2 / M * sin(M * y) * exp(-M * M * c_f * t);
    }
    if (t == 0)
    {
        pressure = 0;
    }
    values[0] = pressure;
}