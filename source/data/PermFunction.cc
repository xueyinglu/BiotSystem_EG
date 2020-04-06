#include "PermFunction.h"
using namespace std;
double PermFunction::value(const Point<dim> &p, const unsigned int component) const
{
    double k = 1e-16;
    double x = p(0);
    double y = p(1);

    /* ARMA paper */ /*
    if (x >= 0.2 && x <= 0.8 && y >= 0.5 - 1. / 64 && y <= 0.5)
    {
        k = 1e-11;
    }
    if (x >= 0.2 && x <= 0.8 && y >= 0.25 - 1. / 64 && y <= 0.25)
    {
        k = 1e-11;
    }
    if (x >= 0.2 && x <= 0.8 && y >= 0.75 - 1. / 64 && y <= 0.75)
    {
        k = 1e-11;
    }
    */

    if (x >= 12. / 64 && x <= 52. / 64 && y >= 21. / 64 && y <= 22. / 64)
    {
        k = 1e-11;
    }
    if (x >= 12. / 64 && x <= 52. / 64 && y >= 42. / 64 && y <= 43. / 64)
    {
        k = 1e-11;
    }

    if (x >= 21. / 64 && x <= 22. / 64 && y >= 16. / 64 && y <= 27. / 64)
    {
        k = 1e-11;
    }
    if (x >= 42. / 64 && x <= 43. / 64 && y >= 16. / 64 && y <= 27. / 64)
    {
        k = 1e-11;
    }
    if (x >= 21. / 64 && x <= 22. / 64 && y >= 37. / 64 && y <= 48. / 64)
    {
        k = 1e-11;
    }
    if (x >= 42. / 64 && x <= 43. / 64 && y >= 37. / 64 && y <= 48. / 64)
    {
        k = 1e-11;
    }
    /*
    if (y<=0.5){
        k=1e-11;
    }
    */
    return k;
}