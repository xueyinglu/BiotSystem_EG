#include "PermFunction.h"
using namespace std;
double PermFunction::value(const Point<dim> &p, const unsigned int component) const
{
    double k = 1e-16;
    double x = p(0);
    double y = p(1);
    
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
    /*
    if (y<=0.5){
        k=1e-11;
    }
    */
    return k;
}