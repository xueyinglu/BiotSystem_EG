#include "PermFunction.h"
using namespace std;
double PermFunction::value(const Point<dim> &p, const unsigned int component) const
{
    double k =1e-4;
    double x = p(0);
    double y = p(1);

    if (x >= 0.1 && x <=0.9 && y >= 0.5- 1./64 && y <= 0.5+ 1./64){
        k=10;
    }
        if (x >= 0.1 && x <=0.9 && y >= 0.25- 1./64 && y <= 0.25+ 1./64){
        k=10;
    }
        if (x >= 0.1 && x <=0.9 && y >= 0.75- 1./64 && y <= 0.75+ 1./64){
        k=10;
    }
    return k;
}