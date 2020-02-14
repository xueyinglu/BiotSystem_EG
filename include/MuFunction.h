#ifndef MU_FUNCTION_H_
#define MU_FUNCTION_H_
#include "DealiiHeader.h"
using namespace std;
using namespace dealii;
class MuFunction: public Function<dim>{
    public:
        MuFunction():
        Function<dim>(1)
        {}

        virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const;
};
#endif 