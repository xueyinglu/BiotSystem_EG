#ifndef LAMBDA_FUNCTION_H_
#define LAMBDA_FUNCTION_H_
#include "DealiiHeader.h"
using namespace std;
using namespace dealii;
class LambdaFunction: public Function<dim>{
    public:
        LambdaFunction():
        Function<dim>(1)
        {}

        virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const;
};

#endif