#include "InitialPressure.h"
using namespace std;

double InitialPressure::value(const dealii::Point<dim> &p,
                                    const unsigned int component) const
{

    double t = this->get_time();

    double PI = atan(1)*4;
    /*
    if (component == 0)
    {
        return_value = cos(x - y);
        //return_value = exp(-x-y*y);
    }
    else if (component == 1)
    {
        return_value = cos(x - y);
        //return_value = exp(-x-y*y);
    }
    else
    {
        cout << "ERROR " << component << endl;
        exit(0);
    }
    */
    return sin(PI*p(0))*sin(PI*p(1));
}

void InitialPressure::grad_list(const std::vector<Point<dim>> & points, vector<Tensor<1, dim>> &values) const{
    
    Assert (values.size() == points.size(),
            ExcDimensionMismatch (values.size(), points.size()));

    const unsigned int n_points = points.size();
    double PI = atan(1) * 4;
    for (unsigned int p=0; p<n_points; ++p){
        values[p][0] = PI * cos(PI * points[p](0)) * sin(PI * points[p](1));
        values[p][1] = PI * sin(PI * points[p](0)) * cos(PI * points[p](1));
    }
}