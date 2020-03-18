#include "MandelPressure.h"
using namespace std;
double MandelPressure::value(const Point<dim> &p, const unsigned int component) const
{
    // solve tan(a_n) = (1-nu)/(nu_u-nu)*a_n numerically
    vector<double> vector_a;
    vector_a.push_back(0);
    double n_terms = 50;
    double left = 0;
    double right = PI / 2;
    double middle = (left + right) / 2;
    double error_mod = 1;
    double x0;
    for (int n = 1; n < n_terms; n++)
    {

        while (error_mod > 1e-8)
        {
            double error = cc * middle - tan(middle);
            if (error > 0.0)
            {
                left = middle;
            }
            else if (error < 0.0)
            {
                right = middle;
            }
            else
            {
                left = middle;
                right = middle;
            }
            middle = (left + right) / 2.0;
            error_mod = abs(tan(middle) - cc * middle);
            x0 = middle;
        }
        vector_a.push_back(x0);
        error_mod = 1.0;
        left = x0 + PI;
        right = (n + 0.5) * PI;
        middle = (left + right) / 2.0;
    }
    double x = p(0);
    double pressure = 0;
    for (int i = 1; i < n_terms; i++)
    {
        double a_n =vector_a[i];
        pressure += 2 * F * B * (1 + nu_u) / 3 / a * sin(a_n) / (a_n - sin(a_n) * cos(a_n)) * (cos(a_n * x / a) - cos(a_n)) * exp(-a_n * a_n * c_f * t / a / a);
    }
    if (t==0){
        pressure =0;
    }
    return pressure;
}

void MandelPressure::gradient_value(const Point<dim> & point, Tensor<1,dim> & gradient) const{
    // solve tan(a_n) = (1-nu)/(nu_u-nu)*a_n numerically
    vector<double> vector_a;
    vector_a.push_back(0);
    double n_terms = 50;
    double left = 0;
    double right = PI / 2;
    double middle = (left + right) / 2;
    double error_mod = 1;
    double x0;
    for (int n = 1; n < n_terms; n++)
    {

        while (error_mod > 1e-8)
        {
            double error = cc * middle - tan(middle);
            if (error > 0.0)
            {
                left = middle;
            }
            else if (error < 0.0)
            {
                right = middle;
            }
            else
            {
                left = middle;
                right = middle;
            }
            middle = (left + right) / 2.0;
            error_mod = abs(tan(middle) - cc * middle);
            x0 = middle;
        }
        vector_a.push_back(x0);
        error_mod = 1.0;
        left = x0 + PI;
        right = (n + 0.5) * PI;
        middle = (left + right) / 2.0;
    }
    double x = point(0);
     gradient[0] =0;
     gradient[1] = 0;
    for (int i = 1; i < n_terms; i++)
    {
        double a_n =vector_a[i];
        gradient[0] += 2 * F * B * (1 + nu_u) / 3 / a * (-a_n/a)* sin(a_n) / (a_n - sin(a_n) * cos(a_n)) * sin(a_n * x / a) * exp(-a_n * a_n * c_f * t / a / a);
    }
    if (t==0){
        gradient[0] = 0;
    }

}