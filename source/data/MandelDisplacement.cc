#include "MandelDisplacement.h"

using namespace std;
void MandelDisplacement::vector_value(const Point<dim> &p, Vector<double> &values) const
{

    // solve tan(a_n) = (1-nu)/(nu_u-nu)*a_n numerically
    cout <<"cc = " << cc << endl;
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
            // cout << "error = " << error << endl;
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
        left = n* PI;
        right = (n + 0.5) * PI;
        middle = (left + right) / 2.0;
    }

    double x = p(0);
    double y = p(1);
    values(0) = F * nu / 2 / mu / a * x;
    values(1) = -F * (1 - nu) / 2 / mu / a * y;
    for (int i = 1; i < n_terms; i++)
    {
        double a_n = vector_a[i];
        values(0) += -F * nu_u / mu / a * sin(a_n) * cos(a_n) / (a_n - sin(a_n) * cos(a_n)) * exp(-a_n * a_n * c_f * t / a / a) * x + F / mu * cos(a_n) / (a_n - sin(a_n) * cos(a_n)) * sin(a_n * x / a) * exp(-a_n * a_n * c_f * t / a / a);
        values(1) += F * (1 - nu_u) / mu / a * sin(a_n) * cos(a_n) / (a_n - sin(a_n) * cos(a_n)) * exp(-a_n * a_n * c_f * t / a / a) * y;
    }
    if (t == 0)
    {
        values(0) = 0;
        values(1) = 1;
    }
    // cout <<"Mandel displacement ux=" <<values(0) <<endl;
    // cout <<"Mandel displacement uy=" <<values(1) <<endl;
}
void MandelDisplacement::gradient_value(const Point<dim> &point, Tensor<2, dim> &gradient) const
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
            // cout << "error = " << error << endl;
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
    gradient[0][0] = F * nu / 2 / mu / a ;
    gradient[0][1] = 0;
    gradient[1][0] = 0;
    gradient[1][1]= -F * (1 - nu) / 2 / mu / a ;
    for (int i = 1; i < n_terms; i++)
    {
        double a_n = vector_a[i];
        gradient[0][0] += -F * nu_u / mu / a * sin(a_n) * cos(a_n) / (a_n - sin(a_n) * cos(a_n)) * exp(-a_n * a_n * c_f * t / a / a)
                        + F / mu * cos(a_n) / (a_n - sin(a_n) * cos(a_n)) * a_n/a * cos(a_n * x / a) * exp(-a_n * a_n * c_f * t / a / a) ;
        gradient[1][1] += F * (1 - nu_u) / mu / a * sin(a_n) * cos(a_n) / (a_n - sin(a_n) * cos(a_n)) * exp(-a_n * a_n * c_f * t / a / a) ;
    }
    if (t == 0)
    {
        gradient[0][0] = 0;
        gradient[1][1] = 1;
    }
}

void MandelDisplacement::vector_value_list(const vector<Point<dim>> &points,
                                           vector<Vector<double>> &value_list) const
{
    Assert(value_list.size() == points.size(),
           ExcDimensionMismatch(value_list.size(), points.size()));

    const unsigned int n_points = points.size();

    for (unsigned int p = 0; p < n_points; ++p)
        vector_value(points[p], value_list[p]);
}