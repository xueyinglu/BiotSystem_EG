#include "MandelDisplacement.h"

using namespace std;
void MandelDisplacement::vector_value(const Point<dim> &p, Vector<double> &values) const
{

    double PI = atan(1) * 4;
    // Data from Bin Wang
    // computational domain is axb
    double a = 100;
    double b = 10;
    double k = 1e-13;
    double mu_f = 1e-3;
    double F = 5.94e8;
    double E = 5.94e9;
    double nu = 0.2;
    double alpha = 1;
    // c_0 = inv_M;
    double c_0 = 1.0 / 1.65e10;
    double lambda = E * nu / (1 + nu) / (1 - 2 * nu);
    double mu = E / 2 / (1 + nu);
    double K_bulk = lambda + 2. / 3 * mu;
    double K_u = K_bulk + alpha * alpha / c_0;
    // fluid diffusivity coefficient
    double c_f = 1. / c_0 * k / mu_f * (K_bulk + 4. / 3 * mu) / (K_u + 4. / 3 * mu);
    // Skempton's coefficient
    // double B = alpha / c_0 / K_u;
    double B = 5.0 / 6;
    // undrained Poisson's ratio
    // double nu_u = (3*nu + alpha* B *(1-2*nu))/(3-alpha*B*(1-2*nu));
    double nu_u = 0.44;
    double cc = (1 - nu) / (nu_u - nu);
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
    double PI = atan(1) * 4;
    // Data from Bin Wang
    // computational domain is axb
    double a = 100;
    double b = 10;
    double k = 1e-13;
    double mu_f = 1e-3;
    double F = 5.94e8;
    double E = 5.94e9;
    double nu = 0.2;
    double alpha = 1;
    // c_0 = inv_M;
    double c_0 = 1.0 / 1.65e10;
    double lambda = E * nu / (1 + nu) / (1 - 2 * nu);
    double mu = E / 2 / (1 + nu);
    double K_bulk = lambda + 2. / 3 * mu;
    double K_u = K_bulk + alpha * alpha / c_0;
    // fluid diffusivity coefficient
    double c_f = 1. / c_0 * k / mu_f * (K_bulk + 4. / 3 * mu) / (K_u + 4. / 3 * mu);
    // Skempton's coefficient
    // double B = alpha / c_0 / K_u;
    double B = 5.0 / 6;
    // undrained Poisson's ratio
    // double nu_u = (3*nu + alpha* B *(1-2*nu))/(3-alpha*B*(1-2*nu));
    double nu_u = 0.44;
    double cc = (1 - nu) / (nu_u - nu);
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