#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_cblas.h>


double pi = 4.0 * std::atan(1.0);

double s(double n, double x) {
    return std::sin(pi*(n+1.0)*x/10.0);
}

double ds(double n, double x) {
    return pi*(n+1.0)*std::cos(pi*(n+1.0)*x/10.0)/10.0;
}

double d2s(double n, double x) {
    return -pi*pi*(n+1.0)*(n+1.0)*std::sin(pi*(n+1.0)*x/10.0)/100.0;
}

double c(double n, double x) {
    return std::cos(pi*(n+0.5)*x/10.0);
}

double dc(double n, double x) {
    return -pi*(n+0.5)*std::sin(pi*(n+0.5)*x/10.0)/10.0;
}

double d2c(double n, double x) {
    return -pi*pi*(n+0.5)*(n+0.5)*std::cos(pi*(n+0.5)*x/10.0)/100.0;
}

double function_to_integrate_O_sin(double x, void* params) {
    int* parameters = (int*) params;
    int l = parameters[0];
    int n1 = parameters[1];
    int n2 = parameters[2];
    return s(n1, x)*(- x*x*d2s(n2, x) - (2.0*x + x*x)*ds(n2, x) + l*(l+1)*s(n2, x)) * std::exp(x);
}

double function_to_integrate_S_sin(double x, void* params) {
    int* parameters = (int*) params;
    int n1 = parameters[0];
    int n2 = parameters[1];
    return s(n1, x) * s(n2, x) * x*x;
}

double function_to_integrate_O_cos(double x, void* params) {
    int* parameters = (int*) params;
    int l = parameters[0];
    int n1 = parameters[1];
    int n2 = parameters[2];
    return c(n1, x)*(- x*x*d2c(n2, x) - (2.0*x + x*x)*dc(n2, x) + l*(l+1)*c(n2, x)) * std::exp(x);
}

double function_to_integrate_S_cos(double x, void* params) {
    int* parameters = (int*) params;
    int n1 = parameters[0];
    int n2 = parameters[1];
    return c(n1, x) * c(n2, x) * x*x;
}

int main() {
    int l_max = 50;
    int n_max_big = 200;

    double* O = new double[(l_max+1)*n_max_big*n_max_big];
    double* S = new double[(l_max+1)*n_max_big*n_max_big];

    #pragma omp parallel for
    for (int l = 0; l < l_max + 1; l++) {
        gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(10000);
        for (int n1 = 0; n1 < n_max_big; n1++) {
            for (int n2 = 0; n2 < n_max_big; n2++) {
                double result, error;

                gsl_function function_to_integrate_O_gsl;
                if (l == 0) {
                    function_to_integrate_O_gsl.function = &function_to_integrate_O_cos;
                } else {
                    function_to_integrate_O_gsl.function = &function_to_integrate_O_sin;
                }
                int parameters_O[3] = {l, n1, n2};
                function_to_integrate_O_gsl.params = parameters_O;
                gsl_integration_qags(&function_to_integrate_O_gsl, 0.0, 10.0, 1e-8, 1e-8, 10000, workspace, &result, &error);
                O[l*n_max_big*n_max_big + n1*n_max_big + n2] = result;

                gsl_function function_to_integrate_S_gsl;
                if (l == 0) {
                    function_to_integrate_S_gsl.function = &function_to_integrate_S_cos;
                } else {
                    function_to_integrate_S_gsl.function = &function_to_integrate_S_sin;
                }
                int parameters_S[2] = {n1, n2};
                function_to_integrate_S_gsl.params = parameters_S;
                gsl_integration_qags(&function_to_integrate_S_gsl, 0.0, 10.0, 1e-8, 1e-8, 10000, workspace, &result, &error);
                S[l*n_max_big*n_max_big + n1*n_max_big + n2] = result;
            }
        }
        gsl_integration_workspace_free(workspace);
    }

    // Save O and S to file
    std::ofstream file_O("O.txt");
    std::ofstream file_S("S.txt");

    file_O << std::setprecision(16);
    file_S << std::setprecision(16);

    for (int l = 0; l < l_max + 1; l++) {
        for (int n1 = 0; n1 < n_max_big; n1++) {
            for (int n2 = 0; n2 < n_max_big; n2++) {
                file_O << O[l*n_max_big*n_max_big + n1*n_max_big + n2] << " ";
                file_S << S[l*n_max_big*n_max_big + n1*n_max_big + n2] << " ";
            }
        }
    }

    file_O.close();
    file_S.close();

    delete[] O;
    delete[] S;
    return 0;
}
