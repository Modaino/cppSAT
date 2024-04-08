#include "Cash_Karp.hpp"
#include <iostream>

// Lorenz system's right-hand side function
void lorenz_rhs(double t, const std::vector<double>& state, std::vector<double>& result) {
    double sigma = 10.0;
    double rho = 28.0;
    double beta = 2.66667;

    result[0] = sigma * (state[1] - state[0]);
    result[1] = state[0] * (rho - state[2]) - state[1];
    result[2] = state[0] * state[1] - beta * state[2];
}

int main() {
    // Initial conditions
    std::vector<double> y = {1.0, 1.0, 1.0};
    std::vector<double> dydx = {0.0, 0.0, 0.0};
    double x = 0.0;
    double htry = 0.01;
    double eps = 1e-6;
    std::vector<double> yscal = {1.0, 1.0, 1.0};
    double hdid, hnext;

    // Create instance of CashKarp integrator
    CashKarp<double> integrator(lorenz_rhs);

    // Perform one step using rkck
    std::vector<double> yout(3), yerr(3);
    integrator.rkck(y, dydx, x, htry, yout, yerr);

    // Print the result
    std::cout << "Result after one step:" << std::endl;
    std::cout << "x: " << x << std::endl;
    std::cout << "yout: ";
    for (const auto& val : yout)
        std::cout << val << " ";
    std::cout << std::endl;

    return 0;
}
