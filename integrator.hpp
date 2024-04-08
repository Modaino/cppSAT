#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <vector>

template <typename T>
class Integrator {
public:
    using RHSFunction = void (*)(double, const std::vector<T>&, std::vector<T>&);

    Integrator(RHSFunction rhs_function) : rhs_function(rhs_function) {}

    //virtual void solve_ivp(/* arguments */) = 0;

protected:
    RHSFunction rhs_function;
};

#endif // INTEGRATOR_HPP