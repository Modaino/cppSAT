#ifndef CASH_KARP_HPP
#define CASH_KARP_HPP

#include "integrator.hpp"
#include <cmath>
#include <stdexcept>

template <typename T>
class CashKarp : public Integrator<T> {
public:
    using Integrator<T>::Integrator;

    void rkck(std::vector<T>& y, std::vector<T>& dydx, double x, double h, std::vector<T>& yout, std::vector<T>& yerr);
    void rkqs(std::vector<T>& y, std::vector<T>& dydx, int n, double& x, double htry, double eps, std::vector<T>& yscal, double& hdid, double& hnext);


private:
    static constexpr double a2 = 0.2;
    static constexpr double a3 = 0.3;
    static constexpr double a4 = 0.6;
    static constexpr double a5 = 1.0;
    static constexpr double a6 = 0.875;
    static constexpr double b21 = 0.2;
    static constexpr double b31 = 3.0 / 40.0;
    static constexpr double b32 = 9.0 / 40.0;
    static constexpr double b41 = 0.3;
    static constexpr double b42 = -0.9;
    static constexpr double b43 = 1.2;
    static constexpr double b51 = -11.0 / 54.0;
    static constexpr double b52 = 2.5;
    static constexpr double b53 = -70.0 / 27.0;
    static constexpr double b54 = 35.0 / 27.0;
    static constexpr double b61 = 1631.0 / 55296.0;
    static constexpr double b62 = 175.0 / 512.0;
    static constexpr double b63 = 575.0 / 13824.0;
    static constexpr double b64 = 44275.0 / 110592.0;
    static constexpr double b65 = 253.0 / 4096.0;
    static constexpr double c1 = 37.0 / 378.0;
    static constexpr double c3 = 250.0 / 621.0;
    static constexpr double c4 = 125.0 / 594.0;
    static constexpr double c6 = 512.0 / 1771.0;
    static constexpr double dc1 = c1-2825.0/27648.0;
    static constexpr double dc3 = c3-18575.0/48384.0;
    static constexpr double dc4 = c4-13525.0/55296.0;
    static constexpr double dc5 = -277.00 / 14336.0;
    static constexpr double dc6=c6-0.25;

    static constexpr double SAFETY = 0.9;
    static constexpr double PGROW = -0.2;
    static constexpr double PSHRNK = -0.25;
    static constexpr double ERRCON = 1.89e-4;
};

#endif // CASH_KARP_HPP
