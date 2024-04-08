#include "Cash_Karp.hpp"

/* 
 * Function originaly from Numerical Recipies in C by William H. Press & Saul A. Teukolsky
 * impemented in c++ by Áron Vízkeleti on 2024-04-07
 * 
 * Given values for n variables y[1..n] and their derivatives dydx[1..n] known at x,use the fifth-order Cash-Karp Runge-Kutta method to advance the solution over an interval h and return the incremented variables as yout[1..n]. Also return an estimate of the local truncation error in yout using the embedded fourth-order method. The user supplies the routine derivs(x,y,dydx), which returns derivatives dydx at x.
 * 
 */


template<typename T>
void CashKarp<T>::rkck(std::vector<T>& y, std::vector<T>& dydx, double x, double h, std::vector<T>& yout, std::vector<T>& yerr) {
    int n = y.size();
    std::vector<T> ak2(n), ak3(n), ak4(n), ak5(n), ak6(n), ytemp(n);

    for (int i = 0; i < n; ++i)
        ytemp[i] = y[i] + b21 * h * dydx[i];
    this->rhs_function(x + a2 * h, ytemp, ak2);

    for (int i = 0; i < n; ++i)
        ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
    this->rhs_function(x + a3 * h, ytemp, ak3);

    for (int i = 0; i < n; ++i)
        ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
    this->rhs_function(x + a4 * h, ytemp, ak4);

    for (int i = 0; i < n; ++i)
        ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
    this->rhs_function(x + a5 * h, ytemp, ak5);

    for (int i = 0; i < n; ++i)
        ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
    this->rhs_function(x + a6 * h, ytemp, ak6);

    for (int i = 0; i < n; ++i)
        yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);

    for (int i = 0; i < n; ++i)
        yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);
}


template<typename T>
void CashKarp<T>::rkqs(std::vector<T>& y, std::vector<T>& dydx, int n, double& x, double htry, double eps, std::vector<T>& yscal, double& hdid, double& hnext) {
    //void rkck(std::vector<T>& y, std::vector<T>& dydx, double x, double h, std::vector<T>& yout, std::vector<T>& yerr);

    std::vector<T> yerr(n), ytemp(n);
    double errmax, h, htemp, xnew;

    h = htry;
    for (;;) {
        rkck(y, dydx, x, h, ytemp, yerr); // Take a step

        errmax = 0.0;
        for (int i = 0; i < n; ++i)
            errmax = std::max(errmax, fabs(yerr[i] / yscal[i]));
        errmax /= eps;

        if (errmax <= 1.0)
            break; // Step succeeded. Compute size of next step

        htemp = SAFETY * h * pow(errmax, PSHRNK);
        h = (h >= 0.0) ? std::max(htemp, 0.1 * h) : std::min(htemp, 0.1 * h); // No more than a factor of 10

        xnew = x + h;
        if (xnew == x)
            throw std::runtime_error("Stepsize underflow in rkqs");
    }

    if (errmax > ERRCON)
        hnext = SAFETY * h * pow(errmax, PGROW);
    else
        hnext = 5.0 * h; // No more than a factor of 5 increase

    x += (hdid = h);
    for (int i = 0; i < n; ++i)
        y[i] = ytemp[i];
}


// Explicit instantiation for double
template class CashKarp<double>;
