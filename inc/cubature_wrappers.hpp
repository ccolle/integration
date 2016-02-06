#ifndef CUBATURE_WRAPPERS_HPP
#define CUBATURE_WRAPPERS_HPP

#include "integrator.hpp"
#include "cubature.h"
#include <algorithm> // std::fill

struct cubature_integrand_wrapper{
    static int f(unsigned ndim,const double* x,void* fdata,unsigned ncomp,double* fval){
        struct cubature_integrand_wrapper c = *(struct cubature_integrand_wrapper*) fdata;
        int nd = static_cast<int>(ndim);
        int nc = static_cast<int>(ncomp);
        return c.integrand(&nd,x,&nc,fval,c.params);
    }
    void* params;
    integrand_t integrand;
};

struct Integrator_hcubature : Integrator {
    virtual int exec(double* integral,double* error){
        cubature_integrand_wrapper c;
        c.params    = params;
        c.integrand = integrand;

        double xl[ndim];
        double xu[ndim];
        std::fill(xl,xl+ndim,0.); // unit cube!
        std::fill(xu,xu+ndim,1.); // unit cube!

        return hcubature(ncomp,cubature_integrand_wrapper::f,(void*) &c,
                        ndim,xl,xu,maxEval,epsabs,epsrel,err_norm,integral,error);
        }
    size_t maxEval;
    double epsabs,epsrel;
    error_norm err_norm; //<** options are: [ERROR_L1,ERROR_L2,ERROR_LINF,ERROR_INDIVIDUAL,ERROR_PAIRED] 
};

struct Integrator_pcubature : Integrator {
    virtual int exec(double* integral,double* error){
        cubature_integrand_wrapper c;
        c.params    = params;
        c.integrand = integrand;

        double xl[ndim];
        double xu[ndim];
        std::fill(xl,xl+ndim,0.); // unit cube!
        std::fill(xu,xu+ndim,1.); // unit cube!
        
        return pcubature(ncomp,cubature_integrand_wrapper::f,(void*) &c,
                        ndim,xl,xu,maxEval,epsabs,epsrel,err_norm,integral,error);
        }
    size_t maxEval;
    double epsabs,epsrel;
    error_norm err_norm; //<** options are: [ERROR_L1,ERROR_L2,ERROR_LINF,ERROR_INDIVIDUAL,ERROR_PAIRED] 
};





#endif // CUBATURE_WRAPPERS_HPP
