#ifndef GSL_INTEGRATOR_WRAPPERS
#define GSL_INTEGRATOR_WRAPPERS

#include "integrator.hpp"
#include <cassert>
#include <cstdlib>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <algorithm>

struct gsl_monte_integrand_wrapper {
    static double f(double* x,size_t nd, void* p){
        struct gsl_monte_integrand_wrapper m = *(struct gsl_monte_integrand_wrapper*) p;
        double res;
        int ndim = static_cast<int>(nd);
        m.integrand(&ndim,x,&m.ncomp,&res,m.params);
        return res;
    }
    int ncomp;
    void* params;
    integrand_t integrand;
};

struct gsl_1d_integrand_wrapper {
    static double f(double x,void* p){
        struct gsl_1d_integrand_wrapper m = *(struct gsl_1d_integrand_wrapper*) p;
        double res;
        int ncomp = 1; // always 1
        int ndim  = 1; // always 1
        m.integrand(&ndim,&x,&ncomp,&res,m.params);
        return res;
    }
    void* params;
    integrand_t integrand;
};

struct Integrator_gsl_monte_plain : Integrator {
    virtual int exec(double* integral,double* error){
        assert(ncomp==1); // gsl mc integration does only single component yo!
        gsl_monte_integrand_wrapper iw;
        iw.ncomp     = ncomp;
        iw.params    = params;
        iw.integrand = integrand;

        gsl_monte_function G;
        G.f      = gsl_monte_integrand_wrapper::f;
        G.dim    = ndim;
        G.params = (void*) &iw;
        
        double xl[ndim];
        double xu[ndim];
        std::fill(xl,xl+ndim,0.); // unit cube!
        std::fill(xu,xu+ndim,1.); // unit cube!
        
        const gsl_rng_type* T;
        gsl_rng* r;

        gsl_rng_env_setup();
        
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        
        gsl_monte_plain_state* s = gsl_monte_plain_alloc(ndim);
        gsl_monte_plain_integrate(&G,xl,xu,ndim,calls,r,s,integral,error);
        
        gsl_monte_plain_free(s);
        gsl_rng_free(r);
        return 0; // succes!
        }
    size_t calls;
};

struct Integrator_gsl_monte_vegas : Integrator {
    virtual int exec(double* integral,double* error){
        assert(ncomp==1); // gsl mc integration does only single component yo!
        gsl_monte_integrand_wrapper iw;
        iw.ncomp     = ncomp;
        iw.params    = params;
        iw.integrand = integrand;

        gsl_monte_function G;
        G.f      = gsl_monte_integrand_wrapper::f;
        G.dim    = ndim;
        G.params = (void*) &iw;
        
        double xl[ndim];
        double xu[ndim];
        std::fill(xl,xl+ndim,0.); // unit cube!
        std::fill(xu,xu+ndim,1.); // unit cube!
        
        const gsl_rng_type* T;
        gsl_rng* r;

        gsl_rng_env_setup();
        
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        
        gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(ndim);
        //gsl_monte_vegas_integrate(&G,xl,xu,ndim,calls,r,s,integral,error);
        
        // integrate with low number of calls to map the function
        gsl_monte_vegas_integrate(&G, xl, xu, ndim, warm_up_calls, r, s, integral, error);

        // iterate integration until convergence is reached
        while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > delta_chi_squared) {
            gsl_monte_vegas_integrate(&G, xl, xu, ndim, in_iteration_calls, r, s, integral, error);
        }

        gsl_monte_vegas_free(s);
        gsl_rng_free(r);
        return 0; // succes!
        }
    size_t warm_up_calls; // small number of calls to map the function
    size_t in_iteration_calls; // number of calls in integration during iteration to convergence
    double delta_chi_squared; // convergence criterium: if the calculated chi-squared per dof 
                              // lies within this distance of the value 1, the iteration is stopped.
};

struct Integrator_gsl_monte_miser : Integrator {
    virtual int exec(double* integral,double* error){
        assert(ncomp==1); // gsl mc integration does only single component yo!
        gsl_monte_integrand_wrapper iw;
        iw.ncomp     = ncomp;
        iw.params    = params;
        iw.integrand = integrand;

        gsl_monte_function G;
        G.f      = gsl_monte_integrand_wrapper::f;
        G.dim    = ndim;
        G.params = (void*) &iw;
        
        double xl[ndim];
        double xu[ndim];
        std::fill(xl,xl+ndim,0.); // unit cube!
        std::fill(xu,xu+ndim,1.); // unit cube!
        
        const gsl_rng_type* T;
        gsl_rng* r;

        gsl_rng_env_setup();
        
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        
        gsl_monte_miser_state* s = gsl_monte_miser_alloc(ndim);
        gsl_monte_miser_integrate(&G,xl,xu,ndim,calls,r,s,integral,error);
        
        gsl_monte_miser_free(s);
        gsl_rng_free(r);
        return 0; // succes!
        }
    size_t calls;
};

struct Integrator_gsl_integration_qag : Integrator {
    virtual int exec(double* integral,double* error){
        assert(ncomp==1 && ndim == 1); // gsl 1d integration...
        gsl_1d_integrand_wrapper iw;
        iw.integrand = integrand;
        iw.params = params;

        gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
        gsl_function F;
        F.function= gsl_1d_integrand_wrapper::f;
        F.params = (void*) &iw;

        gsl_integration_qag(&F,0,1,epsabs,epsrel,limit,key,w,integral,error);
        gsl_integration_workspace_free(w);
        return 0; // succes!
    }
    size_t limit;
    double epsabs,epsrel;
    int key;
};

#endif // GSL_INTEGRATOR_WRAPPERS
