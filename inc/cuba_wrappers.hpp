#ifndef CUBA_INTEGRATOR_WRAPPERS_HPP
#define CUBA_INTEGRATOR_WRAPPERS_HPP

#include "integrator.hpp"
#include "cuba.h"

/** A cuba common struct which holds the common
  * settings for all four integration routines in
  * the cuba library. This way you can more easily
  * try out different integration routines, while
  * keeping the common settings the same.
  *
  * c++11 is required for the default values setting
  * without constructor...
  */
struct Cuba_common : Integrator {
    // common settings for the cuba integrators...
    int nvec        = 1 ;
    double epsrel   = 0.;
    double epsabs   = 0.;
    int flags       = 0x00;;
    int mineval     = -1;
    int maxeval     = -1;
    char* statefile = nullptr;
    void* spin      = nullptr;
    int neval       = 0;
    int fail        = 0;
    double* prob    = nullptr;
    virtual int exec(double*,double*){ throw; } /**< this should never get called **/
};

struct Integrator_vegas : Cuba_common  {
    Integrator_vegas() : Cuba_common () {}
    Integrator_vegas(Cuba_common &c) : Cuba_common(c) {}  
    virtual int exec(double* integral,double* error){
        Vegas(ndim,ncomp,integrand,params,nvec,
              epsrel,epsabs,flags,seed,mineval,
              maxeval,nstart,nincrease,nbatch,gridno,
              statefile,spin,&neval,&fail,integral,error,prob);
        return 0; // "succes"
    }
    
    int seed,nstart,nincrease,nbatch,gridno;
};

struct Integrator_suave : Cuba_common {
    Integrator_suave() : Cuba_common() {}
    Integrator_suave(Cuba_common& c) : Cuba_common(c) {}

    int nnew,nmin,nregions,seed;
    double flatness;

    virtual int exec(double* integral,double* error){
        Suave(ndim,ncomp,integrand,params,nvec,
            epsrel,epsabs,flags,seed,mineval,
            maxeval,nnew,nmin,flatness,statefile,
            spin,&nregions,&neval,&fail,integral,
            error,prob);
            return 0;
    }
};

struct Integrator_cuhre : Cuba_common {
    Integrator_cuhre() : Cuba_common() {}
    Integrator_cuhre(Cuba_common& c) : Cuba_common(c) {}

    int key,nregions;

    virtual int exec(double* integral,double* error){
        Cuhre(ndim,ncomp,integrand,params,nvec,
            epsrel,epsabs,flags,
            mineval,maxeval,key,
            statefile,spin,&nregions,&neval,&fail,
            integral,error,prob);
            return 0;
    }
};

struct Integrator_divonne : Cuba_common {
    Integrator_divonne() : Cuba_common() {}
    Integrator_divonne( Cuba_common& c) : Cuba_common(c) {}

    int key1,key2,key3,maxpass,ngiven,ldxgiven,nextra,seed,nregions;
    double border,maxchisq,mindeviation;
    double* xgiven;
    peakfinder_t peakfinder;

    virtual int exec(double* integral,double* error){
        Divonne(ndim,ncomp,integrand,params,nvec,
            epsrel,epsabs,flags,seed,
            mineval,maxeval,key1,key2,key3,
            maxpass,border,maxchisq,mindeviation,
            ngiven,ldxgiven,xgiven,nextra,peakfinder,
            statefile,spin,&nregions,&neval,&fail,
            integral,error,prob);
        return 0;
    }
};
#endif // CUBA_INTEGRATOR_WRAPPERS
