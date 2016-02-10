#ifndef CUBA_INTEGRATOR_WRAPPERS_HPP
#define CUBA_INTEGRATOR_WRAPPERS_HPP

#include "integrator.hpp"
#include "cuba.h"
#include <iostream>
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
    int flags       = 0x00;
    int mineval     = -1;
    int maxeval     = -1;
    char* statefile = nullptr;
    void* spin      = nullptr;
    double* prob    = nullptr; // you have to set this to an array of size ncomp (number of components)

    int neval       = 0; // you don't need to set this
    int fail        = 0; // you don't need to set this
    /**
     * @brief init, this makes sure prob array is allocated with size ncomp (not done @ construction time as ncomp is unknown then)
     */
    void init(){ delete prob; prob = new double[ncomp]; }// this makes sure the prob array is allocated with size ncomp
    virtual int exec(double*,double*){ throw; }          // this should never get called
    /**
      * @brief destructor, make sure the prob array is deleted correctly.
      */
    ~Cuba_common(){ delete prob; }
};

struct Integrator_vegas : Cuba_common  {
    Integrator_vegas() : Cuba_common () {}
    Integrator_vegas(Cuba_common &c) : Cuba_common(c) {}  
    virtual int exec(double* integral,double* error){
        Cuba_common::init();
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

    int nnew,nmin,seed; /**< input parameters **/
    double flatness;    /**< input parameters **/
    int nregions;       /**< output parameters **/
    virtual int exec(double* integral,double* error){
        Cuba_common::init();
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
        Cuba_common::init();
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
        Cuba_common::init();
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
