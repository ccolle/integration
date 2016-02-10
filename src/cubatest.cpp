#include <cstdlib>
#include <cassert>
#include <cstdio>
#include "integrator.hpp"
#include "cuba_wrappers.hpp"
#include "gsl_wrappers.hpp"
#include "cubature_wrappers.hpp"
#include <iostream>


int myintegrand_multicomp(const int* ndim, const double* x, const int* ncomp, double* f,void* param){
    //assert(*ndim==2 && *ncomp==3);
    for (int i=0; i<*ncomp;i++){
        f[i] = 1.;
    }
    return 0; // succes
}

void integrate_mcomp(struct Integrator &i){
    i.ndim  = 2;
    i.ncomp = 8;
    double f[i.ncomp];
    double e[i.ncomp];
    i.params = nullptr;
    i.integrand = myintegrand_multicomp;
    i.exec(f,e);
    for (int c=0;c<i.ncomp;c++){
        printf("[multicomp] component[%d] : result,error : (% .2e) , (% .2e)\n",c,f[c],e[c]);
    }
}

void test_cuba_raw(){ // test cuba without using the common integratior interface...
    int ndim  = 2;
    int ncomp = 8;
    integrand_t integr = myintegrand_multicomp;
    void* param = nullptr;
    int nvec  = 1;
    double epsrel = 1e-8;
    double epsabs = 0.;
    int flags = 0x00;
    int seed = 1234;
    int mineval = 10000;
    int maxeval = 100000;
    int nstart  = 1000;
    int nincrease = 1000;
    int nbatch  = 1000;
    int gridno  = 0;
    char* statefile = nullptr;
    void* spin = nullptr;
    int neval,fail;
    double integral[ncomp],error[ncomp],prob[ncomp];
    
    Vegas(ndim,ncomp,integr,param,nvec,epsrel,epsabs,flags,seed,
            mineval,maxeval,nstart,nincrease,nbatch,gridno,statefile,spin,
            &neval,&fail,integral,error,prob);
    
    for (int c=0;c<ncomp;c++){
        printf("[multicomp] component[%d] : result,error : (% .2e) , (% .2e)\n",c,integral[c],error[c]);
    }
}

void test_cuba(){
    Cuba_common c;
    c.nvec      = 1;
    c.epsabs    = 1e-5;
    c.epsrel    = 1e-5;
    c.flags     = 0x00;
    c.mineval   = 10000;
    c.maxeval   = 100000;
    c.statefile = nullptr;
    c.spin      = nullptr;
    
    Integrator_vegas vi(c);
    vi.seed   = 1234;
    vi.nstart = 1000;
    vi.nincrease = 100;
    vi.nbatch = 100;
    vi.gridno = 0;
    
    integrate_mcomp(vi);
    
    Integrator_cuhre ci(c);
    ci.key = 9;

    integrate_mcomp(ci);
}

int main(){
    test_cuba_raw();
    test_cuba();
    return 0;
}
