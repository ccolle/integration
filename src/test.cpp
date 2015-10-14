#include <cstdlib>
#include <cassert>
#include <cstdio>
#include "integrator.hpp"
#include "cuba_wrappers.hpp"
#include "gsl_wrappers.hpp"
#include "cubature_wrappers.hpp"


/********************************************/
/** this is the applied part of your code  **/
/** it should not use anything outside of  **/
/** Integrator for the integrations...     **/
/********************************************/

struct myintegrand_params{
    int a;
};
int myintegrand1d(const int* ndim, const double* x,const int* ncomp,double* f,void* param){
    assert(*ndim==1 && *ncomp==1);
    myintegrand_params p = *(struct myintegrand_params*) param;
    f[0] = p.a*x[0];
    return 0; // succes!
}

int myintegrand(const int* ndim,const double* x,const int* ncomp,double* f,void* param){
    assert(*ndim==2 && *ncomp==1);
    myintegrand_params p = *(struct myintegrand_params*) param;
    f[0] = p.a*x[0]*x[1]*exp(-sin(x[0])); // should give a/4
    return 0; // succes
}
void integrate1d(struct Integrator& i){
    double f;
    double e;
    i.ndim=1;
    i.ncomp=1;
    i.integrand = myintegrand1d;
    struct myintegrand_params p;
    p.a = 2;
    i.params = (void*) &p;
    i.exec(&f,&e);
    printf("result,error :   %11.6e, %8.3e\n",f,e);
}

void integrate(struct Integrator& i){
    double f;
    double e;
    i.ndim = 2;
    i.ncomp = 1;
    i.integrand = myintegrand;
    struct myintegrand_params p;
    p.a = 2;
    i.params = (void*) &p;
    i.exec(&f,&e);
    printf("result,error :   %11.6e, %8.3e\n",f,e);
}
/********************************************/
/**          this is the end of it         **/
/********************************************/

void test_cuba(){
    Cuba_common c;
    c.nvec    = 1;
    c.epsabs  = 1e-5;
    c.epsrel  = 1e-5;
    c.flags   = 0x00;
    c.mineval = 10000;
    c.maxeval = 100000;
    c.statefile = nullptr;
    c.spin      = nullptr;
    c.prob      = new double(0.);
    
    Integrator_vegas vi(c);
    vi.seed   = 1234;
    vi.nstart = 1000;
    vi.nincrease = 100;
    vi.nbatch = 100;
    vi.gridno = 0;

    Integrator_cuhre ci(c);
    ci.key = 9;


    integrate(vi);
    integrate(ci);
    delete c.prob;
}

void test_gsl(){
    Integrator_gsl_monte_plain p;     
    p.calls = 100000; // max number of evals
    integrate(p);

    Integrator_gsl_monte_vegas v;
    v.calls = 100000;
    integrate(v);

    Integrator_gsl_monte_miser m;
    m.calls = 100000;
    integrate(m);
}

void test_cubature(){
    Integrator_hcubature hc;
    hc.maxEval = 100000;
    hc.epsabs  = 0.;
    hc.epsrel  = 1e-4;
    
    integrate(hc);

    Integrator_pcubature pc;
    pc.maxEval = 100000;
    pc.epsabs  = 0.;
    pc.epsrel  = 1e-4;

    integrate(pc);
}

void test_1d(){
    Cuba_common c;
    c.nvec    = 1;
    c.epsabs  = 1e-5;
    c.epsrel  = 1e-5;
    c.flags   = 0x00;
    c.mineval = 10000;
    c.maxeval = 100000;
    c.statefile = nullptr;
    c.spin      = nullptr;
    c.prob      = new double(0.);
    
    Integrator_vegas vi(c);
    vi.seed   = 1234;
    vi.nstart = 1000;
    vi.nincrease = 100;
    vi.nbatch = 100;
    vi.gridno = 0;

    Integrator_cuhre ci(c);
    ci.key = 9;
    
    integrate1d(vi);
    integrate1d(ci);

    delete c.prob;

    Integrator_gsl_monte_plain p;
    p.calls = 100000;
    integrate1d(p);

    Integrator_gsl_integration_qag q;
    q.limit  = 1000;
    q.epsabs = 1e-5;
    q.epsrel = 1e-5;
    q.key    = 6;
    integrate1d(q);

    Integrator_hcubature hc;
    hc.maxEval = 100000;
    hc.epsabs = 1e-5;
    hc.epsrel = 1e-5;
    integrate1d(hc);

    Integrator_pcubature pc;
    pc.maxEval = 100000;
    pc.epsabs = 1e-5;
    pc.epsrel = 1e-5;
    integrate1d(pc);

}

int main(){
    printf("cuba     : \n"); test_cuba();
    printf("gsl      : \n"); test_gsl();
    printf("cubature : \n"); test_cubature();
    printf("\n*** 1d test ***\n");
    test_1d();
    return 0;
}
