#include <cstdlib>
#include <cassert>
#include <cstdio>
#include "integrator.hpp"
#include "cuba_wrappers.hpp"
#include "gsl_wrappers.hpp"
#include "cubature_wrappers.hpp"
#include <iostream>
#include <string>
#include <cmath>

struct integrand_params{
    int a;
};

int myintegrand0(const int* ndim, const double* x,const int* ncomp,double* f,void* param){
    assert(*ndim==1 && *ncomp==1);
    integrand_params p = *(struct integrand_params*) param;
    f[0] = p.a*x[0];
    return 0; // succes!
}
/** the analytical result of the integration of myintegrand0 **/
int myintegrand0_integral(const int* ncomp,double* f, void* param){
    assert(*ncomp==1);
    integrand_params p = *(struct integrand_params*) param;
    f[0] = 0.5*p.a;
    return 0; // succes!
} 

int myintegrand1(const int* ndim,const double* x,const int* ncomp,double* f,void* param){
    assert(*ndim==1 && *ncomp==1);
    f[0] = exp(sin(x[0])); // should give 1.63187...
    return 0; // succes
}
/** the analytical result of the integration of myintegrand_1 **/
int myintegrand1_integral(const int* ncomp,double* f,void* param){
    assert(*ncomp==1);
    f[0] = 1.63187;
    return 0; // succes
}


int myintegrand_multicomp(const int* ndim, const double* x, const int* ncomp, double* f,void* param){
    assert(*ndim==2 && *ncomp==8);
    f[0] = x[0]*x[1];           // = 0.25
    f[1] = x[0]-x[1];           // = 0
    f[2] = x[1]+x[0];           // = 1
    f[3] = std::pow(x[0],x[1]); // = ln(2)
    f[4] = 1;                   // = 1
    f[5] = x[0]*x[0]*x[0];      // = 0.25
    f[6] = x[0]*x[0]*x[0]*x[1]; // = 0.125
    f[7] = sin(x[0])*cos(x[1]); // = sin(1) * ( 1 - cos(1))
    return 0; // succes
}
/** the analytical result of the integration of myintegrand_multicomp **/
int myintegrand_multicomp_integral(const int* ncomp, double* f,void* param){
    assert(*ncomp==8);
    f[0] = 0.25;
    f[1] = 0;
    f[2] = 1;
    f[3] = log(2);
    f[4] = 1;
    f[5] = 0.25;
    f[6] = 0.125;
    f[7] = sin(1) * ( 1 - cos(1));
    return 0; // succes
}

bool test_integration(struct Integrator& i,integrand_t myintegrand, int (*solution)(const int*,double*,void*)){
    double f[i.ncomp]; // the result of the integration
    double e[i.ncomp]; // the estimated error of the integration routine
    double a[i.ncomp]; // the analytical result
    i.integrand = myintegrand;
    i.exec(f,e);
    solution(&i.ncomp,a,i.params);
    bool pass = true; // check if integration was succesful
    printf("#[Info] -- component [--] : (    result,     error) : (analytical), [   diff   ], [pass]\n");
    printf("#=======================================================================================\n");
    for (int c=0;c<i.ncomp;c++){
        printf("#[Info] -- component [%2d] : (%10.3e,%10.3e) : (%10.3e), [%10.3e], ",c,f[c],e[c],a[c],fabs(f[c]-a[c]));
        pass = ( pass && fabs(f[c]-a[c]) <= e[c]);
        if ( fabs(f[c] - a[c]) <= e[c])
            printf("[oooo]\n");
        else
            printf("[xxxx]\n");
    }
    return pass;
}

void run_integration_tests(struct Integrator& i,std::string integrator_name,bool multicomp=true){
    struct integrand_params p;
    p.a = 2.0;
    i.ncomp = 1;
    i.ndim  = 1;
    i.params = (void*) &p;
    printf("#[Info] Integrator: %s\n",integrator_name.c_str());
    printf("\n#[Info] Testing myintegrand0...\n");
    test_integration(i,myintegrand0,myintegrand0_integral);
    printf("\n#[Info] Testing myintegrand1...\n");
    test_integration(i,myintegrand1,myintegrand1_integral);
    if (multicomp){
        i.ncomp = 8;
        i.ndim  = 2;
        printf("\n#[Info] Testing myintegrand_multicomp...\n");
        test_integration(i,myintegrand_multicomp,myintegrand_multicomp_integral);
    }
    printf("\n");
}

void test_cuba(){
    Cuba_common c;
    c.nvec    = 1;
    c.epsabs  = 0;
    c.epsrel  = 1e-4;
    c.flags   = 0x00;
    c.mineval = 10000;
    c.maxeval = 1000000;
    c.statefile = nullptr;
    c.spin      = nullptr;
    
    Integrator_vegas vi(c);
    vi.seed   = 1234;
    vi.nstart = 1000;
    vi.nincrease = 100;
    vi.nbatch = 100;
    vi.gridno = 0;
    
    Integrator_cuhre ci(c);
    ci.key = 0;
    
    run_integration_tests(vi,"Vegas");
    run_integration_tests(ci,"Cuhre");
}
void test_gsl(){
    Integrator_gsl_monte_plain p;     
    p.calls = 100000; // max number of evals
    run_integration_tests(p,"gsl_monte_plain",false);

    Integrator_gsl_monte_vegas v;
    v.warm_up_calls = 1000;
    v.in_iteration_calls = 10000;
    v.delta_chi_squared = 0.5;
    run_integration_tests(v,"gsl_monte_vegas",false);

    Integrator_gsl_monte_miser m;
    m.calls = 100000;
    run_integration_tests(m,"gsl_monte_miser",false);
}


void test_cubature(){
    Integrator_hcubature hc;
    hc.maxEval = 100000;
    hc.epsabs  = 0.;
    hc.epsrel  = 1e-4;
    hc.err_norm= ERROR_INDIVIDUAL;
    
    run_integration_tests(hc,"hcubature");

    Integrator_pcubature pc;
    pc.maxEval = 100000;
    pc.epsabs  = 0.;
    pc.epsrel  = 1e-4;
    pc.err_norm= ERROR_INDIVIDUAL;
    
    run_integration_tests(pc,"pcubature");
}

int main(){
    printf("#######################\n");
    printf("#[Info][TESTING] : cuba\n"); 
    printf("#######################\n\n");
    test_cuba();
    printf("#######################\n");
    printf("#[Info][TESTING] : gsl \n"); 
    printf("#######################\n\n");
    test_gsl();
    printf("###########################\n");
    printf("#[Info][TESTING] : cubature \n"); 
    printf("###########################\n\n");
    test_cubature();
    return 0;
}
