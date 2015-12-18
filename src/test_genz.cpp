#include "test_integrands.hpp"
#include "cuba_wrappers.hpp"
#include "cubature_wrappers.hpp"
#include "gsl_wrappers.hpp"
#include "integrator.hpp"
#include <cstdio>

void integrate(Integrator&,integrand_t);

void test_cuba(integrand_t myIntegrand){
    Cuba_common c;
    c.nvec    = 1;
    c.epsabs  = 1e-5;
    c.epsrel  = 1e-5;
    c.flags   = 0x00;
    c.mineval = 10000;
    c.maxeval = 100000;
    c.statefile = nullptr;
    c.spin      = nullptr;
    
    Integrator_vegas vi(c);
    vi.seed   = 1234;
    vi.nstart = 1000;
    vi.nincrease = 100;
    vi.nbatch = 100;
    vi.gridno = 0;

    Integrator_cuhre ci(c);
    ci.key = 9;
    
    Integrator_divonne di(c);
    di.key1 = 9;
    di.key2 = 9;
    di.key3 = 1;
    di.maxpass = 5;
    di.ngiven  = 0;
    di.ldxgiven= 0;
    di.nextra  = 0;
    di.seed    = 1234;
    di.border  = 0;
    di.maxchisq= 5;
    di.mindeviation = .25;
    di.xgiven = nullptr;
    di.peakfinder = nullptr;

    Integrator_suave si(c);
    si.nnew = 1000;
    si.nmin = 100;
    si.seed = 1234;
    si.flatness = 1.;
    
    printf("%-20s","Vegas   :");
    integrate(vi,myIntegrand);
    printf("%-20s","Cuba    :");
    integrate(ci,myIntegrand);
    printf("%-20s","Suave   :");
    integrate(di,myIntegrand);
    printf("%-20s","Divonne :");
    integrate(si,myIntegrand);
    //Integrator_divonne di(c);

}

void test_cubature(integrand_t myIntegrand){
    Integrator_hcubature hi;
    hi.maxEval = 10000;
    hi.epsabs  = 0.;
    hi.epsrel  = 1e-5;
    
    Integrator_pcubature pi;
    pi.maxEval = 10000;
    pi.epsabs  = 0.;
    pi.epsrel  = 1e-5;

    printf("%-20s","hCubature :");
    integrate(hi,myIntegrand);
    printf("%-20s","pCubature :");
    integrate(pi,myIntegrand);
}

void test_gsl(integrand_t myIntegrand){
    Integrator_gsl_monte_plain p;     
    p.calls = 10000; // max number of evals

    Integrator_gsl_monte_vegas v;
    v.calls = 10000;

    Integrator_gsl_monte_miser m;
    m.calls = 10000;
    
    printf("%-20s","gsl_monte_plain : ");
    integrate(p,myIntegrand);
    printf("%-20s","gsl_monte_vegas : ");
    integrate(v,myIntegrand);
    printf("%-20s","gsl_monte_miser : ");
    integrate(m,myIntegrand);
}

void integrate(Integrator &i,integrand_t myIntegrand){ // <-- this function would be burried deep in your code //
    double f,e;
    i.ndim = 3;
    i.ncomp = 1;
    i.integrand = myIntegrand;
    struct genz_params gp;
    /** see cuba.pdf (of the cuba library for more explanation of w,c parameters **/
    double c[] = {.1,.1,.2}; // sum_i abs( c[i] ) determines difficulty of integrand, smaller -> more difficult
    double w[] = {.5,.5,.7}; // should not influence difficulty, changes location of peaks...
    gp.c = c;
    gp.w = w;
    i.params = (void*) &gp;
    i.exec(&f,&e);
    printf("Result, Error: %e    %e\n",f,e);
}

int main(){
    printf("\n[Info] Integrand = genz_oscillatory\n\n");
    test_cuba(&genz_oscillatory);
    test_cubature(&genz_oscillatory);
    test_gsl(&genz_oscillatory);
    printf("\n[Info] Integrand = genz_productPeak\n\n");
    test_cuba(&genz_productPeak);
    test_cubature(&genz_productPeak);
    test_gsl(&genz_productPeak);
    printf("\n[Info] Integrand = genz_cornerPeak\n\n");
    test_cuba(&genz_cornerPeak);
    test_cubature(&genz_cornerPeak);
    test_gsl(&genz_cornerPeak);
    printf("\n[Info] Integrand = genz_c0Continuous\n\n");
    test_cuba(&genz_c0Continuous);
    test_cubature(&genz_c0Continuous);
    test_gsl(&genz_c0Continuous);
    printf("\n[Info] Integrand = genz_disContinuous\n\n");
    test_cuba(&genz_disContinuous);
    test_cubature(&genz_disContinuous);
    test_gsl(&genz_disContinuous);
    return 0;
}
