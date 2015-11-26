#include "test_integrands.hpp"
#include "cuba_wrappers.hpp"
#include "integrator.hpp"
#include <cstdio>

void integrate(Integrator&,integrand_t);

void test_cuba(integrand_t integrand){
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
     
    integrate(vi,integrand);
    integrate(ci,integrand);
    integrate(di,integrand);
    integrate(si,integrand);
    //Integrator_divonne di(c);

    delete c.prob;
}

void integrate(Integrator &i,integrand_t integrand){ // <-- this function would be burried deep in your code //
    double f,e;
    i.ndim = 3;
    i.ncomp = 1;
    i.integrand = integrand;
    struct genz_params gp;
    double c[] = {1.,1.,1.};
    double w[] = {1.,1.,1.};
    gp.c = c;
    gp.w = w;
    i.params = (void*) &gp;
    i.exec(&f,&e);
    printf("Result, Error: %e    %e\n",f,e);
}

int main(){
    printf("\n[Info] Integrand = genz_oscillatory\n\n");
    test_cuba(&genz_oscillatory);
    printf("\n[Info] Integrand = genz_productPeak\n\n");
    test_cuba(&genz_productPeak);
    printf("\n[Info] Integrand = genz_cornerPeak\n\n");
    test_cuba(&genz_cornerPeak);
    printf("\n[Info] Integrand = genz_c0Continuous\n\n");
    test_cuba(&genz_c0Continuous);
    printf("\n[Info] Integrand = genz_disContinuous\n\n");
    test_cuba(&genz_disContinuous);
    return 0;
}
