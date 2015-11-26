#include <math.h>
#include <cassert>
#include <cstdio>

struct genz_params{
    double* c = nullptr;
    double* w = nullptr;
};

int genz_oscillatory(const int* ndim,const double* x,const int* ncomp,double* f,void* param){
    assert(*ncomp==1);
    struct genz_params gp = *(struct genz_params*) param;
    double p=0;
    for( int i=0; i<*ndim; i++){
        p+=gp.c[i]*x[i];
    }
    f[0] = cos(p + 2*M_PI*gp.w[0]);
    return 0; // succes!
}

int genz_productPeak(const int* ndim,const double* x,const int* ncomp,double* f,void* param){
    assert(*ncomp==1);
    struct genz_params gp = *(struct genz_params*) param;
    f[0]=1.; // prepare on 1 for *= operations
    for (int i=0; i<*ndim; i++){
        double diff = x[i]-gp.w[i];
        f[0] *= 1./(diff*diff + 1./(gp.c[i]*gp.c[i]));
    }
    return 0; // succes!
}

int genz_cornerPeak(const int* ndim,const double* x,const int* ncomp,double* f,void* param){
    assert(*ncomp==1);
    struct genz_params gp = *(struct genz_params*) param;
    double p=0;
    for( int i=0; i<*ndim; i++){
        p+=gp.c[i]*x[i];
    }
    f[0] = 1./pow(1.+p,*ndim+1);
    return 0; // succes!
}

int genz_Gaussian(const int* ndim,const double* x,const int* ncomp, double* f, void* param){
    assert(*ncomp==1);
    struct genz_params gp = *(struct genz_params*) param;
    double p=0.;
    for (int i=0;i<*ndim;i++){
        p+=gp.c[i]*gp.c[i]*( x[i] - gp.w[i] )*( x[i] - gp.w[i]);
    }
    f[0] = exp(-p);
    return 0; // succes!
}

int genz_c0Continuous(const int* ndim,const double* x, const int* ncomp, double* f,void* param){
    assert(*ncomp==1);
    struct genz_params gp = *(struct genz_params*) param;
    double p=0.;
    for (int i=0;i<*ndim;i++){
        p+=gp.c[i]*fabs( x[i] - gp.w[i] );
    }
    f[0] = exp(-p);
    return 0; // succes!
}

int genz_disContinuous(const int *ndim,const double* x,const int* ncomp,double* f,void* param){
    assert(*ncomp==1 && *ndim >= 2);
    struct genz_params gp = *(struct genz_params*) param;
    if ( x[0] < gp.w[0] || x[1] < gp.w[1]){
        double p=0.;
        for (int i=0;i<*ndim;i++){
            p += gp.c[i]*x[i];
        }
        f[0] = exp(p);
    } else {
        f[0] = 0.;
    }
    return 0; // succes!
}


    
