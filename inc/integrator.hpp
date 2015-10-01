#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

/* this is the same typedef as in the cuba libraries...
 * no proz if you try to typedef twice exactly the same thing
 **/
typedef int (*integrand_t)(const int*,const double*,const int*,double*,void*);


/** All integrators should inherit from this
  * base class. It has a very simple structure
  * reflecting what all integrators should look like.
  * Integrator interfaces should only expose:
  *   - the integration result ( double* f)
  *   - the integration error  ( double* e)
  * 
  * Only thing that the integrator should know is the
  * integrand to integrate. integrand_t integrand
  * and the number of dimensions/components.
  *
  * All the rest (integrator specific settings) like
  * error tolerance, min/max number of evaluations, etc
  * should be set externally.
  *
  *
  * So a typical example might look like this. Deep inside your
  * code you want to perform some integration.
  * 
  * doTheIntegration(...,Integrator& i,...) { 
  *    ....
  *    integrand_t myIntegrand;
  *    struct myLoadOfParametersForMyIntegrand;
  *    double result[ncomp],error[ncomp];
  *    
  *    i.ndim      = nofDimensions;
  *    i.ncomp     = nofComponents;
  *    i.integrand = myIntegrand;
  *    i.params    = (void*) myLoadOfParametersForMyIntegrand;
  *    i.exec(&f,&e);
  * }
  * ...
  * int myIntegrand(const int*,const double* x,const int*,double* f, void* p){
  * {
  *     struct myLoadOfParameters param = *(struct myLoadOfParameters*) p;
  *     f[0] = p.a*exp( -p.b*cos(x[0]));
  *     f[1] = ...
  *     return succes; // succes is 0 if integrand evaluated succesfully
  * }
  *
  * Where Integrator is constructed externally! So your code should
  * be very flexible in changing integration routines!
  * 
  * external code might look something like:
  *
  * int main() {
  *    ....
  *    myComplicatedClass myclass(a,b,c,....);
  *    
  *    myCustomIntegrator v;
  *    v.maxeval = 10000000;
  *    v.mineval = 1999;
  *    v.epsabs  = 1e-5;
  *    ....
  *    myclass.doTheIntegration(...,&v,...);
  *    
  *    return 0;
  * }
  * 
  * Of course myCostumIntegrator has to inherit from Integrator...
  *
  * Note that no integration boundaries are passed. In this I follow
  * the filosofy of the Cuba integration libraries. All integrations 
  * are performed in the unit hypercube and a jacobian is introduced to 
  * compensate for this. If you want explicit boundaries you can always
  * sneak them in into the void* params field where your integrand parameters go.
  * 
  *
  * 
  */
struct Integrator {
    virtual int exec(double* f,double* e) = 0;
    integrand_t integrand;
    int ndim,ncomp;
    void* params;
};


#endif // INTEGRATOR_HPP
