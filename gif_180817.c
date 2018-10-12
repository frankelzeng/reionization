#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/*Delta N_HI*/
#define DNHI 2.5e16
/*N_grid*/
#define NGRID 2000
/*Delta t^prime*/
#define DTIMESTEP 2.5e15
/*N t^prime ???*/
#define NTIMESTEP 12000

/*Number of Frequency bins*/
#define N_NU 128
/*Rydberg temperature*/
#define RYD_K 157800
#define ION_HE 1.809
/*f_He*/
#define ABUND_HE 0.079

/* Sets blackbody incident spectrum at given T */
void set_bb(double *frac, double T) {
  int i;
  double nu, tot=0.;
  double myxs = 0, epsilon = 0;

  T /= RYD_K;
  for(i=0;i<N_NU;i++) {
    nu = pow(4., (i+0.5)/N_NU);
#ifdef THIS_HCOL
    epsilon = sqrt(nu-1.);
    myxs = THIS_HCOL * 6.30e-18*pow(nu,-4.)*exp(4.-4.*atan2(epsilon,1.)/epsilon)/(1.-exp(-2*M_PI/epsilon));
#endif
    tot += frac[i] = nu*nu*nu/(exp(nu/T)-1.) * exp(-myxs);
  }
  for(i=0;i<N_NU;i++) frac[i]/=tot;
}

/* Sets cross sections [in cm^2] */
void set_sigma(double *nu, double *sigH, double *sigHe) {
  int i;
  double epsilon, x, y;

  for(i=0;i<N_NU;i++) {
    nu[i] = pow(4., (i+0.5)/N_NU);
    epsilon = sqrt(nu[i]-1.);
    sigH[i] = 6.30e-18*pow(nu[i],-4.)*exp(4.-4.*atan2(epsilon,1.)/epsilon)/(1.-exp(-2*M_PI/epsilon));
    sigHe[i] = 0.;
    if (nu[i]>ION_HE) {

      /* Verner et al. 1996
       * sigma = sigma_0 [ (x-1)^2 + y_w^2 ] y^{0.5P-5.5} ( 1 + sqrt{y/y_a} )^{-P}
       * y = sqrt(x^2 + y_1^2)
       * x = E/E_0-y_0
       */
      x = nu[i]/1.0007-0.4434;
      y = sqrt(x*x+2.136*2.136);
      sigHe[i] = 9.492e-16 * ( (x-1.)*(x-1.)+2.039*2.039 ) * pow(y,0.5*3.188-5.5) * pow(1.+sqrt(y/1.469), -3.188);
    }
#if 0
    sigH[i] = 6.2e-18*pow(nu[i],-3.);
    sigHe[i] = 0.;
    if (nu[i]>ION_HE) sigHe[i] = 7.9e-18*pow(nu[i]/ION_HE,-3.);
#endif
  }
}

/* Takes in neutral fractions y1H[] and y1He[], and the incident photon spectrum.
 * Outputs rates: delta(y1H), delta(y1He), delta(EH) [Ryd/H]
 *  are done with a "time variable" that is 1 photon/cm^2
 *  [i.e. scaled time = F * t', where F = ionizing photon flux in ph/cm^2/s]
 */
void get_ion_rate(double *y1H, double *y1He, double *fracflux, double *dy1H, double *dy1He, double *dEH) {
  static int is_initialized = 0;
  static double *nu, *sigH, *sigHe;
  long i,j;
  double tauH, tauHe, flux, tautot, wt;

  /* Build internal tables */
  if (!is_initialized) {
    is_initialized = 1;
    sigH = (double*)malloc((size_t)(3*N_NU*sizeof(double)));
    sigHe = sigH + N_NU;
    nu = sigHe + N_NU;
    set_sigma(nu,sigH,sigHe);
  }

  for(j=0;j<NGRID;j++) dy1H[j] = dy1He[j] = dEH[j] = 0.;

  /* Now get contribution from each bin */
  for(i=0;i<N_NU;i++) {
    flux = fracflux[i];
    for(j=0;j<NGRID;j++) {
      /* Optical depth in this slice */
      tauH = DNHI * sigH[i] * y1H[j];
      tauHe = ABUND_HE * DNHI * sigHe[i] * y1He[j];
      tautot = tauH + tauHe;
      /* weight = (mean flux in this slice) */
      wt = flux * (tautot>1e-5? (1.-exp(-tautot))/tautot: 1.-tautot/2.);

      dy1H[j] -= wt * sigH[i] * y1H[j];
      dy1He[j] -= wt * sigHe[i] * y1He[j];
      dEH[j] += wt * sigH[i] * y1H[j] * (nu[i]-1.);
      dEH[j] += wt * sigHe[i] * y1He[j] * (nu[i]-ION_HE) * ABUND_HE;

      flux *= exp(-tautot);
    }
  }
}

/* Cooling rate [divided by n_H] in Ry cm^3/s at given neutral fractions
 */
double get_cooling_rate(double Te, double y1H, double y1He) {

  double q12, q13, gamma;

  /* 1->2 */
  gamma = 0.531+2.71e-5*Te-3.22e-10*Te*Te+3.88e-16*Te*Te*Te;
  if (Te>2.5e4)
    gamma = 0.637+1.47e-5*Te-5.92e-12*Te*Te-3.78e-18*Te*Te*Te;
  q12 = 8.6287e-6/2./sqrt(Te)*gamma*exp(-0.75*RYD_K/Te);

  /* 1->3 */
  gamma = 0.350-2.62e-7*Te-8.15e-11*Te*Te+6.19e-15*Te*Te*Te;
  if (Te>2.5e4)
    gamma = 0.276+4.99e-6*Te-8.85e-12*Te*Te+7.18e-18*Te*Te*Te;
  q13 = 8.6287e-6/2./sqrt(Te)*gamma*exp(-0.8888888889*RYD_K/Te);

  return( (3./4.*q12 + 8./9.*q13)*y1H*(1.-y1H+ABUND_HE*(1.-y1He)) );
}

int main(int argc, char **argv) {
  int i;
  long j, istep;
  double fracflux[N_NU], nu[N_NU], sigH[N_NU], sigHe[N_NU];
  double y1H[NGRID], y1He[NGRID], EH[NGRID], Te[NGRID];
  double dy1H[NGRID], dy1He[NGRID], dEH[NGRID];
  double T, U;

  sscanf(argv[1], "%lf", &T);
  set_bb(fracflux,T);
  set_sigma(nu,sigH,sigHe);
#if 0
  for(i=0;i<N_NU;i++)
    printf("%2d %8.6lf %8.6lf %11.5lE %11.5lE\n", i, nu[i], fracflux[i], sigH[i], sigHe[i]);
#endif
  i=0; /* to avoid complaining */

  /* Set up initial conditions */
  for(j=0; j<NGRID; j++) {
    y1H[j] = y1He[j] = 1.;
    EH[j] = 1.e-10;
  }

  sscanf(argv[2], "%lf", &U);
  for(istep=0;istep<NTIMESTEP;istep++) {
    if (istep==0)
      for(j=0;j<NGRID;j++)
        Te[j] = EH[j]/1.5/(2.-y1H[j]+ABUND_HE*(2.-y1He[j]))*RYD_K;
    get_ion_rate(y1H,y1He,fracflux, dy1H,dy1He,dEH);
    for(j=0;j<NGRID;j++)
      dEH[j] -= get_cooling_rate(Te[j], y1H[j], y1He[j])/(1.+ABUND_HE)/U;
    for(j=0;j<NGRID;j++) {
      y1H[j] += DTIMESTEP * dy1H[j];
      y1He[j] += DTIMESTEP * dy1He[j];
      EH[j] += DTIMESTEP * dEH[j];
      Te[j] = EH[j]/1.5/(2.-y1H[j]+ABUND_HE*(2.-y1He[j]))*RYD_K;
    }
  }

  for(j=0; j<NGRID; j++) {
    printf("%4ld %11.5lE %8.6lf %8.6lf %8.6lf %7.1lf\n",
      j, (j+.5)*DNHI, y1H[j], y1He[j], EH[j], Te[j]); 
  }

  return(0);
}
