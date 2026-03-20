/* 
*
* This Script contains the UDFs required for simulation of Hydrogen Degassing  using a Vacuum arc degasser.
* The following functions are present
* 1. DEFINE_MASS_TRANSFER(){} - used to model the source term for mass transfer 
* 2. DEFINE_PROPERTY(){} - used to model density of phases- argon and liquid steel.
*
*/


#include "udf.h"


//The following parameters can be changed based on composition

#define f_h 1.09
#define Diff_coef 1e-07
#define rho_l 7000
#define R_idlg 8.314
#define M_ar 0.039948
#define M_hyd 0.002016
#define nu 7.14e-07



DEFINE_MASS_TRANSFER( hydrogen_transport, cell, thread_mix, index_ls, species_h2_ls, index_gas, species_h2_gas){
 Thread *thread_ls = THREAD_SUB_THREAD(thread_mix, index_ls);
 Thread *thread_gas = THREAD_SUB_THREAD(thread_mix, index_gas);
 real xc[ND_ND];
 C_CENTROID(xc,cell,thread_mix);
 real temp = 0;


 /* The following code can be used to limit the mass transfer reaction to regions where VOFs of argon or steel are higher than 
 * limit
 *
 *
 * //if volume fraction of cell for argon is tending to zero, the mass transfer is also 0.
 * //similarly if the volume fraction of liquid steel tends to zero, the mass transfer is 0
 * 
 * if(C_VOF(cell,thread_gas) < 1e-12){
 *    return temp;
 * }
 * if(C_VOF(cell,thread_ls) < 1e-12){
 *   return temp;
 *  }
 */

 /*
 * Due to numerical instability at the start of the simulation, there is a possibility of oscillating pressure.
 * In such cases the pressure can be very high, or negative. As pressure is required to determine the source term, 
 * an unphysical value may cause divergence.Thus we limit our mass transfer to regions where pressure is realistic 
 */
 real pressure_cell;
 pressure_cell = C_P(cell,thread_mix);
 if(pressure_cell > 404650){
   return temp;
 }
 if(pressure_cell < 500){
   return temp;
 }

 //calculate partial pressure of hydrogen in aolecular argon bubbles
 real C_h_arg;
 real C_arg_arg;
 C_h_arg = C_YI(cell,thread_gas,species_h2_gas);
 C_arg_arg = C_YI(cell, thread_gas, 1); //argon species index is 1; ensure in fluent that its the last species
 real denominator;
 denominator = C_h_arg/M_hyd + C_arg_arg/M_ar;
 denominator = MAX(denominator,1e-12);
 real partial_pressure_h2;
 partial_pressure_h2 = pressure_cell*(C_h_arg/M_hyd)/(denominator);
 partial_pressure_h2 = MAX(partial_pressure_h2,0);

 //calculate C_h2_equilibrium
 real C_h2_equilibrium;
 real Coeff_ch2eq;
 Coeff_ch2eq = 0.002246; //exp(-del Go / rt)/fh
 C_h2_equilibrium = Coeff_ch2eq*sqrt(partial_pressure_h2)*0.000001; //convert from wt-ppm, to mass fraction.

 //calculate k
 real k;
 real eps;
 eps = C_D(cell,thread_ls);
 eps = MAX(eps,0);
 real internal = Diff_coef*sqrt(eps/nu);
 internal = MAX(internal,0);
 k = 0.3*sqrt(internal);

 //calculation of Area
 real sauter_mean;
 sauter_mean = MAX(C_PHASE_DIAMETER(cell,thread_gas),0.0001);
 real Area;
 Area = 6*C_VOF(cell,thread_gas)/sauter_mean;

 //calculate the source
 real source_term;
 real C_h_ls;
 C_h_ls = C_YI(cell,thread_ls,species_h2_ls);
 source_term = rho_l*Area*k*( C_h_ls - C_h2_equilibrium);
 return source_term;
}
//density udf
DEFINE_PROPERTY(argon_phase_density,c,t){
  real y_cord;
  real y_max = 2.5;
  real temp = 1871;
  real M_arg = 0.03995;
  real M_h = 0.002016;
  real R_no = 8.314;
  real density;
  real xc[ND_ND];
  real C_arg_arg = C_YI(c,t,0);
  real C_h_arg = C_YI(c,t,1);
  real M_av = (C_h_arg/M_h + C_arg_arg/M_arg);
  C_CENTROID(xc,c,t);
  y_cord = xc[1];
  real h = MAX(y_max - y_cord,0.0);
  density =(1013 + (h)*9.81*7000)/(R_no*temp*M_av);
  return density;
}
DEFINE_PROPERTY(steel_phase_density,c,t){
  real density = 7000;
  return density;
}
