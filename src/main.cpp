

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include "omp.h"
#include <math.h>
#include <random>
#include <iostream>
#include <cmath>

#include <stdio.h>

#include <string>
#include <iostream>
#include <sstream>
#include <vector>


#include "c_headers.h"
#include "cpp_headers.h"



#include <sys/resource.h>

#include <omp.h>
#include "omp.h"
int main(int argc, const char **argv)
{
  double start,end;
  double total_start = omp_get_wtime();
  struct GAS gas;
  struct BOX box;
  struct CELLS *cells;
  struct Collision CBA;
  struct NerualNetwork NN;
  struct GPR gpr;

  FILE *fp;
//   double mu_crr;
//   double mu0;
//   double a0;
  int i,j,k, ii, iii;
  if (argc > 1){
    box.file_name = argv[1];
  } else {
    std::cout << " problem reading comand line \n";
    return 1;
  };


//solver( );
  ////					 reading input script
  readSettingsFile(&gas, &box);
// reading NN
//  if(gas.model == "Hybrid"){
double yy[3];
double yy_NNpred[3] = {-4.3389455617270652e-02,  5.0189721301572154e-01,  1.4463151872451448e-02};
double xx = -8.7450418019444184e-02;
    read_NN(&NN);
    evaluation(&NN, xx, yy);
    if( abs(yy[0]-yy_NNpred[0])+abs(yy[1]-yy_NNpred[1])+abs(yy[2]-yy_NNpred[2])<1e-6 ){
      printf("The NN is set properly\n");
    }
    read_GPR(&gpr);
    evaluation_GPR(&gpr, xx, yy);
    if( abs(yy[0]-yy_NNpred[0])+abs(yy[1]-yy_NNpred[1])+abs(yy[2]-yy_NNpred[2])<1e-6 ){
      printf("The GPR is set properly\n");
    }
//  }
// Memory Allocation
if(box.problem == dcones){
  cells = (struct CELLS *) malloc( box.Nx[0]*box.Nx[1]*box.Nx[2]*box.Ny * sizeof(struct CELLS) );
  make_cells(&box, cells, &gas);
}
else if(box.problem == flatnose){
  cells = (struct CELLS *) malloc( (box.N1[0]+box.N1[1])*(box.N2[0]+box.N2[1]) * sizeof(struct CELLS) );
  make_cells(&box, cells, &gas);
}
else{
  cells = (struct CELLS *) malloc( box.N[0]*box.N[1]*box.N[2] * sizeof(struct CELLS) );
  make_cells(&box, cells, &gas);
  for (int j=0; j < box.N[0]; j++)
    for (int k=0; k < box.N[1]; k++)
        for (int l=0; l < box.N[2]; l++){
	  int num = l*box.N[0]*box.N[1] + k*box.N[0] +j;
	  cells[num].indices_inside = (int *) malloc( 1 * sizeof(int) );
	  for(int i=0; i<1; i++){
	    cells[num].indices_inside[i] = -1;}
	}
}
////////			allocating memory for variables
  double *U1, *U2, *U3, *x1, *x2;
  int *index, *flag, done;
  double *F1, *F2, *F3, *F1_old, *F2_old, *F3_old, *x3, *x3_old;
  double *xi_1f, *xi_2f, *xi_3f, *xi_1, *xi_2, *xi_3;
  double *V, *V0, *Ek, *Ek0;
  double *Mp1, *Mp2, *Mp3;
  int *color;
  int *n_ratio;
  double *T, *rho;
  U1 = (double *) malloc( gas.N * sizeof(double) );
  U2 = (double *) malloc( gas.N * sizeof(double) );
  U3 = (double *) malloc( gas.N * sizeof(double) );
  x1 = (double *) malloc( gas.N * sizeof(double) );
  x2 = (double *) malloc( gas.N * sizeof(double) );
  index = (int *) malloc( gas.N * sizeof(int) );
  flag = (int *) malloc( gas.N * sizeof(int) );
  n_ratio = (int *) malloc( gas.N * sizeof(int) );
  color = (int *) malloc( gas.N * sizeof(int) );
if(box.problem != dcones && box.problem != flatnose){
 	F1 = (double *) malloc( gas.N * sizeof(double) );
 	F2 = (double *) malloc( gas.N * sizeof(double) );
 	F3 = (double *) malloc( gas.N * sizeof(double) );

 	Mp1 = (double *) malloc( gas.N * sizeof(double) );
 	Mp2 = (double *) malloc( gas.N * sizeof(double) );
 	Mp3 = (double *) malloc( gas.N * sizeof(double) );

 	F1_old = (double *) malloc( gas.N * sizeof(double) );
 	F2_old = (double *) malloc( gas.N * sizeof(double) );
 	F3_old = (double *) malloc( gas.N * sizeof(double) );

	x3 = (double *) malloc( gas.N * sizeof(double) );
	x3_old = (double *) malloc( gas.N * sizeof(double) );

	xi_1f = (double *) malloc( gas.N * sizeof(double) );
	xi_2f = (double *) malloc( gas.N * sizeof(double) );
	xi_3f = (double *) malloc( gas.N * sizeof(double) );
	xi_1 = (double *) malloc( gas.N * sizeof(double) );
	xi_2 = (double *) malloc( gas.N * sizeof(double) );
	xi_3 = (double *) malloc( gas.N * sizeof(double) );


	V = (double *) malloc( (box.N[1]+2) * sizeof(double) );
	V0 = (double *) malloc( (box.N[1]+2) * sizeof(double) );
	Ek = (double *) malloc( (box.N[1]) * sizeof(double) );
	Ek0 = (double *) malloc( (box.N[1]) * sizeof(double) );
}

//if( gas.model == "SPH"){
  T = (double *) malloc( gas.N * sizeof(double) );
  rho = (double *) malloc( gas.N * sizeof(double) );
  //for(i=0; i<gas.N; i++){
  //  if(i<gas.Nv1)
  //    T[i] = gas.T1;
  //  else
  //    T[i] = gas.T3;
  //}
//}
 if(gas.model == "MD" ){
      for(int id=0; id<gas.N; id++){
             F1_old[id] = 0.0;
             F2_old[id] = 0.0;
             F3_old[id] = 0.0;
             F1[id] = 0.0;
             F2[id] = 0.0;
             F3[id] = 0.0;
      }
 }


  double * x2_old = (double *) malloc( gas.N * sizeof(double) );
  double * x1_old = (double *) malloc( gas.N * sizeof(double) );

/////////////////////////////////////////////////////////////////////


double pi = acos(-1);
double X = increase_collision_rate(gas.n, &gas);

  double p_aim = gas.n*gas.kb*gas.T*(1.0 + X*gas.b*gas.n);
  gas.p =  p_aim;

  printf("nkT(1+nbY)-patt = %e\n",p_aim-2.0*gas.phi*gas.n*gas.n*gas.b);
  printf("nkT(1+Xnb) = %e\n",p_aim);
  printf("nkT(1+nb) = %e\n",gas.n*gas.kb*gas.T*(1.0 + gas.b*gas.n));
  printf("nkT = %e\n",gas.n*gas.kb*gas.T);
  printf("nb = %e\n",gas.n*gas.b);
  printf("b = %e\n",gas.b);
  printf("n sigma^3 = %e\n",gas.n*gas.sigma*gas.sigma*gas.sigma);
 for( i=0; i<gas.N; i++){
	if(box.problem==dcones || box.problem == flatnose){
		flag[i] = 1;
		// in dcone, flag means weather the particle is in or out
	}
	else{
		flag[i] = 0;
	}
	color[i] = 0;// color is used for computing evaporation coefficient
 }


 double a, b, c, y, x;
  /// realization
for(int realization = 1; realization <= box.num_realizations; realization++){

//if(box.reset == 1)
//	read_pos_vel(x1, x2, U1, U2, U3, gas.N);
//else
initialization(U1, U2, U3, x1, x2, x3, &gas, &box, cells);

//shuffle(x1, x2, x3, &gas);
//BC_specular_reflection(U1, U2, U3, x1, x2,x3, &gas, &box);
     double std = sqrt(gas.kb*gas.T/gas.m);
     double mean;
  std::cout<<"std of U suppoesed to be = "<<std<<"\n";
  std::cout<<"std of U3 = "<<standard_deviation(U1,gas.N,&mean)<<"\n";

  printf("T0 = %e\n", gas.m*standard_deviation(U1, gas.N, &mean)*standard_deviation(U1, gas.N, &mean)/gas.kb);
  printf("T0 = %e\n", gas.m*standard_deviation(U2, gas.N, &mean)*standard_deviation(U2, gas.N, &mean)/gas.kb);
  printf("T0 = %e\n", gas.m*standard_deviation(U3, gas.N, &mean)*standard_deviation(U3, gas.N, &mean)/gas.kb);
  printf("gas.T0 = %e\n", gas.T0);
   printf("gas.T = %e\n", gas.T);
  printf("gas.crref = %e\n", gas.crref);
  int nsamples_after_moving = 1;//floor( 0.05*box.reset_every );
  printf("\n\n -PostProc:\nevery %d\nafter %d\n%d steps\nmove_every %d\nsamples_after_moving %d\n",box.every, box.after, box.num_steps,  box.reset_every, nsamples_after_moving);

  double thermal;
  thermal =  sqrt(gas.kb*gas.T/gas.m);
  double t, total_time;// time_till_now = 0.0;
  int step;
  int N_cell;
  int num;

  step = 0;
  box.step = step;
  t = 0.0;
  gas.times = 1.0;
  gas.N_abs = 0.0;
  if( box.problem == shock )
    cell_update_hybrid(x1, x2, U1, U2, U3, &gas, cells, &box, index, T, rho, color);
  else
    cell_update(x1, x2, U1, U2, U3, &gas, cells, &box, index, &T, rho, color);

  // initilize Fn for particles
  for (i=0; i<gas.N; i++){
	num = index[i];
	n_ratio[i] = cells[num].n_ratio;
  }
if(box.problem == dcones)
	N_cell = (box.Nx[0]+box.Nx[1]+box.Nx[2])*box.Ny;
else if(box.problem == flatnose)
	N_cell = (box.N1[0]+box.N1[1])*(box.N2[0]+box.N2[1]);
else
	N_cell = box.N[0]*box.N[1]*box.N[2];
  // set <M> and <MM> to zero
  for( num=0; num<N_cell; num++){
	//cells[num].T=0.0;
  thermal =  sqrt(gas.kb*cells[num].T/gas.m);
 	cells[num].crm=10.0*( (thermal)) ;//*pow(gas.crref/thermal,2.0*gas.visp-1.0) );
  cells[num].omega_max = 4.0*pi*gas.sigma*gas.sigma*cells[num].crm*cells[num].n*gas.delta_t;
	cells[num].U_space[0]=0.0;
	cells[num].U_space[1]=0.0;
	cells[num].U_space[2]=0.0;
	for(ii=0; ii<9; ii++)
	  cells[num].alpha[ii] = 0.0;
	for(ii=0; ii<3; ii++)
	  cells[num].beta[ii] = 0.0;
 	for(ii=0; ii<6; ii++)
	  cells[num].sum_MM[ii] = 0.0;
 	for(ii=0; ii<3; ii++)
	  cells[num].sum_M[ii] = 0.0;
        cells[num].sum_weight = 0.0;
  	if(gas.model == "ESMC")
        	cells[num].omega_max =cells[num].omega_max*increase_collision_rate(gas.n, &gas);
 }


for(i=0; i<gas.N; i++)
  index[i] = 0;
cells[0].num_inside = gas.N;
  //generating_Python_paraview_script(box.num_steps, &box);
  total_time = box.num_steps*gas.delta_t;

  std::cout << "====== Start of Simulation =============================================================" << std::endl;

  box.dummy_nom = 0.0;//, dummy_denom = 0.0;
  box.area = 0.0;
  if(box.direction[0]==1)
      for(i=0; i<2; i++)
            box.area += box.Area[i];
  if(box.direction[1]==1)
      for(i=2; i<4; i++)
            box.area += box.Area[i];
   if(box.direction[2]==1)
       for(i=4; i<6; i++)
            box.area += box.Area[i];
done = 0;
if(box.reset == 1)
	done = 1;
double tau;
double st0 = 0.0;

/////////////////////
/////////////////////
/////////////////////
double F = 0.0, G=0.0;
int id, id2;
double r, epsilon;
//if( box.direction[0] == 1 ){
	for(i=0; i<gas.N; i++)
    		x1_old[i] = x1[i];
//}
//if( box.direction[1] == 1 ){
	for(i=0; i<gas.N; i++)
    		x2_old[i] = x2[i];
//}
if( box.direction[2] == 1 && box.problem != dcones && box.problem != flatnose){
	for(i=0; i<gas.N; i++)
    		x3_old[i] = x3[i];
}
double sqrt_fr1fr2;
int N_tot_genenrate, N_in;
FILE *gg;
      gg = fopen("relaxation.txt", "w");
      fprintf(gg, "%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s\n",
      "coll",
      "sum_weight",
      "M[0]",
      "M[1]",
      "M[2]",
      "PIJ[0]",
      "PIJ[1]",
      "PIJ[2]",
      "PIJ[3]",
      "PIJ[4]",
      "PIJ[5]");
fclose(gg);

//calculate_potential(V, &gas, &box, cells);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


              /////////////////////////////////////////////////////////////////////////////////
      //##########################       Boundary Condition     ##########################
          /////////////////////////////////////////////////////////////////////////////////
          start = omp_get_wtime();
          if(box.problem == closed_box || box.problem == relaxation)
            BC_specular_reflection(U1, U2, U3, x1, x2,x3, &gas, &box);
            //BC_periodic(U1, U2, U3, x1, x2,x3, &gas, &box);
          else if (box.problem == couette_flow)
            BC_thermal_wall_3D(U1, U2, U3, x1, x2, x3, &gas, &box, index, cells, flag, x1_old, x2_old, x3_old);
            //BC_thermal_wall( U1, U2, U3, x1, x2, &gas, &box, index, cells, flag, x1_old, x2_old);
	  else if ( box.problem == vacuum){
		//BC_specular_reflection(U1, U2, U3, x1, x2,x3, &gas, &box);
		//BC_thermal_wall( U1, U2, U3, x1, x2, &gas, &box, index, cells, flag, x1_old, x2_old);
		BC_evaporation(U1,U2,U3,x1,x2,x3,&gas,&box,index,cells,flag,x1_old,x2_old,x3_old);
	  }
	else if (box.problem == evaporation ){
		//BC_specular_reflection(U1, U2, U3, x1, x2,x3, &gas, &box);
		//BC_thermal_wall( U1, U2, U3, x1, x2, &gas, &box, index, cells, flag, x1_old, x2_old);
		BC_evaporation2(U1,U2,U3,x1,x2,x3,&gas,&box,index,cells,x1_old,x2_old,x3_old);
	  }
	  else if(box.problem == inverted){
		BC_inverted(U1,U2,U3,x1,x2,x3,&gas,&box,index,cells,x1_old,x2_old,x3_old);
	  }

    if( box.problem == shock )
      cell_update_hybrid(x1, x2, U1, U2, U3, &gas, cells, &box, index, T, rho, color);
    else
      cell_update(x1, x2, U1, U2, U3, &gas, cells, &box, index, &T, rho, color);

//writting_position_to_file2D(x1, x2, gas.N, &t, 1, step);

//	post processing: console + writing to file
post_processing(step, cells,  &box, &done, &gas, U1, U2, U3, x1, x2, flag, index, color,  n_ratio, x2_old);
write_post_processing(step, cells,  &box, &done, &gas, U1, U2, U3, x1, x2, flag, index, color,  n_ratio, x2_old);

step = 1;
box.step = step;

double *x1n, *x2n, *U1n, *U2n, *U3n, *x1_oldn, *x2_oldn, *dtn;
int *flagn, *indexn;
double *x1new, *x2new, *U1new, *U2new, *U3new, *x1_oldnew, *x2_oldnew;
int *flagnew, *indexnew;
int done = 1;
//post_vtk(step, cells,  &box, &done, &gas, U1, U2, U3, x1, x2, flag, index);
int post;
int sample = 0;
while(t<total_time){
	st0 = omp_get_wtime();
  post = 0;
	box.step = step;
	if (step % box.every == 0)
		printf("%%%%%%%%%%  Step= %d, Np=%ld   Ex. time till now = %lf min\n",step, gas.N, (omp_get_wtime()-total_start)/60.0);
  if ( box.problem == inverted){
    if ( step % box.reset_every == 0 && step > 1000  ){
        post = reset_inverted_addremove_particles(&U1, &U2, &U3, &x1, &x1_old,&x2, &x2_old, &gas, &box, cells, &index, &flag,&n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3, &color);
        if(post ==1)
          sample = 0;
    }
  }
    if ( box.problem == inverted ){
      reset_inverted_ghost(&U1, &U2, &U3, &x1, &x1_old,&x2, &x2_old, &gas, &box, cells, &index, &flag,&n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3, &color);
    }

   //else{
  //
   //}

	//     store old positions




  if( box.problem == shock )
    cell_update_hybrid(x1, x2, U1, U2, U3, &gas, cells, &box, index, T, rho, color);
  else
    cell_update(x1, x2, U1, U2, U3, &gas, cells, &box, index, &T, rho, color);


  //cell_update(x1, x2, U1, U2, U3, &gas, cells, &box, index, T, rho);
  //update_index_number_inside(x1, x2,  U1, U2, U3, Mp1, Mp2, Mp3, index, cells,  &box, &gas);
  if ( (box.problem == inverted && sample < nsamples_after_moving && step >box.after) || box.problem != inverted ){
    post_processing(step, cells,  &box, &done, &gas, U1, U2, U3, x1, x2, flag, index, color,  n_ratio, x2_old);
    //printf("sample taken for step = %d\n", step);
  }
  write_post_processing(step, cells,  &box, &done, &gas, U1, U2, U3, x1, x2, flag, index, color,  n_ratio, x2_old);

  if ( (box.problem == evaporation || box.problem == vacuum  || box.problem == inverted) && box.step%box.reset_every == 0){
    //if (step % 100 == 0 && step > 10){
    //  reset_inverted_addremove_particles(&U1, &U2, &U3, &x1, &x1_old,&x2, &x2_old, &gas, &box, cells, &index, &flag,&n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3, &color);
    //}
    //ghost_inverted_addremove_particles(&U1, &U2, &U3, &x1, &x1_old,&x2, &x2_old, &gas, &box, cells, &index, &flag,&n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3);
		thermostat( U1, U2, U3, x2, &gas, &box, cells, index);
    cell_update(x1, x2, U1, U2, U3, &gas, cells, &box, index, &T, rho, color);
	}
  //thermostat( U1, U2, U3, x2, &gas, &box, cells, index);
  //cell_update(x1, x2, U1, U2, U3, &gas, cells, &box, index, T, rho);
	//update_index_number_inside(x1, x2, U1, U2, U3, Mp1, Mp2, Mp3, index, cells, &box, &gas);

	//	generate new samples for FP
  if ( gas.model == "FP_dense" || gas.model == "FP_ideal_gas"){
  	random_generator(xi_1, xi_2, xi_3, xi_1f, xi_2f, xi_3f, cells, &box,  &gas, index);
  }
	//	velocity update: collision + long range
	velocity_update(U1, U2, U3, x1, x2, x3, &gas, &box, cells, index, &CBA, xi_1f, xi_2f, xi_3f, Mp1, Mp2, Mp3, F1, F2, F3, rho, T, color);

	// Thermostat
  /*
  if ( gas.model == "ESMC" && (box.problem == evaporation || box.problem == vacuum  || box.problem == inverted)){
    //ghost_inverted_addremove_particles(&U1, &U2, &U3, &x1, &x1_old,&x2, &x2_old, &gas, &box, cells, &index, &flag,&n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3);
    //if (step % 100 == 0 && step > 10){
      //reset_inverted_addremove_particles(&U1, &U2, &U3, &x1, &x1_old,&x2, &x2_old, &gas, &box, cells, &index, &flag,&n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3, &color);
    //}

    cell_update(x1, x2, U1, U2, U3, &gas, cells, &box, index, T, rho);
		thermostat( U1, U2, U3, x2, &gas, &box, cells, index);
	}
  */
	//	cell update
	//cell_update(x1, x2, U1, U2, U3, &gas, cells, &box, index, T, rho);

	//	post processing: console + writing to file
	//post_processing(step, cells,  &box, &done, &gas, U1, U2, U3, x1, x2, flag, index);


	//	apply BC in case of CBA mehtod
	BC_for_CBA(U1, U2, U3, x1, x2, x3, x1_old, x2_old, x3_old, &gas, &box, cells, index, flag);

  //
	//cell_update(x1, x2, U1, U2, U3, &gas, cells, &box, index, T, rho);
	//	stream particles
  if( box.problem == shock && gas.model == "Hybrid" ){
    // for ghost cells, generate SPH/DSMC particles, after Streaming
    // if they're still in the wrong cell, will be deleted!
    add_ghost_particles(&U1, &U2, &U3, &x1, &x1_old, &x2, &x2_old, &gas, &box, cells, &index, &flag,  &n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3, &color,  &T, &rho, &gpr, &NN);
  }
  if ( gas.model == "FP_dense"){
    DFP_extra_streaming(U1, U2, U3, x1, x2, &gas, &box, cells, index, gas.delta_t);
  }
  if( box.problem != shock ){
	 apply_BC(U1, U2, U3, x1, x2, x3, x1_old, x2_old, x3_old, &gas, &box, cells, index, flag, xi_1, xi_2, xi_3, xi_1f, xi_2f, xi_3f, Mp1, Mp2, Mp3, T);
  }
  store_old_positions(x1, x2, x3, x1_old, x2_old, x3_old, &gas, &box);

	stream_particles(U1, U2, U3, x1, x2, x3, x1_old, x2_old, x3_old, &gas, &box, cells, index);
  if( box.problem != shock ){
	 apply_BC(U1, U2, U3, x1, x2, x3, x1_old, x2_old, x3_old, &gas, &box, cells, index, flag, xi_1, xi_2, xi_3, xi_1f, xi_2f, xi_3f, Mp1, Mp2, Mp3, T);
  }
  if ( step > gas.times*box.after){
    measure_pressure(U2, x2, x2_old, &gas, &box, cells);
  }

  if ( step >= box.after ){
    if( box.problem == evaporation )
      compute_evap_coeff_evaporation(U2, x2, x2_old, &gas, &box, color, n_ratio, cells);
  }
  //update_index_number_inside(x1, x2,  U1, U2, U3, Mp1, Mp2, Mp3, index, cells,  &box, &gas);

  if ( (box.problem == inverted && sample < nsamples_after_moving &&  step > box.after) ){
    update_face_info(U1, U2, U3, x1, x2, x1_old, x2_old, &gas, &box, cells, index);
  }
  sample ++;






  //update_face_info(U1, U2, U3, x1, x2, x1_old, x2_old, &gas, &box, cells, index);
	//	apply BC
  if( box.problem != shock ){
	 apply_BC(U1, U2, U3, x1, x2, x3, x1_old, x2_old, x3_old, &gas, &box, cells, index, flag, xi_1, xi_2, xi_3, xi_1f, xi_2f, xi_3f, Mp1, Mp2, Mp3, T);
  }
  else{
    if(gas.model == "SPH"){
      //shock_BC(U2, x2, T, &gas, &box, cells);
      BC_specular_reflection(U1, U2, U3, x1, x2,x3, &gas, &box);
    }
    else{
      //BC_specular_reflection(U1, U2, U3, x1, x2,x3, &gas, &box);
      // add remove for shock BC of Boltzmann-SPH
      addrem_BC_hybrid(&U1, &U2, &U3, &x1, &x1_old, &x2, &x2_old, &gas, &box, cells, &index, &flag,  &n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3, &color,  &T, &rho);
    }
  }
  /*
  if ( gas.model == "FP_dense" && (box.problem == evaporation || box.problem == vacuum  || box.problem == inverted)){
    //ghost_inverted_addremove_particles(&U1, &U2, &U3, &x1, &x1_old,&x2, &x2_old, &gas, &box, cells, &index, &flag,&n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3);
    //if (step % 100 == 0 && step > 10){
      //reset_inverted_addremove_particles(&U1, &U2, &U3, &x1, &x1_old,&x2, &x2_old, &gas, &box, cells, &index, &flag,&n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3, &color);
    //}

    cell_update(x1, x2, U1, U2, U3, &gas, cells, &box, index, T, rho);
		thermostat( U1, U2, U3, x2, &gas, &box, cells, index);
	}
  */
  //if( box.problem == shock )
  //  cell_update_hybrid(x1, x2, U1, U2, U3, &gas, cells, &box, index, T, rho, color);

  if( box.problem == shock && gas.model == "Hybrid" ){
    add_remove_hybrid(&U1, &U2, &U3, &x1, &x1_old, &x2, &x2_old, &gas, &box, cells, &index, &flag,  &n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3, &color,  &T, &rho, &NN, &gpr);
    //hybrid_surface_flux(&U1, &U2, &U3, &x1, &x1_old, &x2, &x2_old, &gas, &box, cells, &index, &flag,  &n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3, &color,  &T, &rho, &gpr);
  }
  if( box.problem == zoomed_inverted){
    zoomed_inverted_addremove_particles(&U1, &U2, &U3, &x1, &x1_old,&x2, &x2_old, &gas, &box, cells, &index, &flag,&n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3);
    update_index_number_inside(x1, x2,  U1, U2, U3, Mp1, Mp2, Mp3, index, cells,  &box, &gas);
  }



	//if ( step >  gas.times*box.after){
	//      if( box.problem == inverted )
	//		compute_evap_coeff(U2, x2, x2_old, &gas, &box, color, n_ratio);
	//}

	// addremove particles
	if( (box.problem == inverted || box.problem == evaporation) && gas.model != "MD" ){
		if(gas.n_ratio > 1)
			inverted_addremove_particles(&U1, &U2, &U3, &x1, &x1_old, &x2, &x2_old, &gas, &box, cells, &index, &flag, &n_ratio, &xi_1, &xi_2, &xi_3, &xi_1f, &xi_2f, &xi_3f, &Mp1, &Mp2, &Mp3, &color);

	}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(box.problem == dcones || box.problem == flatnose){
		start = omp_get_wtime();
		N_tot_genenrate = new_particles(&gas, &box, cells);


		x1n = (double *) malloc( N_tot_genenrate * sizeof(double) );
		x2n = (double *) malloc( N_tot_genenrate * sizeof(double) );
		U1n = (double *) malloc( N_tot_genenrate * sizeof(double) );
		U2n = (double *) malloc( N_tot_genenrate * sizeof(double) );
		U3n = (double *) malloc( N_tot_genenrate * sizeof(double) );
		x1_oldn = (double *) malloc( N_tot_genenrate * sizeof(double) );
		x2_oldn = (double *) malloc( N_tot_genenrate * sizeof(double) );
		dtn = (double *) malloc( N_tot_genenrate * sizeof(double) );
		flagn = (int *) malloc( N_tot_genenrate * sizeof(int) );
		for(i=0; i<N_tot_genenrate; i++)
			flagn[i] = 1;
		generate_new_particles(x1n, x2n, U1n, U2n, U3n, x1_oldn, x2_oldn, &gas, &box, dtn, cells);
		if(box.problem == dcones)
			N_in = dcones_BC(U1, U2, U3, x1, x2, x1_old, x2_old, flag, &gas, &box,  cells, U1n, U2n, U3n, x1n, x2n,x1_oldn, x2_oldn, flagn, dtn);
		else if(box.problem == flatnose)
			N_in = flatnose_BC(U1, U2, U3, x1, x2, x1_old, x2_old, flag, &gas, &box,  cells, U1n, U2n, U3n, x1n, x2n,x1_oldn, x2_oldn, flagn, dtn);
		int N_old, N_new;
/*
		double mn[3], st[3], mnx2n;
		//mn[0] = mean(U1n, gas.n_inflow);
		//mn[1] = mean(U2n, gas.n_inflow);
		//mn[2] = mean(U3n, gas.n_inflow);
		st[0] = gas.m*standard_deviation(U1n, gas.n_inflow, &mn[0])*standard_deviation(U1n, gas.n_inflow, &mn[0])/gas.kb;
		st[1] = gas.m*standard_deviation(U2n, gas.n_inflow, &mn[1])*standard_deviation(U2n, gas.n_inflow, &mn[1])/gas.kb;
		st[2] = gas.m*standard_deviation(U3n, gas.n_inflow, &mn[2])*standard_deviation(U3n, gas.n_inflow, &mn[2])/gas.kb;
*/
		x1new = (double *) malloc( N_in * sizeof(double) );
		x2new = (double *) malloc( N_in * sizeof(double) );
		U1new = (double *) malloc( N_in * sizeof(double) );
		U2new = (double *) malloc( N_in * sizeof(double) );
		U3new = (double *) malloc( N_in * sizeof(double) );
		x1_oldnew = (double *) malloc( N_in * sizeof(double) );
		x2_oldnew = (double *) malloc( N_in * sizeof(double) );
		flagnew = (int *) malloc( N_in * sizeof(int) );
		indexnew = (int *) malloc( N_in * sizeof(int) );

		j=0;
		for(i=0; i<gas.N; i++){
			if(flag[i] == 1){
				x1new[j] = x1[i];
				x2new[j] = x2[i];
				U1new[j] = U1[i];
				U2new[j] = U2[i];
				U3new[j] = U3[i];
				x1_oldnew[j] = x1_old[i];
				x2_oldnew[j] = x2_old[i];
				flagnew[j] = flag[i];
				j++;
			}
		}
		N_old = j;
		for(i=0; i<box.N_tot_genenrate; i++){
			if(flagn[i] == 1){
				x1new[j] = x1n[i];
				x2new[j] = x2n[i];
				U1new[j] = U1n[i];
				U2new[j] = U2n[i];
				U3new[j] = U3n[i];
				x1_oldnew[j] = x1_oldn[i];
				x2_oldnew[j] = x2_oldn[i];
				flagnew[j] = flagn[i];
				j++;
			}
		}
		N_new = j-N_old;
/*
		st[0] = gas.m*pow( standard_deviation(&U1new[N_old], gas.n_inflow, &mn[0]),2.0 )/gas.kb;
		st[1] = gas.m*pow( standard_deviation(&U2new[N_old], gas.n_inflow, &mn[1]),2.0 )/gas.kb;
		st[2] = gas.m*pow( standard_deviation(&U3new[N_old], gas.n_inflow, &mn[2]),2.0 )/gas.kb;
*/
		free(x1);free(x2);free(U1);free(U2);free(U3);
		free(flag); free(index);
		free(x1_old); free(x2_old);

		free(x1n);free(x2n);free(U1n);free(U2n);free(U3n);
		free(flagn);free(dtn);
		free(x1_oldn); free(x2_oldn);

		gas.N = N_in;
		x1 = x1new;
		x2 = x2new;
		U1 = U1new;
		U2 = U2new;
		U3 = U3new;
		x1_old = x1_oldnew;
		x2_old = x2_oldnew;
		flag = flagnew;
		index = indexnew;
/*
		st[0] = gas.m*pow( standard_deviation(&U1[N_old], gas.n_inflow, &mn[0]),2.0 )/gas.kb;
		st[1] = gas.m*pow( standard_deviation(&U2[N_old], gas.n_inflow, &mn[1]),2.0 )/gas.kb;
		st[2] = gas.m*pow( standard_deviation(&U3[N_old], gas.n_inflow, &mn[2]),2.0 )/gas.kb;
*/
		end = omp_get_wtime();
		if (step % box.every == 0)
	      		printf("BC of dcones/flatnose done in %e sec\n", end - start);
	  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	   t = t + gas.delta_t;
     step++;
}
// end of time evolution
}

// change the factor
// gas.factor = gas.factor/( 1.0 - (p_aim - p_mean_avg)/p_aim );

  printf("nkT(1+Xnb) = %e\n",p_aim);
  printf("nkT(1+nb) = %e\n",gas.n*gas.kb*gas.T*(1.0 + (gas).b*gas.n));
  printf("nkT = %e\n",gas.n*gas.kb*gas.T);


// writting_data_to_file(box.num_steps, gas.delta_t, p, p_kin, T, &gas, E_kin, 1, p_mod);



printf("execution time = %e sec\n", end - total_start);
free(x1);
free(x2);
free(U1);
free(U2);
free(U3);
free(flag);
free(index);
if(box.problem != dcones){
	free(xi_1f);free(xi_2f); free(xi_3f);
}

for(int num=0; num < box.N[0]*box.N[1]*box.N[2]; num++)
	free(cells[num].indices_inside);
free(cells);
return 0;
}
