#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <memory>
#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include "H5Cpp.h"
using namespace H5;

#define NMAX 16384
#define SQR(x) ((x)*(x))

void set_caculation_parameters(int *n, double *eps2, double *dt, double *t_end, double *t_out, double *r_v);
void set_from_param_txt(int *n, double *eps2, double *dt, double *t_end, double *t_out, double *r_v);
double gaussian(void);
void make_spherical_df(int n, double m[], double x[][3], double v[][3], double r_v, double eps2);
void check_partcile_distribution(int n, double x[][3]);
void calc_force(int n, double m[], double x[][3], double a[][3], double eps2);
void leap_frog(int n, double m[], double x[][3], double v[][3], double a[][3], double dt, double eps2);
double calc_energy(int n, double m[], double x[][3], double v[][3], double eps2);
void animated_snapshot(int n, double t, double x[][3], FILE **gp);
void open_window(FILE **gp);
void close_window(FILE **gp);
void virial_ratio(int n, double m[], double x[][3], double v[][3], double eps2, double t);

struct Snapshot {
    int index;
    double t;
    std::vector<std::array<double, 3>> x;
    std::vector<std::array<double, 3>> v;
    std::vector<double> m;
};

void WriteSnapshotsToHDF5(const std::string& filename, const std::vector<Snapshot>& snapshots) {
    try {
        H5File file(filename, H5F_ACC_TRUNC);  // Create new HDF5 file

        for (const auto& snapshot : snapshots) {
            std::ostringstream oss;
            oss << "/snapshot_" << std::setw(3) << std::setfill('0') << snapshot.index;
            std::string group_name = oss.str();
            Group group(file.createGroup(group_name));

            // Save time 't'
            hsize_t dim_scalar = 1;
            DataSpace scalar_space(1, &dim_scalar);
            DataSet dataset_t = file.createDataSet(group_name + "/t", PredType::NATIVE_DOUBLE, scalar_space);
            dataset_t.write(&snapshot.t, PredType::NATIVE_DOUBLE);

            // Save positions 'x'
            hsize_t dims_x[2] = { snapshot.x.size(), 3 };
            DataSpace dataspace_x(2, dims_x);
            DataSet dataset_x = file.createDataSet(group_name + "/x", PredType::NATIVE_DOUBLE, dataspace_x);
            dataset_x.write(snapshot.x.data(), PredType::NATIVE_DOUBLE);

            // Save velocities 'v'
            hsize_t dims_v[2] = { snapshot.v.size(), 3 };
            DataSpace dataspace_v(2, dims_v);
            DataSet dataset_v = file.createDataSet(group_name + "/v", PredType::NATIVE_DOUBLE, dataspace_v);
            dataset_v.write(snapshot.v.data(), PredType::NATIVE_DOUBLE);

            // Save masses 'm'
            hsize_t dims_m = snapshot.m.size();
            DataSpace dataspace_m(1, &dims_m);
            DataSet dataset_m = file.createDataSet(group_name + "/m", PredType::NATIVE_DOUBLE, dataspace_m);
            dataset_m.write(snapshot.m.data(), PredType::NATIVE_DOUBLE);
        }

        std::cout << "HDF5 file created: " << filename << std::endl;
    } catch (FileIException& error) {
        error.printErrorStack();
    } catch (DataSetIException& error) {
        error.printErrorStack();
    } catch (DataSpaceIException& error) {
        error.printErrorStack();
    }
}





// main function -----------------------


int main(void){
   /* number of partciles */
   int n;

   /* particle data */
   static double m[NMAX], x[NMAX][3], v[NMAX][3], a[NMAX][3];

   /* system energy, Virial ratio */
   double e, e_init, r_v;

   /* current time, timestep, end time, data interval */
   double t, dt, t_end, t_out;

   /* squared softning parameter */
   double eps2;

   /* Set the parameters for the calculation */ 
   set_from_param_txt(&n, &eps2, &dt, &t_end, &t_out, &r_v);
   // set_caculation_parameters(&n, &eps2, &dt, &t_end, &t_out, &r_v);
   eps2 = SQR(eps2);

   make_spherical_df(n, m, x, v, r_v, eps2);

   // check_partcile_distribution(n, x);

   e_init = calc_energy(n, m, x, v, eps2);

   calc_force(n, m, x, a, eps2);
   
   std::vector<Snapshot> snapshots;
   Snapshot snapshot;
   int i = 0;

   while(t < t_end){

      leap_frog(n, m, x, v, a, dt, eps2);

      if(fmod(t, t_out) == 0.0){
         
         snapshot.index = i;
         snapshot.t = t;
         snapshot.x.resize(n);
         snapshot.v.resize(n);
         snapshot.m.resize(n);

         std::copy(&x[0][0], &x[n][0], &snapshot.x[0][0]);
         std::copy(&v[0][0], &v[n][0], &snapshot.v[0][0]);
         std::copy(&m[0], &m[n], snapshot.m.begin());

        snapshots.push_back(snapshot);

        printf("%d %f\n", i, t);
        i++;
      }

      t += dt;
   }


   WriteSnapshotsToHDF5("snapshots.hdf5", snapshots);
   return 0;

}



// define functions -----------------------

void set_from_param_txt(int *n, double *eps2, double *dt, double *t_end, double *t_out, double *r_v){
   FILE *file = fopen("param.txt", "r");
   if (file == NULL) {
      fprintf(stderr, "Error opening param.txt\n");
      exit(EXIT_FAILURE);
   }

   fscanf(file, "n %d\n", n);
   fscanf(file, "eps2 %lf\n", eps2);
   fscanf(file, "dt %lf\n", dt);
   fscanf(file, "t_end %lf\n", t_end);
   fscanf(file, "t_out %lf\n", t_out);
   fscanf(file, "r_v %lf\n", r_v);

   fclose(file);
   
}


void set_caculation_parameters(int *n, double *eps2, double *dt, double *t_end, double *t_out, double *r_v){

   // Set the parameters for the calculation
   fprintf(stderr, "number of particles(n): ");
   scanf("%d", n);

   fprintf(stderr, "softning(eps2): ");
   scanf("%lf", eps2);
   
   fprintf(stderr, "timestep(dt): ");
   scanf("%lf", dt);
   
   fprintf(stderr, "end time(t_end): ");
   scanf("%lf", t_end);
   
   fprintf(stderr, "data interval(t_out): ");
   scanf("%lf", t_out);

   fprintf(stderr, "virial ratio(r_v): ");
   scanf("%lf", r_v);

   // Print the parameters
   fprintf(stderr, "number of particles(n) = %d\n", *n);
   fprintf(stderr, "softning(eps2) = %g\n", *eps2);
   fprintf(stderr, "timestep(dt) = %g\n", *dt);
   fprintf(stderr, "end time(t_end) = %g\n", *t_end);
   fprintf(stderr, "data interval(t_out) = %g\n", *t_out);
   fprintf(stderr, "virial ratio(r_v) = %g\n", *r_v);
}

double gaussian(void)
/* Gaussian with mean = 0.0 and dispersion = 0.0 by Box-Muller method */
{
  double x, y, r2;
  double z;

  do{
    x = 2.0*drand48() - 1.0;
    y = 2.0*drand48() - 1.0;
    r2 = x*x + y*y;
  }while(r2 >= 1.0 || r2 == 0.0);
  z = sqrt(-2.0*log(r2)/r2)*x; /* discard another Gaussian */

  return(z);
}


void make_spherical_df(int n, double m[], double x[][3], double v[][3], double r_v, double eps2){

   double W = 0.0;
   double r2, sigma;

   for (int i = 0; i < n; i++) {
      double _x = drand48()*2 - 1;
      double _y = drand48()*2 - 1;
      double _z = drand48()*2 - 1;
      m[i] = 1.0 / n; // Assign equal mass to each particle
      double r = sqrt(_x * _x + _y * _y + _z * _z); // Radial distance

      if (r < 1.0) {
         x[i][0] = _x;
         x[i][1] = _y;
         x[i][2] = _z;
      }else{
         i--;
      }
   }

   for (int i = 0; i < n-1; i++) {
      for (int j=i+1; j<n; j++) {
         r2 = SQR(x[j][0]-x[i][0]) + SQR(x[j][1]-x[i][1]) + SQR(x[j][2]-x[i][2]);
         W = W - m[i]*m[j]/sqrt(r2 + eps2);
      }
   }

   sigma = sqrt(2.0*r_v*fabs(W)/3.0);

   for (int i = 0; i < n; i++) {
      for (int j = 0; j < 3; j++) {
         v[i][j] = sigma*gaussian();
      }
   }
}

void check_partcile_distribution(int n, double x[][3]){
   for (int i = 0; i < n; i++) {
      printf("%f %f %f\n", x[i][0], x[i][1], x[i][2]);
   }
}

void calc_force(int n, double m[], double x[][3], double a[][3], double eps2){

   for(int i=0;i<n;i++){
      for(int k=0;k<3;k++){
         a[i][k] = 0.0;
      }
   }

   double r[3], r3inv;
    for (int i = 0; i < n; i++){
      for (int j = i + 1; j < n; j++){
         for (int k = 0; k < 3; k++){
            r[k] = x[j][k] - x[i][k];
         }
            
         r3inv = 1.0/pow(sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + eps2), 3.0);
         for (int k = 0; k < 3; k++){
            a[i][k] += m[j]*r[k]*r3inv;
            a[j][k] += -m[i]*r[k]*r3inv;
         }
        }
    }
}

void leap_frog(int n, double m[], double x[][3], double v[][3], double a[][3], double dt, double eps2){

   double v_half[n][3];

   for(int i=0;i<n;i++){
      for(int k=0;k<3;k++){
         v_half[i][k] = v[i][k] + 0.5*a[i][k]*dt;
         x[i][k] += v_half[i][k]*dt;
      }
   }

   calc_force(n, m, x, a, eps2);

   for(int i=0;i<n;i++){
      for(int k=0;k<3;k++){
         v[i][k] = v_half[i][k] + 0.5*a[i][k]*dt;
      }
   }
}

double calc_energy(int n, double m[], double x[][3], double v[][3], double eps2){

   /* kinetic energy */
   double K = 0.0;

   for(int i=0; i<n; i++){
      K += m[i]*(SQR(v[i][0]) + SQR(v[i][1]) + SQR(v[i][2]));
   }
   K *= 0.5;

   /* potential energy */
   double U = 0.0;
   double r2;

   for(int i=0; i<n-1; i++){
      for(int j=i+1; j<n; j++){
         r2 = SQR(x[j][0]-x[i][0]) + SQR(x[j][1]-x[i][1]) + SQR(x[j][2]-x[i][2]);
         U -= m[i]*m[j]/sqrt(r2 + eps2);
      }
   }

   return (K + U);
}

void open_window(FILE **gp)
{

  if ((*gp = popen("gnuplot -persist", "w")) == NULL) {
    printf("gnuplot open error!!\n");
    exit(EXIT_FAILURE);
  }
  fprintf(*gp, "set size square\n");
//   fprintf(*gp, "set terminal gif animate\n");
//   fprintf(*gp, "set output \"plot.gif\"\n");
  fprintf(*gp, "set xrange [-%f:%f]\n", VMAX, VMAX);
  fprintf(*gp, "set yrange [-%f:%f]\n", VMAX, VMAX);
  
}
  
  
void close_window(FILE **gp)
{

  pclose(*gp);

}


void animated_snapshot(int n, double t, double x[][3], FILE **gp)
{

  int i;
  fprintf(*gp, "set key title \"%f\"\n", t);
  fprintf(*gp, "plot '-' with points pointtype 0 notitle \n");
  for(i=0; i<=n; i++){
    fprintf(*gp,"%f\t%f\n", x[i][0], x[i][1]);
  }
  fprintf(*gp,"e\n");

}

void virial_ratio(int n, double m[], double x[][3], double v[][3], double eps2, double t)
{
   /* kinetic energy */
   double K = 0.0;

   for(int i=0; i<n; i++){
      K += m[i]*(SQR(v[i][0]) + SQR(v[i][1]) + SQR(v[i][2]));
   }
   K *= 0.5;

   /* potential energy */
   double U = 0.0;
   double r2;

   for(int i=0; i<n-1; i++){
      for(int j=i+1; j<n; j++){
         r2 = SQR(x[j][0]-x[i][0]) + SQR(x[j][1]-x[i][1]) + SQR(x[j][2]-x[i][2]);
         U -= m[i]*m[j]/sqrt(r2 + eps2);
      }
   }
   double rr = -K/U;
   printf("%f %f\n", t, rr);
}