//#include <stdio.h>
//#include <omp.h>
//#define NUM_THREADS 16
//static long num_steps = 100000000;
//double step;
//int main ()
//{
//    double pi=0;
//    double start_time, run_time;
//    omp_set_num_threads(NUM_THREADS);
//
//    step = 1.0/(double) num_steps;
//    start_time = omp_get_wtime();
////    omp_set_schedule(omp_sched_static, 1);
//#pragma omp parallel for reduction(+:pi) schedule(auto)
//    for (int i=1;i<= num_steps; i++){
//        double x = (i-0.5)*step;
//        pi += 4.0/(1.0+x*x);
//    }
//    pi *= step;
//    run_time = omp_get_wtime() - start_time;
//    printf("\n pi with %ld steps is %lf in %lf seconds\n ",num_steps,pi,run_time);
//}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

# define NPOINTS 1000
# define MAXITER 1000



struct d_complex{
    double r;
    double i;
};

struct d_complex c;

void testpoint(d_complex);

int numoutside = 0;

int main(){
    int i, j;
    double area, error, eps  = 1.0e-5;


//   Loop over grid of points in the complex plane which contains the Mandelbrot set,
//   testing each point to see whether it is inside or outside the set.

#pragma omp parallel for default(none) private(c, i, j) firstprivate(eps)
    for (i=0; i<NPOINTS; i++) {
        for (j=0; j<NPOINTS; j++) {
            c.r = -2.0+2.5*(double)(i)/(double)(NPOINTS)+eps;
            c.i = 1.125*(double)(j)/(double)(NPOINTS)+eps;
            testpoint(c);
        }
    }

// Calculate area of set and error estimate and output the results

    area=2.0*2.5*1.125*(double)(NPOINTS*NPOINTS-numoutside)/(double)(NPOINTS*NPOINTS);
    error=area/(double)NPOINTS;

    printf("Area of Mandlebrot set = %12.8f +/- %12.8f\n",area,error);
    printf("Correct answer should be around 1.510659\n");

}

void testpoint(d_complex c){

// Does the iteration z=z*z+c, until |z| > 2 when point is known to be outside set
// If loop count reaches MAXITER, point is considered to be inside the set

    d_complex z;
    int iter;
    double temp;

    z=c;
    for (iter=0; iter<MAXITER; iter++){
        temp = (z.r*z.r)-(z.i*z.i)+c.r;
        z.i = z.r*z.i*2+c.i;
        z.r = temp;
        if ((z.r*z.r+z.i*z.i)>4.0) {
#pragma omp atomic
            numoutside++;
            break;
        }
    }

}