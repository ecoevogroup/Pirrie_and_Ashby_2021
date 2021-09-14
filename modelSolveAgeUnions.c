/*****************************************************************************************************
 * Solve the model using an adaptive step size RK method
 * Tested 11/05/21. Needs to be compiled in MATLAB as a MEX file
 * For example using the commands:
 * "mex -setup C; mex modelSolveAgeUnions.c -R2018a"
 ****************************************************************************************************/

#include <mex.h>
#include <math.h>

/***********************************
 * CONSTANT PARAMATER VALUES
 ***********************************/
#define MAXSTEPS 1e6 /* Maximum number of steps for ODP solver */
#define EPS 1e-9 /* ODE solver tolerance */
#define TINY 1e-30 /* Constant value for solver */
#define b21 0.2
#define b31 3.0/40.0
#define b32 9.0/40.0
#define b41 0.3
#define b42 -0.9
#define b43 1.2
#define b51 -11.0/54.0
#define b52 2.5
#define b53 -70.0/27.0
#define b54 35.0/27.0
#define b61 1631.0/55296
#define b62 175.0/512.0
#define b63 575.0/13824.0
#define b64 44275.0/110592
#define b65 253.0/4096.0
#define c1 37.0/378.0
#define c3 250.0/621.0
#define c4 125.0/594.0
#define c6 512.0/1771.0
#define dc5 -277.00/14336
// #define isnan(x) ((x) != (x))

// Define all the parameters 
struct PARAM{
    double t_max;
    double MJ;
    double FJ;
    double MA;
    double FA;
    double U;
    double b;
    double q;
    double s;
    double uMJ;
    double uFJ;
    double uMA;
    double uFA;
    double tau;
    double delta;
    double d;
};

/*************************************
 * Function prototypes
 *************************************/
int my_rungkut (double *T, double *xOut, struct PARAM *p);
void rkqs(double *x, double *dxdt, double *h, double *hnext, double *xScale, struct PARAM *p);
void rkck(double *x, double *dxdt, double *xOut, double *xErr, double h, struct PARAM *p);
void dynamic(double *x, double *dxdt, struct PARAM *p);
double FMAX(double, double);
double FMIN(double, double);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Define local variables */
    double *t, *x, *parameter, *tTemp, *xTemp;
    int i, j, colLen, maxsteps;
    struct PARAM p;
    
    /* Allocate inputs */
    if(nrhs!=16){
        mexErrMsgTxt("Incorrect number of input arguments!\n");
    }
    else{
        /* Converts input structure into local structure */
        parameter= mxGetPr(prhs[0]);
        p.t_max= *parameter;
        parameter= mxGetPr(prhs[1]);
        p.MJ= *parameter;
        parameter= mxGetPr(prhs[2]);
        p.FJ= *parameter;
        parameter= mxGetPr(prhs[3]);
        p.MA= *parameter;
        parameter= mxGetPr(prhs[4]);
        p.FA= *parameter;
        parameter= mxGetPr(prhs[5]);
        p.U= *parameter;
        parameter= mxGetPr(prhs[6]);
        p.b= *parameter;
        parameter= mxGetPr(prhs[7]);
        p.q= *parameter;
        parameter= mxGetPr(prhs[8]);
        p.s= *parameter;
        parameter= mxGetPr(prhs[9]);
        p.uMJ= *parameter;
        parameter= mxGetPr(prhs[10]);
        p.uFJ= *parameter;
        parameter= mxGetPr(prhs[11]);
        p.uMA= *parameter;
        parameter= mxGetPr(prhs[12]);
        p.uFA= *parameter;
        parameter= mxGetPr(prhs[13]);
        p.tau= *parameter;
        parameter= mxGetPr(prhs[14]);
        p.delta= *parameter;
        parameter= mxGetPr(prhs[15]);
        p.d= *parameter;
    }
    maxsteps = (int)MAXSTEPS; // Casts MAXSTEPS into a local integer
    
    
    
    /* Allocate memory */
    /* sizeof(double) because there are a specific number of memory blocks
    for a double */
    tTemp = malloc(maxsteps*sizeof(double)); 
    xTemp = malloc(maxsteps*5*sizeof(double));
    
    /* Call ODP solver */
    colLen = my_rungkut(tTemp, xTemp, &p);
    // colLen is the number of timesteps
    /* Create outputs */
    plhs[0] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(colLen, 5, mxREAL);
    
    t = mxGetPr(plhs[0]);
    x = mxGetPr(plhs[1]);
    
    /* Copy data to outputs */
    for (i=0;i<colLen;i++){
        t[i] = tTemp[i];
        for (j=0;j<5;j++) {
            x[i + j*colLen] = xTemp[i + j*maxsteps];
        }
    }
    
    /* Free memory */
    free(tTemp);
    free(xTemp);
    
    return;
}

/*****************************************
 * ODP solver
 ****************************************/
int my_rungkut (double *T, double *xOut, struct PARAM *p){
    
    double t, x[5], dxdt[5], xScale[5], hnext[1], h[1];
    int i, j, k, exitflag, count, maxsteps;
    
    /* Other parameters */
    exitflag = 1;
    count=0; // Number of timesteps
    k=1;
    h[0] = 1e-3;
    hnext[0] = 1e-3;
    t=0;
    maxsteps = (int)MAXSTEPS;
    
    /* Initialise populations */
    x[0] = p->MJ;
    x[1] = p->FJ;
    x[2] = p->MA;
    x[3] = p->FA;
    x[4] = p->U;
    
    /* Update output */
    T[0]=t;
    for (i=0; i<5; i++) {
        xOut[i*maxsteps] = x[i]; // Initialising the output
    }
    
    /* Main loop: */
    do{
        /* This ensures the final step lands us on the final time point 
         because of the adaptive step size */
        if(1.1*hnext[0]>(p->t_max-t)){
            // If the next step takes you past t_max shrink the step
            hnext[0] = p->t_max-t;
            h[0] = p->t_max-t;
            t=p->t_max;
            exitflag=0; // Updates exit flag to indicate it's at t_max
        }
        else{ // Add next timestep on
            h[0] = hnext[0];
            t+=h[0];
        }
        if(t>=p->t_max) {
            t=p->t_max;
            exitflag=0; // Updates exit flag to indicate it's at t_max
        }
        /* This is where the equations are first solved */
        dynamic(x, dxdt, p);

        /* Adjust the step size to maintain accuracy */
        for (i=0; i<5; i++){
            x[i] = FMAX(x[i],0); // Stops population going negative
            // Gets an idea of the next population size
            xScale[i]=fabs(x[i])+fabs(dxdt[i]*(*h))+TINY; 
        }
        
        rkqs(x, dxdt, h, hnext, xScale, p);
        
        /* Stops population going below zero */ 
        for (i=0; i<5; i++){
            x[i] = FMAX(x[i],0);
        }
        
        /* Update output */
        count++;
        T[count] = t;
        for (i=0; i<5; i++) {
            xOut[count + i*maxsteps] = x[i];
        }
    }while(count<(maxsteps-1) && t<=p->t_max && exitflag);
    count++;
    
    return count;
}

/***************************************
 * This generates the adaptive step-size (using Cash-Karp method)
 **************************************/
void rkqs(double *x, double *dxdt, double *h, double *hnext, double *xScale, struct PARAM *p)
{
    double xTemp[5], xErr[5], htemp, errmax;
    int i, j, count;
    
    count = 0;
    while(count<1e5)
    {
        rkck(x, dxdt, xTemp, xErr, *h, p);
        
        errmax= 0.0;
        for(i=0;i<5;i++){
            errmax= FMAX(errmax, fabs(xErr[i]/(xScale[i])));
        }
        errmax/= EPS;
        if(errmax<=1.0) break;
        htemp= 0.9*(*h)*pow(errmax, -0.25);
        *h= (*h>=0.0 ? FMAX(htemp, 0.1*(*h)) : FMIN(htemp, 0.1*(*h)));
        count++;
        if(count>1e4){
            mexErrMsgTxt("stuck in loop!\n");
            break;
        }
    }    
    if(errmax > 1.89E-4) {
        *hnext= 0.9*(*h)*pow(errmax, -0.2);
    }
    else {
        *hnext= 5.0*(*h);
    }
    
    for(i=0;i<5;i++){
        x[i] = xTemp[i];
    }
}

/**************************************
 * Standard RK solver
 **************************************/
void rkck(double *x, double *dxdt, double *xOut, double *xErr, double h, struct PARAM *p){
    
    int i, j;
    double xk1[5], xk2[5], xk3[5], xk4[5], xk5[5], xk6[5], xTemp[5];
    double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0, dc4=c4-13525.0/55296.0,
            dc6=c6-0.25;
    
    for(i=0;i<5;i++){
        xTemp[i] = x[i] + b21*h*dxdt[i];
    }
    dynamic(xTemp, xk2, p);
    
    for(i=0;i<5;i++){
        xTemp[i] = x[i]+h*(b31*dxdt[i]+b32*xk2[i]);
    }
    dynamic(xTemp, xk3, p);
    
    for(i=0;i<5;i++){
        xTemp[i] = x[i]+h*(b41*dxdt[i]+b42*xk2[i]+b43*xk3[i]);
    }
    dynamic(xTemp, xk4, p);
    
    for(i=0;i<5;i++){
        xTemp[i] = x[i]+h*(b51*dxdt[i]+b52*xk2[i]+b53*xk3[i]+b54*xk4[i]);
    }
    dynamic(xTemp, xk5, p);
    
    for(i=0;i<5;i++){
        xTemp[i] = x[i]+h*(b61*dxdt[i]+b62*xk2[i]+b63*xk3[i]+b64*xk4[i]+b65*xk5[i]);
    }
    dynamic(xTemp, xk6, p);
    
    for(i=0;i<5;i++){
        xOut[i]= x[i]+h*(c1*dxdt[i]+c3*xk3[i]+c4*xk4[i]+c6*xk6[i]);
        xErr[i]= h*(dc1*dxdt[i]+dc3*xk3[i]+dc4*xk4[i]+dc5*xk5[i]+dc6*xk6[i]);
    }
}

/**************************************
 * Population and evolutionary dynamics
 **************************************/
void dynamic(double *x, double *dxdt, struct PARAM *p){
    
//     Define parameters
    double MJ, FJ, MA, FA, U, births, NJ, Nu, NA, N;
   
    
//     Define populations
    MJ = x[0];
    FJ = x[1];
    MA = x[2];
    FA = x[3];
    U = x[4];
    
    NJ = MJ + FJ;
    Nu = MA + FA;
    NA = Nu + 2*U;
    N = NJ + NA;
    
//     Note: if you have a birth term, ensure it is greater than zero
    births = FMAX(p->b*(1 - p->q*N)*U,TINY);
    
//     Calculate derivatives
    dxdt[0] = births*p->s - (p->uMJ + p->tau)*MJ;
    dxdt[1] = births*(1-p->s) - (p->uFJ + p->tau)*FJ;
    dxdt[2] = p->tau*MJ - p->uMA*MA - (p->delta*MA*FA)/FMAX(Nu,TINY) + (p->uFA + p->d)*U;
    dxdt[3] = p->tau*FJ - p->uFA*FA - (p->delta*MA*FA)/FMAX(Nu,TINY) + (p->uMA + p->d)*U;
    dxdt[4] = (p->delta*MA*FA)/FMAX(Nu,TINY) - (p->uMA + p->uFA + p->d)*U;
}

/***************************************
 * Return maximum of two inputs
 ***************************************/
double FMAX(double l, double r)
{
    if(l>r)return l;
    else   return r;
}

/***************************************
 * Return minimum of two inputs
 ***************************************/
double FMIN(double l, double r)
{
    if(l<r)return l;
    else   return r;
}
