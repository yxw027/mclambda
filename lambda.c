/*------------------------------------------------------------------------------
* lambda.c : integer ambiguity resolution
*
*          Copyright (C) 2007-2008 by T.TAKASU, All rights reserved.
*
* reference :
*     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
*         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
*         1995
*     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
*         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/13 1.0 new
*           2015/05/31 1.1 add api lambda_reduction(), lambda_search()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include "mclambda\mclambda.h"
#include "mclambda\mclambda_emxAPI.h"

/* constants/macros ----------------------------------------------------------*/

#define LOOPMAX     10000           /* maximum count of search loop */

#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)

/* LD factorization (Q=L'*diag(D)*L) -----------------------------------------*/
static int LD(int n, const double *Q, double *L, double *D)
{
    int i,j,k,info=0;
    double a,*A=mat(n,n);
    
    memcpy(A,Q,sizeof(double)*n*n);
    for (i=n-1;i>=0;i--) {
        if ((D[i]=A[i+i*n])<=0.0) {info=-1; break;}
        a=sqrt(D[i]);
        for (j=0;j<=i;j++) L[i+j*n]=A[i+j*n]/a;
        for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[j+k*n]-=L[i+k*n]*L[i+j*n];
        for (j=0;j<=i;j++) L[i+j*n]/=L[i+i*n];
    }
    free(A);
    if (info) fprintf(stderr,"%s : LD factorization error\n",__FILE__);
    return info;
}
/* integer gauss transformation ----------------------------------------------*/
static void gauss(int n, double *L, double *Z, int i, int j)
{
    int k,mu;
    
    if ((mu=(int)ROUND(L[i+j*n]))!=0) {
        for (k=i;k<n;k++) L[k+n*j]-=(double)mu*L[k+i*n];
        for (k=0;k<n;k++) Z[k+n*j]-=(double)mu*Z[k+i*n];
    }
}
/* permutations --------------------------------------------------------------*/
static void perm(int n, double *L, double *D, int j, double del, double *Z)
{
    int k;
    double eta,lam,a0,a1;
    
    eta=D[j]/del;
    lam=D[j+1]*L[j+1+j*n]/del;
    D[j]=eta*D[j+1]; D[j+1]=del;
    for (k=0;k<=j-1;k++) {
        a0=L[j+k*n]; a1=L[j+1+k*n];
        L[j+k*n]=-L[j+1+j*n]*a0+a1;
        L[j+1+k*n]=eta*a0+lam*a1;
    }
    L[j+1+j*n]=lam;
    for (k=j+2;k<n;k++) SWAP(L[k+j*n],L[k+(j+1)*n]);
    for (k=0;k<n;k++) SWAP(Z[k+j*n],Z[k+(j+1)*n]);
}
/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
static void reduction(int n, double *L, double *D, double *Z)
{
    int i,j,k;
    double del;
    
    j=n-2; k=n-2;
    while (j>=0) {
        if (j<=k) for (i=j+1;i<n;i++) gauss(n,L,Z,i,j);
        del=D[j]+L[j+1+j*n]*L[j+1+j*n]*D[j+1];
        if (del+1E-6<D[j+1]) { /* compared considering numerical error */
            perm(n,L,D,j,del,Z);
            k=j; j=n-2;
        }
        else j--;
    }
}
/* modified lambda (mlambda) search (ref. [2]) -------------------------------*/
static int search(int n, int m, const double *L, const double *D,
                  const double *zs, double *zn, double *s)
{
    int i,j,k,c,nn=0,imax=0;
    double newdist,maxdist=1E99,y;
    double *S=zeros(n,n),*dist=mat(n,1),*zb=mat(n,1),*z=mat(n,1),*step=mat(n,1);
    
    k=n-1; dist[k]=0.0;
    zb[k]=zs[k];
    z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
    for (c=0;c<LOOPMAX;c++) {
        newdist=dist[k]+y*y/D[k];
        if (newdist<maxdist) {
            if (k!=0) {
                dist[--k]=newdist;
                for (i=0;i<=k;i++)
                    S[k+i*n]=S[k+1+i*n]+(z[k+1]-zb[k+1])*L[k+1+i*n];
                zb[k]=zs[k]+S[k+k*n];
                z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
            }
            else {
                if (nn<m) {
                    if (nn==0||newdist>s[imax]) imax=nn;
                    for (i=0;i<n;i++) zn[i+nn*n]=z[i];
                    s[nn++]=newdist;
                }
                else {
                    if (newdist<s[imax]) {
                        for (i=0;i<n;i++) zn[i+imax*n]=z[i];
                        s[imax]=newdist;
                        for (i=imax=0;i<m;i++) if (s[imax]<s[i]) imax=i;
                    }
                    maxdist=s[imax];
                }
                z[0]+=step[0]; y=zb[0]-z[0]; step[0]=-step[0]-SGN(step[0]);
            }
        }
        else {
            if (k==n-1) break;
            else {
                k++;
                z[k]+=step[k]; y=zb[k]-z[k]; step[k]=-step[k]-SGN(step[k]);
            }
        }
    }
    for (i=0;i<m-1;i++) { /* sort by s */
        for (j=i+1;j<m;j++) {
            if (s[i]<s[j]) continue;
            SWAP(s[i],s[j]);
            for (k=0;k<n;k++) SWAP(zn[k+i*n],zn[k+j*n]);
        }
    }
    free(S); free(dist); free(zb); free(z); free(step);
    
    if (c>=LOOPMAX) {
        fprintf(stderr,"%s : search loop count overflow\n",__FILE__);
        return -1;
    }
    return 0;
}
/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
extern int lambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s)
{
    //FILE *f = fopen("file.txt", "w");
    int info, i;
    double *L,*D,*Z,*z,*E;
    
    if (n<=0||m<=0) return -1;
    L=zeros(n,n); D=mat(n,1); Z=eye(n); z=mat(n,1); E=mat(n,m);
    
    /* LD factorization */
    if (!(info=LD(n,Q,L,D))) {
        
        /* lambda reduction */
        reduction(n,L,D,Z);
        matmul("TN",n,1,n,1.0,Z,a,0.0,z); /* z=Z'*a */
        
        /* mlambda search */
        if (!(info=search(n,m,L,D,z,E,s))) {
            
            info=solve("T",Z,E,n,m,F); /* F=Z'\E */
        }
    }
    free(L); free(D); free(Z); free(z); free(E);

    // Print only for testing
    /*for (i=0; i<n; i++){
        fprintf(f, "Amb: %f\n", a[i]);
        fprintf(f, "SolAmb: %f\n", F[i]);
    }
    fprintf(f, "S1: %f\n", s[0]);
    fprintf(f, "S2: %f\n", s[1]);
    fclose(f);*/

    return info;
}
/* lambda reduction ------------------------------------------------------------
* reduction by lambda (ref [1]) for integer least square
* args   : int    n      I  number of float parameters
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *Z     O  lambda reduction matrix (n x n)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_reduction(int n, const double *Q, double *Z)
{
    double *L,*D;
    int i,j,info;
    
    if (n<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    for (i=0;i<n;i++) for (j=0;j<n;j++) {
        Z[i+j*n]=i==j?1.0:0.0;
    }
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* lambda reduction */
    reduction(n,L,D,Z);
     
    free(L); free(D);
    return 0;
}
/* mlambda search --------------------------------------------------------------
* search by  mlambda (ref [2]) for integer least square
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_search(int n, int m, const double *a, const double *Q,
                         double *F, double *s)
{
    double *L,*D;
    int info;
    
    if (n<=0||m<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* mlambda search */
    info=search(n,m,L,D,a,F,s);
    
    free(L); free(D);
    return info;
}

/* MC-LAMBDA RUN FUNCTION ----------------------------------------------------
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
// --------------------------------------------------------------------------
//                           MC-LAMBDA RUN FUNCTION
// --------------------------------------------------------------------------
extern int mclambda_exec(rtk_t *rtk, int n, int m, const double *a, const double *Q, double *F,
                  double *s)
{

    //FILE *f = fopen("file.txt", "w");
    // Return variable definition
    int info;
    int i, j;
    int rtk_opt;

    // Definiton of variables MCLAMBDA
    double ahat[999]; 
    double Qahat[999];
    double method;
    double param;
	emxArray_char_T *type;
    double value;

    // Output variables MCLAMBDA
    emxArray_real_T *afixed;
    emxArray_real_T *sqnorm;
    double Ps;
    double Qzhat[999];
    double Z1[999];
    double nfixed;
	double mu;

    // Initialize input argument 'type'.
    static int iv1[2] = { 2, 2 };
    int idx0;
    int idx1;
    // Set the size of the array.
    // Change this size to the value that the application requires.
    type = emxCreateND_char_T(2, *(int (*)[2])&iv1[0]);

    // Format input data
    //fprintf(f, "NAmb: %d\n", n);
    //fprintf(f, "Ambiguities\n"); 
    for (i=0; i<n; i++){ // Building ambiguity vector
        ahat[i] = a[i];
        //fprintf(f, "Amb: %f\n", a[i]);
    }

    //fprintf(f, "Variance-Covariance\n");
    for (i=0; i<n; i++){ // Building VC ambiguity matrix
        for (j=0; j<n; j++){
            Qahat[i + j*n] = Q[i + j*n];
        }
        //fprintf(f, "VC: %f\n", Qahat[i]);
    }
    /*for (i=0; i<n*n; i++){
        Qahat[i] = Q[i];
    }*/ 

    // Obtain options from data input
    rtk_opt = rtk->opt.mcopt;
    if (rtk_opt == MCLAMBDA_ILS_SAS){
        method = 1.0;
        param = 1.0;
        value = 10;
        //fprintf(f, "Method ILS SAS\n");
    } else if (rtk_opt == MCLAMBDA_ILS_ES){
        method = 2.0;
        param = 1.0;
        value = 10;
        //fprintf(f, "Method ILS ES\n");
    } else if (rtk_opt == MCLAMBDA_IRM){
        method = 3.0;
        param = 0.0;
        value = 0;
        //fprintf(f, "Method IRM\n");
    } else if (rtk_opt == MCLAMBDA_IB){
        method = 4.0;
        param = 0.0;
        value = 0;
        //fprintf(f, "Method IB\n");
    } else if (rtk_opt == MCLAMBDA_ILS_PO){
        method = 5.0;
        param = 1.0;
        value = 0.01;
        //fprintf(f, "Method ILS PO\n");
    } else if (rtk_opt == MCLAMBDA_ILS_MU){
        method = 5.0;
        param = 1.0;
        value = 0.5;
        //fprintf(f, "Method ILS MU\n");
    } else{ 
        method = 1.0;
        param = 1.0;
        value = 10;
        //fprintf(f, "Method ILS SAS default\n");
    }
    
    // Comment uncomment in order to change the method of calculation
    //type->data[0] = 'P';
    //type->data[1] = 'O';
    type->data[0] = 'M';
    //type->data[1] = 'U';
    //type->data[0] = 'N';
    //type->data[1] = 'C';
    //type->data[2] = 'A';
    //type->data[3] = 'N';
    //type->data[4] = 'D';
    //type->data[5] = 'S';

	emxInitArray_real_T(&afixed, 2);
    emxInitArray_real_T(&sqnorm, 2);

    // Testing
    if (n<=0||m<=0) return -1;

    // Execute MCLAMBDA Routine
    mclambda(n, ahat, Qahat, method, param, type, value, afixed, sqnorm, &Ps, Qzhat, Z1, &nfixed, &mu);
    
    //free(L); free(D); free(Z); free(z); free(E);

    // Preparing data to return (Print in file only for tests)
    if (method==3.0 || method==4.0){
        s[0] = 0;
        s[1] = 0;
    } else{
        s = sqnorm->data;
    }

    for (i=0; i<n; i++){
        F[i] = afixed->data[i];
        //fprintf(f, "SolAmb: %f\n", F[i]);
    }
    //fprintf(f, "S1: %f\n", s[0]);
    //fprintf(f, "S2: %f\n", s[1]);
    //fclose(f);
    info = 0;

    return info;
}
// --------------------------------------------------------------------------
//                         END MC-LAMBDA RUN FUNCTION
// --------------------------------------------------------------------------
