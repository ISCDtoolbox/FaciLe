#include "warping.h"

extern Info   info;

/* compute triangle area and unit normal in 3d */
static double area_3d(double *a,double *b,double *c,double *n) {
  double    ux,uy,uz,vx,vy,vz,dd,dd1;

  ux = b[0] - a[0];
  uy = b[1] - a[1];
  uz = b[2] - a[2];

  vx = c[0] - a[0];
  vy = c[1] - a[1];
  vz = c[2] - a[2];

  n[0] = uy*vz - uz*vy;
  n[1] = uz*vx - ux*vz;
  n[2] = ux*vy - uy*vx;
  dd   = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if ( dd > EPSD ) {
    dd1 = 1.0 / dd;
    n[0] *= dd1;
    n[1] *= dd1;
    n[2] *= dd1;
  }

  return(0.5*dd);
}


/* compute volume of tetra */
static double volume(double *a,double *b,double *c,double *d) {
  double  ax,ay,az,bx,by,bz,vol;

  ax = b[0] - a[0];
  ay = b[1] - a[1];
  az = b[2] - a[2];

  bx = c[0] - a[0];
  by = c[1] - a[1];
  bz = c[2] - a[2];

  vol = (d[0]-a[0]) * (ay*bz - az*by) + (d[1]-a[1]) * (az*bx - ax*bz) \
      + (d[2]-a[2]) * (ax*by - ay*bx);
  return(fabs(vol) / 6.0);
}


/* invert 3x3 non-symmetric matrix */
static int invmatg(double m[9],double mi[9]) {
  double  aa,bb,cc,det,vmin,vmax,maxx;
  int     k;

  /* check ill-conditionned matrix */
  vmin = vmax = fabs(m[0]);
  for (k=1; k<9; k++) {
    maxx = fabs(m[k]);
    if ( maxx < vmin )  vmin = maxx;
    else if ( maxx > vmax )  vmax = maxx;
  }
  if ( vmax == 0.0 )  return(0);

  /* compute sub-dets */
  aa = m[4]*m[8] - m[5]*m[7];
  bb = m[5]*m[6] - m[3]*m[8];
  cc = m[3]*m[7] - m[4]*m[6];
  det = m[0]*aa + m[1]*bb + m[2]*cc;
  if ( fabs(det) < EPSD )  return(0);
  det = 1.0 / det;

  mi[0] = aa*det;
  mi[3] = bb*det;
  mi[6] = cc*det;
  mi[1] = (m[2]*m[7] - m[1]*m[8])*det;
  mi[4] = (m[0]*m[8] - m[2]*m[6])*det;
  mi[7] = (m[1]*m[6] - m[0]*m[7])*det;
  mi[2] = (m[1]*m[5] - m[2]*m[4])*det;
  mi[5] = (m[2]*m[3] - m[0]*m[5])*det;
  mi[8] = (m[0]*m[4] - m[1]*m[3])*det;

  return(1);
}

static int setTGV_3d(pMesh mesh,pSol sol,pCsr A) {
  pPoint   ppt;
  int      k;

  /* at vertices */
for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if(ppt->ref ==25){
          csrSet(A,3*(k-1)+0,3*(k-1)+0,LS_TGV);
          csrSet(A,3*(k-1)+1,3*(k-1)+1,LS_TGV);
          csrSet(A,3*(k-1)+2,3*(k-1)+2,LS_TGV);
        }
      }
  
  return(1);
}


static pCsr matA_P1_3d(pMesh mesh,pSol sol) {
  pCsr     A;
  pTetra   pt;
  double  *a,*b,*c,*d,DeD[81],m[9],im[9],Ae[12][12],mm[9][12],nn[9][12],Dp[3][4];
  double   vol;
  int      i,j,k,s,ia,ja,il,ic,ig,jg,nr,nc,nbe;

	/* memory allocation (rough estimate) */
	nr  = nc = 3*mesh->np;
  nbe = 20*mesh->np;
  A   = csrNew(nr,nc,nbe,CS_UT+CS_SYM);

  memset(DeD,0,81*sizeof(double));

  /* Dp */
  Dp[0][0]=1;  Dp[0][1]=0;  Dp[0][2]=0;  Dp[0][3]=-1; 
  Dp[1][0]=0;  Dp[1][1]=1;  Dp[1][2]=0;  Dp[1][3]=-1; 
  Dp[2][0]=0;  Dp[2][1]=0;  Dp[2][2]=1;  Dp[2][3]=-1; 

  /* Fill stiffness matrix A */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->v[0] )  continue;

    /* tD E D */
    DeD[0]  = DeD[40] = DeD[80] = 2*info.mu + info.lambda;
    DeD[4]  = DeD[8]  = DeD[36] = DeD[44] = DeD[72] = DeD[76] = info.lambda;
    DeD[10] = DeD[12] = DeD[20] = DeD[24] = DeD[28] = DeD[30] = info.mu;
    DeD[50] = DeD[52] = DeD[56] = DeD[60] = DeD[68] = DeD[70] = info.mu;

    /* measure of K */
    a = &mesh->point[pt->v[0]].c[0]; 
    b = &mesh->point[pt->v[1]].c[0]; 
    c = &mesh->point[pt->v[2]].c[0]; 
    d = &mesh->point[pt->v[3]].c[0]; 

    /* mm = tB^-1 */
    for (i=0; i<3; i++) {
      m[i+0] = a[i] - d[i];
      m[i+3] = b[i] - d[i];
      m[i+6] = c[i] - d[i];
    }
    if ( !invmatg(m,im) )  return(0);
    vol = volume(a,b,c,d);

    /* mm = (tBt^-1) Dp */
    memset(mm,0,9*12*sizeof(double));
    for (i=0; i<3; i++) {
      for (j=0; j<4; j++) {
        for (s=0; s<3; s++)
          mm[i][j]   += im[i*3+s] * Dp[s][j];
        mm[i+3][j+4] = mm[i][j];
        mm[i+6][j+8] = mm[i][j];
      }
    }

    /* nn = DeD mm */
    for (i=0; i<9; i++) {
      for (j=0; j<12; j++) {
        nn[i][j] = 0.0;
        for (s=0; s<9; s++)
          nn[i][j] += DeD[i*9+s] * mm[s][j];
      }
    }

    /* Ae = vol tmm nn */
    memset(Ae,0,12*12*sizeof(double));
    for (i=0; i<12; i++) {
      for (j=i; j<12; j++) {
        for (s=0; s<9; s++)
          Ae[i][j] += vol * mm[s][i] * nn[s][j];
      }
    }

    /* stifness matrix */
    for (i=0; i<12; i++) {
      ig = pt->v[i % 4];
      ia = 3*(ig-1) + (i / 4);
      for (j=i; j<12; j++) {
        if ( fabs(Ae[i][j]) < EPSD )  continue;
        jg = pt->v[j % 4];
        ja = 3*(jg-1) + (j / 4);
        if ( ia < ja ) {
          il = ia;
          ic = ja;
        }
        else {
          il = ja;
          ic = ia;
        }
				csrPut(A,il,ic,Ae[i][j]);
      }
    }
  }
	setTGV_3d(mesh,sol,A);
	csrPack(A);
  
  return(A);
}

/* build right hand side vector and set boundary conds. */
static double *rhsF_P1_3d(pMesh mesh,pSol sol) {
  pTria    ptt;
  pPoint   ppt;
  double  *F,*vp,aire,n[3],w[3],*a,*b,*c;
  int      k,ig,size;
  char     i;

  size = sol->dim*sol->np;
  F = (double*)calloc(size,sizeof(double));
  assert(F);
  
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !ppt->ref )  continue;
      if( ppt->ref ==25 ) {
        vp = &sol->u[3*(k-1)];
          
        F[3*(k-1)+0] = LS_TGV * vp[0];
        F[3*(k-1)+1] = LS_TGV * vp[1];
        F[3*(k-1)+2] = LS_TGV * vp[2];
      }
    }

    for (k=1; k<=mesh->nt; k++) {
      ptt = &mesh->tria[k];
      if ( !ptt->v[0] )  continue;
      if( ptt->ref == mesh->ref){
     	a = &mesh->point[ptt->v[0]].c[0];
        b = &mesh->point[ptt->v[1]].c[0];
        c = &mesh->point[ptt->v[2]].c[0];
        aire = area_3d(a,b,c,n) / 3.0;
      
        for (i=0; i<3; i++) {
        ig = ptt->v[i];
        ppt = &mesh->point[ig];
      
        w[0] = info.load * n[0];
        w[1] = info.load * n[1];
        w[2] = info.load * n[2];
        for (i=0; i<3; i++) {
          ig = ptt->v[i];
          ppt = &mesh->point[ig];
          F[3*(ig-1)+0] += aire * w[0];
          F[3*(ig-1)+1] += aire * w[1]; 
          F[3*(ig-1)+2] += aire * w[2]; 
        }
       }
      }
    }
  
	return(F);
}

/* linear elasticity */
int elasti1_3d(pMesh mesh,pSol sol) {
  pCsr     A;
	double  *F,err;
	int      ier,nit;
	char     stim[32];
  
  if ( info.imprim )  fprintf(stdout,"   1.1 Assembly \n");

  /* -- Part I: matrix assembly */
  chrono(ON,&info.ctim[3]);

  /* alloc memory */
  sol->dim = mesh->dim;
  sol->ver = mesh->ver;
  sol->np  = mesh->np;
 if ( !sol->u ) {
 sol->u  = (double*)calloc(3*(mesh->npi+mesh->np2),sizeof(double));
 assert(sol->u);
 }
/* build matrix */
  A = 0;
	F = 0;
    A = matA_P1_3d(mesh,sol);
    F = rhsF_P1_3d(mesh,sol);
  chrono(OFF,&info.ctim[3]);
  printim(info.ctim[3].gdif,stim);
  if ( abs(info.imprim) > 4 )
    fprintf(stdout,"     [Time: %s]\n",stim);

  /* -- Part II: solver */
  chrono(ON,&info.ctim[4]);
	if ( info.ncpu > 1 )  csrInit(info.ncpu);
	err = sol->err;
	nit = sol->nit;
  if ( abs(info.imprim) > 4 )  fprintf(stdout,"  1.2 SOLVING  \n");
  ier = csrPrecondGrad(A,sol->u,F,&err,&nit,0);
  if  ( info.ncpu > 1 )  csrStop();
  if ( ier <= 0 )  fprintf(stdout,"  ## SOL NOT CONVERGED: ier= %d  err= %E  nit= %d\n",ier,err,nit);
  else if ( abs(info.imprim) > 4 )
    fprintf(stdout,"  %%%% CONVERGENCE: err= %E  nit= %d\n",err,nit);
	csrFree(A);
	free(F);

  chrono(OFF,&info.ctim[4]);
  printim(info.ctim[4].gdif,stim);
  if ( abs(info.imprim) > 4 )
    fprintf(stdout,"     [Time: %s]\n",stim);

  return(ier > 0);
}
