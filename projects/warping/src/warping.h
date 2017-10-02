#ifndef _WARPING_H
#define _WARPING_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>

#include "chrono.h"
#include "libmesh5.h"
#include "memory.h"
#include "sparse.h"

#define LS_VER   "1.0a"
#define LS_REL   "Sept, 2017"
#define LS_CPY   "Copyright (c) ISCD "
#define LS_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

#define LS_Ver       (1<<0)
#define LS_Edg       (1<<1)
#define LS_Tri       (1<<2)
#define LS_Tet       (1<<3)

#define LS_LAMBDA     10.0e5
#define LS_MU          8.2e5
#define LS_E          10000
#define LS_NU         0.1
#define LS_MAT        50
#define LS_CL         50
#define LS_RES        10.e-6
#define LS_MAXIT      10000
#define LS_TGV        1.e+30

#define MAXIT         100
#define LOAD          30 

#define LS_MAX(a,b)   ( ((a) < (b)) ? (b) : (a) )
#define LS_MIN(a,b)   ( ((a) < (b)) ? (a) : (b) )


#define PRECI    1.0
#define EPS      1.e-6
#define EPSD     1.e-30
#define EPSA     1.e-200

typedef struct {
  double   c[3];
  int      ref,s;
} Point;
typedef Point * pPoint;

typedef struct {
  int     v[6],ref;
  double g[3];
} Tria;
typedef Tria * pTria;

/* au cas ou ... */
typedef struct {
    int     v[8],ref;
} Quad;
typedef Quad * pQuad;

typedef struct {
  int     v[10],ref;
} Tetra;
typedef Tetra * pTetra;

typedef struct {
  double   delta,min[3],max[3];
  double   lambda, mu, e, nu, load;
  int      ncpu,nit,debug;
  char     imprim,cg,rhs;
  mytime   ctim[TIMEMAX];
} Info;

typedef struct {
  int      np,np2,na,nt,ne,npi,nai,nq,nti,nei,hmax,hcur,ver,dim,mark,ref;
  char     *name,*nameout;
  double   min[3],max[3],Ray[3],o[3];
  
  pPoint   point;
  pTria    tria;
  pQuad    quad;
  pTetra   tetra;
} Mesh;
typedef Mesh * pMesh;

typedef struct {
  int      dim,ver,np,na,ne,nit,iter;
  double  *u,*bc,*u0,err;
  char    *namein,*nameout,cltyp;
} Sol;
typedef Sol * pSol;
 

typedef struct {
  int     size;
  int    *head;
  int    *link;
} Bucket;
typedef Bucket * pBucket;


/* prototypes */
int     loadMesh(pMesh );
int     loadSol(pSol );
int     saveSol(pSol );
int     scaleMesh(pMesh);
int     unscaleMesh(pMesh);
int     distance(pMesh ,pMesh ,pBucket,pSol, int);
int     distancequad(pMesh ,pMesh ,pSol,int);
int     initialization(pMesh ,pMesh);
int     initializationsphere(pMesh ,pMesh);
int     saveMesh(pMesh, int);
int     moveMesh(pMesh,pSol);
int     elasti1_3d(pMesh ,pSol );
pBucket newBucket_3d(pMesh ,int );


#endif