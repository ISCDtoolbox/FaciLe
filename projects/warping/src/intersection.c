#include "warping.h"
extern Info   info;

#define LINKSIZ 800;


/* create bucket structure and store triangles */
pBucket newBucket_3d(pMesh mesh,int nmax) {
  pPoint   pa,pb,pc;
  pTria    pt;
  pBucket  bucket;
  double   dd,xmin,xmax,ymin,ymax,zmin,zmax;
  int      i,j,k,ic,p,jmax,jmin,imax,imin,kmin,kmax,s=0,linksize;
  
  linksize = LINKSIZ;
  
  /* memory alloc */
  bucket = (Bucket*)malloc(sizeof(Bucket));
  assert(bucket);
  bucket->size = nmax;
  bucket->head = (int*)calloc(nmax*nmax*nmax+1,sizeof(int));
  assert(bucket->head);
  bucket->link = (int*)calloc(linksize*(mesh->nt+1),sizeof(int));
  assert(bucket->link);
  
  /* insert triangles */
  dd = nmax;
  for (p=1; p<=mesh->nt; p++) {
    
    pt = &mesh->tria[p];
    pa = &mesh->point[pt->v[0]];
    pb = &mesh->point[pt->v[1]];
    pc = &mesh->point[pt->v[2]];
    
    xmin = LS_MIN(pa->c[0],LS_MIN(pb->c[0],pc->c[0]));
    xmax = LS_MAX(pa->c[0],LS_MAX(pb->c[0],pc->c[0]));
    ymin = LS_MIN(pa->c[1],LS_MIN(pb->c[1],pc->c[1]));
    ymax = LS_MAX(pa->c[1],LS_MAX(pb->c[1],pc->c[1]));
    zmin = LS_MIN(pa->c[2],LS_MIN(pb->c[2],pc->c[2]));
    zmax = LS_MAX(pa->c[2],LS_MAX(pb->c[2],pc->c[2]));
    
    
    imin = LS_MAX(0,(int)(dd * xmin)-1);
    jmin = LS_MAX(0,(int)(dd * ymin)-1);
    kmin = LS_MAX(0,(int)(dd * zmin)-1);
    imax = LS_MAX(0,(int)(dd * xmax)-1);
    jmax = LS_MAX(0,(int)(dd * ymax)-1);
    kmax = LS_MAX(0,(int)(dd * zmax)-1);
    
    for (k=kmin; k<=kmax; k++)
      for (j=jmin; j<=jmax; j++)
        for (i=imin; i<=imax; i++){
          
          ic = (k*nmax + j)*nmax + i;
          
          if ( !bucket->head[ic])
            bucket->head[ic] = p;
          else {
            s=0;
            while(bucket->link[linksize* bucket->head[ic]+s]) s++;
            bucket->link[linksize* bucket->head[ic]+s] = p;
            assert(s<linksize);
          }
        }
  }
  
  return(bucket);
}


static int RayTriaIntersection_3d( pMesh mesh, pTria pt1, double *pt,double *o,double *dm,double *q ){
  
  pPoint  ppta,pptb,pptc;
  double  a,b,abx,aby,abz,acx,acy,acz,A,B,r,C,D,E,F,G,H,det,dd,qapp[3];
  double  n[3];
  int     ier;
  
  ier=0;
  dm[0]=0;
  
  ppta = &mesh->point[pt1->v[0]];
  pptb = &mesh->point[pt1->v[1]];
  pptc = &mesh->point[pt1->v[2]];
  
  abx = pptb->c[0] - ppta->c[0];
  aby = pptb->c[1] - ppta->c[1];
  abz = pptb->c[2] - ppta->c[2];
  
  acx = pptc->c[0] - ppta->c[0];
  acy = pptc->c[1] - ppta->c[1];
  acz = pptc->c[2] - ppta->c[2];
  
  n[0] = aby*acz - abz*acy;
  n[1] = abz*acx - abx*acz;
  n[2] = abx*acy - aby*acx;
  
  dd = sqrt ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );
  
  /* A = dot ( n , PA) */
  A = n[0]*(ppta->c[0]-pt[0]) + n[1]*(ppta->c[1]-pt[1]) + n[2]*(ppta->c[2]-pt[2]);
  
  /* B = dot ( n , PO) */
  B = - (n[0]*(pt[0]-o[0]) + n[1]*(pt[1]-o[1]) + n[2]*(pt[2]-o[2]));
  
  if (fabs(B)>0) {
    /* get intersection ray/plane of the triangle Q */
    r = A/B;
    if ( ( r>= 0.0) && (r <=1)) {
      qapp[0] = pt[0] -r* (pt[0]-o[0]);
      qapp[1] = pt[1] -r* (pt[1]-o[1]);
      qapp[2] = pt[2] -r* (pt[2]-o[2]);
      
      /* G = -dot(cq,ab), H = -dot(cq,ac) */
      G = (pptc->c[0] - qapp[0])*abx + (pptc->c[1] - qapp[1])*aby + (pptc->c[2] - qapp[2])*abz;
      H = (pptc->c[0] - qapp[0])*acx + (pptc->c[1] - qapp[1])*acy + (pptc->c[2] - qapp[2])*acz;
      
      /* C = dot(aq, ab) - dot(cq,ab)
       D = dot(bq, ab) - dot(cq,ab)
       E = dot(aq, ac) - dot(cq,ac)
       F = dot(bq, ac) - dot(cq,ac) */
      C = -((ppta->c[0] - qapp[0])*abx + (ppta->c[1] - qapp[1])*aby + (ppta->c[2] - qapp[2])*abz) +G ;
      D = -((pptb->c[0] - qapp[0])*abx + (pptb->c[1] - qapp[1])*aby + (pptb->c[2] - qapp[2])*abz) +G ;
      E = -((ppta->c[0] - qapp[0])*acx + (ppta->c[1] - qapp[1])*acy + (ppta->c[2] - qapp[2])*acz) +H ;
      F = -((pptb->c[0] - qapp[0])*acx + (pptb->c[1] - qapp[1])*acy + (pptb->c[2] - qapp[2])*acz) +H ;
      
      det = C*F - E*D ;
      if (det != 0.0) {
        a = (F*G-D*H)/det;
        if (a>=0) {
          b = (-E*G + C*H)/det;
          /* case ray triangle intersect */
          if ( (b>=0.0)&&((1-a-b)>=0.0) ) {
            ier ++;
            dm[0] = (qapp[0]-pt[0])*(qapp[0]-pt[0]) + (qapp[1]-pt[1])*(qapp[1]-pt[1]) + (qapp[2]-pt[2])*(qapp[2]-pt[2]);
            q[0] = qapp[0];
            q[1] = qapp[1];
            q[2] = qapp[2];
            
          }
        }
      }
    }
  }
  
  return (ier);
  
}

/* check intersection between triangles in mesh and ray ab */
/* if intersection not found return 0 else store in q the intersection point (closer to a in case of multiple intersection ) */
static int RayMeshIntersection_3d(pMesh mesh,pBucket bucket,double *a, double *b, double *q) {
  
  pTria     pt1;
  double    dm[1],dd,dapp[1],qapp[3];
  int       i,j,k,ii,jj,kk,iib,jjb,kkb,ica,icb,ip,ip1,siz;
  int       imin,imax,jmin,jmax,kmin,kmax,ier,ier2,ic,linksize,s,iptria;
  
  linksize = LINKSIZ;
  siz = bucket->size;
  siz = 64;
  dd  = siz / (double)PRECI;
  
  for(i=0; i<3; i++){
    q[i]    = 0.0;
    qapp[i] = 0.0;
  }
  int cont1, cont2,cont3;
  cont1 = 0;
  cont2=0;
  cont3=0;
  dm[0] = 1e3; // in dm we store the minimal distance
  dapp[0]=0;
  ier=0;
  ier2=0;
  
  /* check for intersecting triangles in the cell containing a */
  ii = LS_MAX(0,(int)(dd * a[0])-1);
  jj = LS_MAX(0,(int)(dd * a[1])-1);
  kk = LS_MAX(0,(int)(dd * a[2])-1);
  ica = (kk*siz + jj)*siz + ii;
  
  iib = LS_MAX(0,(int)(dd * b[0])-1);
  jjb = LS_MAX(0,(int)(dd * b[1])-1);
  kkb = LS_MAX(0,(int)(dd * b[2])-1);
  icb = (kkb*siz + jjb)*siz + iib;
  
  
  if ( bucket->head[ica] ) {
    
    ip1 = bucket->head[ica];
    pt1 = &mesh->tria[ip1];
    ier2=RayTriaIntersection_3d(mesh,pt1,a,b,dapp,qapp);
    if(ier2) {
      if(dapp[0]<dm[0]){
        dm[0] = dapp[0];
       q[0] = qapp[0];
       q[1] = qapp[1];
       q[2] = qapp[2];
       ier = LS_MAX(ier2,ier);
        //return ier;
      }
    }
   
    s=0;
    iptria = (linksize*(bucket->head[ica]));
    while ( bucket->link[iptria+s] ) {
      ip = bucket->link[iptria+s];
      pt1 = &mesh->tria[ip];
      
      ier2=RayTriaIntersection_3d(mesh,pt1,a,b,dapp,qapp);
      if (ier2) {
        cont2++;
        if(dapp[0] < dm[0]){
        dm[0] = dapp[0];
        q[0] = qapp[0];
        q[1] = qapp[1];
        q[2] = qapp[2];
        ier = LS_MAX(ier,ier2);
        //return ier;
      }
      }
      s++;
    }
  }
  
  if(ica==icb) return ier;
  else{
    imin = LS_MIN(ii,iib);
    imax = LS_MAX(ii,iib);
    jmin = LS_MIN(jj,jjb);
    jmax = LS_MAX(jj,jjb);
    kmin = LS_MIN(kk,kkb);
    kmax = LS_MAX(kk,kkb);
    
    for (k=kmin; k<=kmax; k++)
      for (j=jmin; j<=jmax; j++)
        for (i=imin; i<=imax; i++){
          if ( ii == i && jj == j && kk == k )  continue;
          ic  = (k*siz + j)*siz + i;
          if ( bucket->head[ic] ) {
            ip1 = bucket->head[ic];
            pt1 = &mesh->tria[ip1];
            ier2=RayTriaIntersection_3d(mesh,pt1,a,b,dapp,qapp);
            if(ier2) {
              ip  = ip1;
              if(dapp[0] < dm[0]){
              dm[0]=dapp[0];
              q[0] = qapp[0];
              q[1] = qapp[1];
              q[2] = qapp[2];
              ier = LS_MAX(ier,ier2);
              //return ier;
              }
            }
            iptria = (linksize*(bucket->head[ic]));
            s=0;
            while ( bucket->link[iptria+s]) {
              ip = bucket->link[iptria+s];
              pt1 = &mesh->tria[ip];
              ier2=RayTriaIntersection_3d(mesh,pt1,a,b,dapp,qapp);
              if (ier2) {
                if(dapp[0] < dm[0]){
                dm[0] = dapp[0];
                q[0] = qapp[0];
                q[1] = qapp[1];
                q[2] = qapp[2];
                ier = LS_MAX(ier,ier2);
                //return ier;
                }
              }
              s++;
            }
          }
        }
  }
  
  return(ier);
}

int distance( pMesh extmesh, pMesh intmesh, pBucket bucket, pSol sol, int count) {
  
  pPoint  ppt;
  int     k,i,is,l,ier;
  double  b[3],q[3],a[3];
  char    stim[32];
  
  chrono(ON,&info.ctim[4]);
  if ( info.imprim ) fprintf(stdout,"  1.3 DISTANCE \n");
  for (k=1; k<=extmesh->np; k++) {
    ppt = &extmesh->point[k];
    is = 3*(k-1);
    if ( ppt->ref == 30) {
      for (i=0; i<3; i++)    b[i] = ppt->c[i] + sol->u[is+i];
      a[0] = ppt->c[0];
      a[1] = ppt->c[1];
      a[2] = ppt->c[2];
      ier=RayMeshIntersection_3d(intmesh,bucket,a,b,q);
      if(ier) {   /* if intersection found the solution is updated */
        ppt->ref=25;
        count ++;
        for (l=0; l<3; l++) {
          sol->u[is+l] = q[l]-ppt->c[l];
        }
      }
    }
  }
  chrono(OFF,&info.ctim[4]);
  printim(info.ctim[4].gdif,stim);
  if ( info.imprim ) fprintf(stdout,"     [Time: %s]\n",stim);
  return(count);
}
