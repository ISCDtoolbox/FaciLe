#include "warping.h"

#define SCALING  1

extern Info    info;



/* just for translation */
int scaleMesh(pMesh mesh) {
  pPoint    ppt;
  double    dd,o[3];
  int       i,k;
  
  /* compute bounding box */
  for (i=0; i<mesh->dim; i++) {
    info.min[i] =  FLT_MAX;
    info.max[i] = -FLT_MAX;
  }
  
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (i=0; i<mesh->dim; i++) {
      if ( ppt->c[i] > info.max[i] )  info.max[i] = ppt->c[i];
      if ( ppt->c[i] < info.min[i] )  info.min[i] = ppt->c[i];
    }
  }
  info.delta = 0.0;
  for (i=0; i<mesh->dim; i++) {
    dd = fabs(info.max[i]-info.min[i]);
    if ( dd > info.delta )  info.delta = dd;
  }
  if ( info.delta < EPSD ) {
    fprintf(stdout,"  ## Unable to scale mesh\n");
    return(0);
  }
  
  for (i=0; i<mesh->dim; i++) {
    o[i] = 0.5 *( info.max[i] + info.min[i] );
  }
  
  /* normalize coordinates */
  
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (i=0; i<mesh->dim; i++){
      ppt->c[i] = SCALING*( ppt->c[i]-o[i]) + 0.5;
      assert(ppt->c[i]>0);
    }
  }
  
  return(1);
}

int unscaleMesh(pMesh mesh) {
  pPoint     ppt;
  int        k,i;
  
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (i=0; i<mesh->dim; i++)
      ppt->c[i] =  ppt->c[i] + info.min[i];
  }
  
  return(1);
}


int moveMesh(pMesh mesh,pSol sol) {
  pPoint    ppt;
  int       i,k;
  
  for (k=1; k<=mesh->np; k++){
    ppt = &mesh->point[k];
    
    for (i=0; i< mesh->dim; i++){
      ppt->c[i] = ppt->c[i] + sol->u[(mesh->dim)*(k-1)+i];
      
    }
  }
  return(1);
}
