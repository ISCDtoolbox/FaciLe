#include "wrapping.h"

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





int copyMesh(pMesh mesh,pMesh cmesh) {
  pPoint       ppt,pptc;
  pEdge        pa,pac;
  pTria        pt1,ptc,pt;
  pTetra       ptt,pttc;
  int          k,dof,l,m,i;
  char        *ptr,data[128];
  
  cmesh->dim = mesh->dim;
  cmesh->np = mesh->np;
  cmesh->na = mesh->na;
  cmesh->nt = mesh->nt;
  cmesh->ne = mesh->ne;
  cmesh->ver = mesh->ver;
  
  cmesh->name = (char *)calloc(128,sizeof(char));
  strcpy(data,mesh->name);
  ptr = strstr(data,".mesh");
  strcat(data,".int.mesh");
  strcpy(cmesh->name,data);
  
  for (i=0; i<mesh->dim; i++) {
    mesh->min[i] =  FLT_MAX;
    mesh->max[i] =  -FLT_MAX;
  }
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( pt->ref == mesh->ref)   {
      for (i=0; i<3 ; i++) {
        ppt = &mesh->point[pt->v[i]];
        for (l=0; l< mesh->dim; l++) {
          if ( ppt->c[l] > mesh->max[l] )  mesh->max[l] = ppt->c[l];
          if ( ppt->c[l] < mesh->min[l] )  mesh->min[l] = ppt->c[l];
        }
      }
    }
  }
  for (m=0; m<mesh->dim; m++) {
    mesh->o[m] = 0.5*( mesh->max[m] + mesh->min[m] );
    mesh->Ray[m] = mesh->max[0] - mesh->o[0];
  }
  if ( !cmesh->np ) {
    fprintf(stdout,"  ** MISSING DATA\n");
    return(0);
  }
  //memory alloc
  dof = info.typ == P2 ? 9 : 1;  // bound on number of nodes
  cmesh->point = (pPoint)calloc(dof*cmesh->np+1,sizeof(Point));
  assert(cmesh->point);
  if ( cmesh->nt ) {
    cmesh->tria  = (pTria)calloc(cmesh->nt+1,sizeof(Tria));
    assert(cmesh->tria);
  }
  if ( cmesh->ne ) {
    cmesh->tetra  = (pTetra)calloc(cmesh->ne+1,sizeof(Tetra));
    assert(cmesh->tetra);
  }
  
  /*copy and shrink mesh vertices */
  if ( cmesh->dim == 2 ) {
    for (k=1; k<=cmesh->np; k++) {
      ppt = &mesh->point[k];
      pptc = &cmesh->point[k];
      pptc->c[0] = mesh->o[0] + 0.1 *(ppt->c[0]-mesh->o[0]);
      pptc->c[1] = mesh->o[1] + 0.1 *(ppt->c[1]-mesh->o[1]);
    }
    
    /*copy mesh edges */
    if ( cmesh->na ) {
      cmesh->edge  = (pEdge)calloc(cmesh->na+1,sizeof(Edge));
      assert(cmesh->edge);
      
      for (k=1; k<=cmesh->na; k++) {
        pa = &cmesh->edge[k];
        pac = &cmesh->edge[k];
        pac->v[0] = pa->v[0];
        pac->v[1] = pac->v[1];
      }
    }
    /* copy mesh triangles */
    for (k=1; k<=cmesh->nt; k++) {
      pt1 = &mesh->tria[k];
      ptc = &cmesh->tria[k];
      ptc->v[0] = pt1->v[0];
      ptc->v[1] = pt1->v[1];
    }
  }
  else {
    /*copy and shrink mesh vertices */
    for (k=1; k<=cmesh->np; k++) {
      ppt = &mesh->point[k];
      pptc = &cmesh->point[k];
      pptc->ref = ppt->ref;
      pptc->c[0] = mesh->o[0] + 0.1 *(ppt->c[0]-mesh->o[0]);
      pptc->c[1] = mesh->o[1] + 0.1 *(ppt->c[1]-mesh->o[1]);
      pptc->c[2] = mesh->o[2] + 0.1 *(ppt->c[2]-mesh->o[2]);
    }
    /* copy mesh triangles */
    for (k=1; k<=mesh->nt; k++) {
      pt1 = &mesh->tria[k];
      ptc = &cmesh->tria[k];
      ptc->ref = pt1->ref;
      ptc->v[0] = pt1->v[0];
      ptc->v[1] = pt1->v[1];
      ptc->v[2] = pt1->v[2];
    }
    /*copy  tetrahedra */
    for (k=1; k<=cmesh->ne; k++) {
      ptt = &mesh->tetra[k];
      pttc = &cmesh->tetra[k];
      pttc->ref = ptt->ref;
      pttc->v[0] = ptt->v[0];
      pttc->v[1] = ptt->v[1];
      pttc->v[2] = ptt->v[2];
      pttc->v[3] = ptt->v[3];
    }
  }
  return(1);
  
}







