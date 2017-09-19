#include "warping.h"
extern Info   info;

/* initialization for envelope template */
int initialization( pMesh extmesh, pMesh intmesh ) {
  pTria   ptt;
  pPoint  ppt,ppt2;
  double  dist,app;
  int     k,i,l,countmax,nbpt;
  
  /* retrieve internal boundary */
  dist = FLT_MAX;
  
  for (k=1; k<=extmesh->nt; k++) {
    ptt = &extmesh->tria[k];
    
    for (i=0; i<3 ; i++) {
      ppt = &extmesh->point[ptt->v[i]];
      
      for (l=1; l<=intmesh->np; l++) {
        ppt2 = &intmesh->point[l];
        app = (ppt->c[0]-ppt2->c[0])*(ppt->c[0]-ppt2->c[0]) + (ppt->c[1]-ppt2->c[1])*(ppt->c[1]-ppt2->c[1]) + (ppt->c[2]-ppt2->c[2])*(ppt->c[2]-ppt2->c[2]);
        
        if (app < dist ){
          dist = app;
          nbpt = k;
        }
      }
    }
  }
  
  ptt = &extmesh->tria[nbpt];
  extmesh->ref = ptt->ref;

  /* initialize countmax */
  countmax = 0;
  for (k=1; k<=extmesh->nt; k++) {
    ptt = &extmesh->tria[k];
    if ( ptt->ref == extmesh->ref)   {
      for (i=0; i<3 ; i++) {
        ppt = &extmesh->point[ptt->v[i]];
        if (!(ppt->ref == 30)) countmax ++;
        ppt->ref = 30 ;
      }
    }
  }
  
  
  for (k=1; k<=extmesh->nt; k++) {
    ptt = &extmesh->tria[k];
    if ( ptt->ref != extmesh->ref){
      for (i=0; i<3 ; i++) {
        ppt = &extmesh->point[ptt->v[i]];
        ppt->ref = 25;
      }
      break;
    }
  }
  
  return(countmax);
}


/* initialization for spherical template */
int  initializationsphere(pMesh extmesh,pMesh intmesh){
  
  pTria   ptt;
  pPoint  ppt,ppt2;
  double  cal,dist,Ray,app;
  int     k,i,m,nbpt,l,countmax,nbr;
  
  /* retrieve internal reference */
  for (i=0; i<3; i++) {
    extmesh->min[i] =  FLT_MAX;
    extmesh->max[i] =  -FLT_MAX;
  }
  
  for (k=1; k<=extmesh->nt; k++) {
    ptt = &extmesh->tria[k];
    for (i=0; i<3 ; i++) {
      ppt = &extmesh->point[ptt->v[i]];
      
      for (l=0; l<3; l++) {
        if ( ppt->c[l] > extmesh->max[l] )  extmesh->max[l] = ppt->c[l];
        if ( ppt->c[l] < extmesh->min[l] )  extmesh->min[l] = ppt->c[l];
      }
    }
  }
  
  for (m=0; m<3; m++) extmesh->o[m] = 0.5*( extmesh->max[m] + extmesh->min[m]);
  
  dist = FLT_MAX;
  
  for (k=1; k<=extmesh->nt; k++) {
    ptt = &extmesh->tria[k];
    
    for (i=0; i<3 ; i++) {
      ppt = &extmesh->point[ptt->v[i]];
      app = (ppt->c[0]-extmesh->o[0])*(ppt->c[0]-extmesh->o[0]) + (ppt->c[1]-extmesh->o[1])*(ppt->c[1]-extmesh->o[1]) + (ppt->c[2]-extmesh->o[2])*(ppt->c[2]-extmesh->o[2]);
      if (app < dist ){
        dist = app;
        nbpt = k;
      }
    }
  }
  ptt = &extmesh->tria[nbpt];
  extmesh->ref = ptt->ref;
  
  /* compute center and ray extmesh */
  for (i=0; i<3; i++) {
    extmesh->min[i] =  FLT_MAX;
    extmesh->max[i] =  -FLT_MAX;
  }
  countmax = 0;
  for (k=1; k<=extmesh->nt; k++) {
    ptt = &extmesh->tria[k];
    if ( ptt->ref == extmesh->ref)   {
      for (i=0; i<3 ; i++) {
        ppt = &extmesh->point[ptt->v[i]];
        if (!(ppt->ref == 30)) countmax ++;
        ppt->ref = 30 ;
        for (l=0; l<3; l++) {
          if ( ppt->c[l] > extmesh->max[l] )  extmesh->max[l] = ppt->c[l];
          if ( ppt->c[l] < extmesh->min[l] )  extmesh->min[l] = ppt->c[l];
        }
      }
    }
  }
  
  for (m=0; m<3; m++) {
    extmesh->o[m] = 0.5*( extmesh->max[m] + extmesh->min[m]);
    extmesh->Ray[m] = extmesh->max[m] - extmesh->o[m];
  }
  
  /* compute center and ray intmesh */
  for (i=0; i<3; i++) {
    intmesh->min[i] =  FLT_MAX;
    intmesh->max[i] = -FLT_MAX;
  }
  for (k=1; k<=intmesh->np; k++) {
    ppt = &intmesh->point[k];
    for (i=0; i<3; i++) {
      if ( ppt->c[i] > intmesh->max[i] )    intmesh->max[i] = ppt->c[i];
      if ( ppt->c[i] < intmesh->min[i] )  intmesh->min[i] = ppt->c[i];
    }
  }
  for (m=0; m<3; m++) {
    intmesh->o[m] = 0.5*( intmesh->max[m] + intmesh->min[m] );
    intmesh->Ray[m] = intmesh->max[m] - intmesh->o[m];
  }
  fprintf(stdout,"     - ray intmesh %f %f %f \n",intmesh->Ray[0],intmesh->Ray[1],intmesh->Ray[2]);
  
  /*compute maximal distance between the intmesh points and the center*/
  dist = -FLT_MAX;
  for (k=1; k<=intmesh->np; k++) {
    ppt = &intmesh->point[k];
    cal = (ppt->c[0]-intmesh->o[0])*(ppt->c[0]-intmesh->o[0])+(ppt->c[1]-intmesh->o[1])*(ppt->c[1]-intmesh->o[1])+(ppt->c[2]-intmesh->o[2])*(ppt->c[2]-intmesh->o[2]);
    if(cal > dist ) {
      dist = cal ;
      nbr = k;
    }
  }
  /*align bounding box */
  Ray = sqrt(dist);
  for (k=1; k<=extmesh->np; k++) {
    ppt2 = &extmesh->point[k];
    for (l=0; l<3; l++) {
      ppt2->c[l] = (Ray /extmesh->Ray[2]) * (ppt2->c[l]-extmesh->o[l]) + intmesh->o[l];
    }
  }
  
  /* compute center and ray extmesh new */
  for (i=0; i<3; i++) {
    extmesh->min[i] =  FLT_MAX;
    extmesh->max[i] =  -FLT_MAX;
  }
  for (k=1; k<=extmesh->nt; k++) {
    ptt = &extmesh->tria[k];
    if ( ptt->ref == extmesh->ref)   {
      for (i=0; i<3; i++) {
        ppt = &extmesh->point[ptt->v[i]];
        for (l=0; l<3; l++) {
          if ( ppt->c[l] > extmesh->max[l] )  extmesh->max[l] = ppt->c[l];
          if ( ppt->c[l] < extmesh->min[l] )  extmesh->min[l] = ppt->c[l];
        }
      }
    }
  }
  for (m=0; m<3; m++) {
    extmesh->o[m] = 0.5*( extmesh->max[m] + extmesh->min[m] );
    extmesh->Ray[m] = extmesh->max[m] - extmesh->o[m];
  }
  fprintf(stdout,"     - ray extmesh new %f \n",extmesh->Ray[0]);
  
  /* find "contact" point */
  ppt2 = &intmesh->point[nbr];
  dist = FLT_MAX;
  for (k=1; k<=extmesh->nt; k++) {
    ptt = &extmesh->tria[k];
    if( ptt->ref == extmesh->ref) {
      for (l=0; l<3; l++) {
        ppt = &extmesh->point[ptt->v[l]];
        cal = (ppt->c[0]-ppt2->c[0])*(ppt->c[0]-ppt2->c[0]) + (ppt->c[1]-ppt2->c[1])*(ppt->c[1]-ppt2->c[1]) + (ppt->c[2]-ppt2->c[2])*(ppt->c[2]-ppt2->c[2]);
        if (cal < dist ){
          dist = cal;
          nbpt = ptt->v[l];
        }
      }
    }
  }
  for (k=1; k<=extmesh->nt; k++) {
    ptt = &extmesh->tria[k];
    if (ptt->ref == extmesh->ref)   {
      for (i=0; i<3 ; i++) {
        if ( ptt->v[i] == nbpt ) {
          ppt = &extmesh->point[ptt->v[i]];
          ppt->ref = 25;
        }
      }
    }
  }
  return(countmax);
}
