#include "warping.h"

extern Info   info;

/* read mesh */
int loadMesh(pMesh mesh) {
  pPoint       ppt;
  pTria        pt1;
  pTetra       ptt;
  pQuad        pq;
  float        fp1,fp2,fp3;
  int          k,i,dof,inm;
  char        *ptr,data[256];

  strcpy(data,mesh->name);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
    }
  }
  else if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  if ( abs(info.imprim) > 3 )
    fprintf(stdout,"  -- READING DATA FILE %s\n",data);

  mesh->np = GmfStatKwd(inm,GmfVertices);
  mesh->na = GmfStatKwd(inm,GmfEdges);
  mesh->nt = GmfStatKwd(inm,GmfTriangles);
  mesh->nq = GmfStatKwd(inm,GmfQuadrilaterals);
  mesh->ne = GmfStatKwd(inm,GmfTetrahedra);

  if ( !mesh->np ) {
    fprintf(stdout,"  ** MISSING DATA\n");
    return(0);
  }
	mesh->npi = mesh->np;
	mesh->nai = mesh->na;
	mesh->nti = mesh->nt;
	mesh->nei = mesh->ne;

  /* memory alloc */
	dof = 1;
  mesh->point = (pPoint)calloc(dof*mesh->np+1,sizeof(Point));
  assert(mesh->point);
  if ( mesh->nt ) {
    mesh->tria  = (pTria)calloc(mesh->nt+1,sizeof(Tria));
    assert(mesh->tria);
  }
  if ( mesh->ne ) {
    mesh->tetra  = (pTetra)calloc(mesh->ne+1,sizeof(Tetra));
    assert(mesh->tetra);
  }
  if ( mesh->nq ) {
      mesh->quad  = (pQuad)calloc(mesh->nq+1,sizeof(Quad));
      assert(mesh->quad);
    }
 
 if(mesh->dim==2) fprintf(stdout,"  -- 2D CASE NON TREATED \n");
 assert(mesh->dim=3);
  
  /* read vertices */
    GmfGotoKwd(inm,GmfVertices);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( mesh->ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ppt->ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
        ppt->c[2] = fp3;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
    }
    /* read triangles and store barycenters */
    GmfGotoKwd(inm,GmfTriangles);
  for (k=1; k<=mesh->nt; k++) {
    pt1 = &mesh->tria[k];
    for (i=0; i<3; i++) {
    pt1->g[i] = 0.0;
    }
  }
    for (k=1; k<=mesh->nt; k++) {
      pt1 = &mesh->tria[k];
      GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
      for (i=0; i<3; i++) {
        ppt = &mesh->point[pt1->v[i]];
        pt1->g[0] += ppt->c[0];
        pt1->g[1] += ppt->c[1];
        pt1->g[2] += ppt->c[2];
      }
      pt1->g[0] /= 3.0;
      pt1->g[1] /= 3.0;
      pt1->g[2] /= 3.0;

    }
    /* read tetrahedra */
    GmfGotoKwd(inm,GmfTetrahedra);
    for (k=1; k<=mesh->ne; k++) {
      ptt = &mesh->tetra[k];
      GmfGetLin(inm,GmfTetrahedra,&ptt->v[0],&ptt->v[1],&ptt->v[2],&ptt->v[3],&ptt->ref);
      
    }
      /* read quadrilaterals */
      GmfGotoKwd(inm,GmfQuadrilaterals);
      for (k=1; k<=mesh->nq; k++) {
          pq = &mesh->quad[k];
          GmfGetLin(inm,GmfQuadrilaterals,&pq->v[0],&pq->v[1],&pq->v[2],&pq->v[3],&pq->ref);
      }

    if ( info.imprim ) {
      fprintf(stdout,"  %%%% NUMBER OF VERTICES   %8d\n",mesh->np);
      if ( mesh->nt )  fprintf(stdout,"  %%%% NUMBER OF TRIANGLES  %8d\n",mesh->nt);
      if ( mesh->ne )  fprintf(stdout,"  %%%% NUMBER OF TETRAHEDRA %8d\n",mesh->ne);
      if ( mesh->nq )  fprintf(stdout,"  %%%% NUMBER OF QUADRILATERALS %8d\n",mesh->nq);
    }

  GmfCloseMesh(inm);
  return(1);
}

int saveMesh(pMesh mesh,int n) {
  pPoint       ppt;
  pTria        pt1;
  pTetra       pt;
  int          k,inm;
  char        *ptr,data[128],suff[10];

  mesh->ver = GmfDouble;
  strcpy(data,mesh->nameout);
  
  if( info.debug == 1) {
    ptr = strstr(data,".mesh");
    if ( !ptr ) {
    sprintf(data, "%s.%03d",data,n);
      strcpy(suff,".mesh");
      strcat(data,suff);
    }
    else {
      *ptr = '\0';
      sprintf(data, "%s.%03d",data,n);
      strcpy(suff,".mesh");
      strcat(data,suff);
    }
  }
  else {
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".mesh");
    if( !(inm = GmfOpenMesh(data, GmfWrite, mesh->ver,mesh->dim)) ) {
        fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
        return(0);
    }
  }
  else {
    *ptr = '\0';
    strcat(data,".mesh");
    }
  }
    if( !(inm = GmfOpenMesh(data, GmfWrite, mesh->ver,mesh->dim)) ){
    fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
    return(0);
   }
  
  fprintf(stdout,"  %%%% %s OPENED\n",data);
  
  GmfSetKwd(inm,GmfVertices,mesh->np);
      for(k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
     
      GmfSetLin(inm,GmfVertices,ppt->c[0],ppt->c[1], ppt->c[2],ppt->ref);
  }
  /* write triangles */
  GmfSetKwd(inm,GmfTriangles,mesh->nt);
  for (k=1; k<=mesh->nt; k++) {
    pt1 = &mesh->tria[k];
    GmfSetLin(inm,GmfTriangles,pt1->v[0],pt1->v[1],pt1->v[2],pt1->ref);
  }
  
  /* write tetrahedra */
  GmfSetKwd(inm,GmfTetrahedra,mesh->ne);
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->v[0] )  continue;
    GmfSetLin(inm,GmfTetrahedra,pt->v[0],pt->v[1],pt->v[2],pt->v[3],pt->ref);
  }

GmfCloseMesh(inm);
  return(1);
}

int loadSol(pSol sol) {
  float       buf[GmfMaxTyp];
  double      bufd[GmfMaxTyp];
  int         i,k,type,inm,typtab[GmfMaxTyp],offset;
  char       *ptr,data[128];

	if ( !sol->namein )  return(-1);
  strcpy(data,sol->namein);
  ptr = strstr(data,".sol");
  if ( !ptr ) {
    strcat(data,".sol");
    if ( !(inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim)) ) {
      ptr = strstr(data,".sol");
      *ptr = '\0';
      strcat(data,".solb");
      if ( !(inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim)) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
    }
  }
  else if ( !(inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim)) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
    return(0);
  }

  if ( !(inm = GmfOpenMesh(data, GmfRead, &sol->ver,&sol->dim)) ) {
    return(-1);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  if ( abs(info.imprim) > 3 )
    fprintf(stdout,"  -- READING DATA FILE %s\n",data);

  sol->np = GmfStatKwd(inm,GmfSolAtVertices,&type,&offset,&typtab);
  if ( !sol->np || typtab[0] != 2 )  return(-1);

  /* alloc */
  sol->u  = (double*)calloc(sol->dim*sol->np,sizeof(double));
  assert(sol->u);

  /* read mesh solutions */
  GmfGotoKwd(inm,GmfSolAtVertices);
  if ( sol->ver == GmfFloat ) {
    for (k=0; k<sol->np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,&buf);
    for (i=0; i<sol->dim; i++)      
        sol->u[sol->dim*k+i] = buf[i];        
    }
  }
  else {
    for (k=0; k<sol->np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,bufd);
    for (i=0; i<sol->dim; i++)      
        sol->u[sol->dim*k+i] = bufd[i];         
    }
  }

  GmfCloseMesh(inm);
  return(1);
}

int saveSol(pSol sol) {
  double       dbuf[GmfMaxTyp];
  float        fbuf[GmfMaxTyp];
  int          k,ia,i,inm,type,typtab[GmfMaxTyp];
  char        *ptr,data[128];

  strcpy(data,sol->nameout);
  ptr = strstr(data,".mesh");
  if ( ptr )  {
    *ptr = '\0';
    strcat(data,".solb");
  }
  else {
    ptr = strstr(data,".sol");
    if ( !ptr )  strcat(data,".sol");
  }

  if ( !(inm = GmfOpenMesh(data,GmfWrite,sol->ver,sol->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  type = 1;
  typtab[0] = GmfVec;

  /* write sol */
  GmfSetKwd(inm,GmfSolAtVertices,sol->np+sol->na,type,typtab);
  if ( sol->ver == GmfFloat ) {
    for (k=0; k<sol->np+sol->na; k++) {
      ia = sol->dim*k;
      for (i=0; i<sol->dim; i++)      
        fbuf[i] = sol->u[ia+i];
      GmfSetLin(inm,GmfSolAtVertices,fbuf);
    }
  }
  else {
    for (k=0; k<sol->np+sol->na; k++) {
      ia = sol->dim*k;
      for (i=0; i<sol->dim; i++)      
        dbuf[i] = sol->u[ia+i];
      GmfSetLin(inm,GmfSolAtVertices,dbuf);
    }
  }

  GmfCloseMesh(inm);
  return(1);
}
