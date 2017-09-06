#include "wrapping.h"
#include "compil.date"
#define BUCKSIZ 64

Info  info;

static void excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
    case SIGABRT:
      fprintf(stdout,"  Abnormal stop\n");  break;
    case SIGBUS:
      fprintf(stdout,"  Code error...\n");  break;
    case SIGFPE:
      fprintf(stdout,"  Floating-point exception\n"); break;
    case SIGILL:
      fprintf(stdout,"  Illegal instruction\n"); break;
    case SIGSEGV:
      fprintf(stdout,"  Segmentation fault.\n");  break;
    case SIGTERM:
    case SIGINT:
      fprintf(stdout,"  Programm killed.\n");  break;
  }
  fprintf(stdout," No data file saved.\n");
  exit(1);
}

static void usage(char *prog) {
  fprintf(stdout,"\n usage: %s [-v[n]] [-h] [opts..] filein[.mesh]\n",prog);
  
  fprintf(stdout,"\n usage: %s [opts..] template[.mesh] source[.mesh] \n",prog);
  
  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"-h                Print this message \n");
  fprintf(stdout,"-scale            scale template mesh (only for default template) \n");
  
  
  fprintf(stdout,"\n**  File specifications\n");
  fprintf(stdout,"-t template[.mesh]   input triangulation to be deformed (default %s) \n","sphere");
  fprintf(stdout,"-s source[.mesh]     input triangulation to be warped\n");
  
  fprintf(stdout,"\n** Parameters\n");
  fprintf(stdout,"-nit n            iterations (default %d)\n",MAXIT);
  
  exit(1);
}

static int parsar(int argc,char *argv[],pMesh extmesh, pMesh intmesh,pSol sol) {
  int      i;
  char    *ptr;
  
  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
        case 'h':  /* on-line help */
        case '?':
          usage(argv[0]);
          break;
        case 'c':
          if ( !strcmp(argv[i],"-cpu") ) {
            ++i;
            if ( i == argc )
              fprintf(stderr,"  ** missing argument\n");
            else {
              if ( isdigit(argv[i][0]) )
                info.ncpu = atoi(argv[i]);
              else
                fprintf(stderr,"  ** argument [%s] discarded, use default: %d\n",argv[i],info.ncpu);
            }
          }
          break;
        case 'd':  /* debug */
          info.ddebug = 1;
          break;
        case 'e':
          if ( !strcmp(argv[i],"-err") ) {
            ++i;
            if ( isdigit(argv[i][0]) )
              sol->err = strtod(argv[i],NULL);
            else
              --i;
          }
          break;
        case 't':
          if ( !strcmp(argv[i],"-t") ) {
            ++i;
            extmesh->name = argv[i];
            info.scale = 0;
          }
          break;
        case 'n':
          if ( !strcmp(argv[i],"-nit") ) {
            ++i;
            if ( i < argc && isdigit(argv[i][0]) )
              info.nit = atoi(argv[i]);
            else
              --i;
          }
          break;
          
      case 's':
          
//          if ( !strcmp(argv[i],"-scale") ){
//            ++i;
//            info.scale = 1;
//          }
          if ( !strcmp(argv[i],"-s") ){
            ++i;
            intmesh->name = argv[i];
          }

          break;
          
//        case 't':
//          if ( ++i < argc ) {
//            if ( isdigit(argv[i][0]) )
//              info.typ = atoi(argv[i]);
//            else {
//              fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
//              usage(argv[0]);
//              i--;
//            }
//          }
//          else {
//            fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
//            usage(argv[0]);
//          }
//          break;
          
        case 'v':
          if ( ++i < argc ) {
            if ( argv[i][0] == '-' || isdigit(argv[i][0]) )
              info.imprim = atoi(argv[i]);
            else
              i--;
          }
          else {
            fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
            usage(argv[0]);
          }
          break;
        default:
          fprintf(stderr,"  Unrecognized option %s\n",argv[i]);
          usage(argv[0]);
      }
    }
    else {
      if ( extmesh->name == NULL ) {
        extmesh->name = argv[i];
        if ( info.imprim == -99 )  info.imprim = 5;
      }
      else if (intmesh->name == NULL) {
        intmesh->name = argv[i];
        if ( info.imprim == -99 )  info.imprim = 5;
      }
      else if (extmesh->ref==0) {
        extmesh->ref = atoi(argv[i]);
        if ( info.imprim == -99 )  info.imprim = 5;
      }
      
      else {
        fprintf(stdout,"  Argument %s ignored\n",argv[i]);
        usage(argv[0]);
      }
    }
    i++;
  }
  
  /* check params */
//  if ( info.imprim == -99 ) {
//    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
//    fflush(stdin);
//    fscanf(stdin,"%d",&i);
//    info.imprim = i;
//  }
  
  if ( extmesh->name == NULL ) {
    extmesh->name = (char *)calloc(128,sizeof(char));
    assert(extmesh->name);
    fprintf(stdout,"  -- EXTMESH BASENAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",extmesh->name);
  }
  
  if ( intmesh->name == NULL ) {
    intmesh->name = (char *)calloc(128,sizeof(char));
    assert(intmesh->name);
    fprintf(stdout,"  -- INTMESH BASENAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",intmesh->name);
  }
  
  if ( !sol->nameout ) {
    sol->nameout = (char *)calloc(128,sizeof(char));
    assert(sol->nameout);
    strcpy(sol->nameout,extmesh->name);
    ptr = strstr(sol->nameout,".mesh");
    if ( ptr ) *ptr = '\0';
    strcat(sol->nameout,".0.sol");
  }
  
  return(1);
}

static int parsop(pMesh mesh,pSol sol) {
  Cl         *pcl;
  Mat        *pm;
  float       fp1,fp2;
  int         i,j,ncld,ret;
  char       *ptr,buf[256],data[256];
  FILE       *in;
  
  strcpy(data,mesh->name);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".elas");
  in = fopen(data,"r");
  if ( !in ) {
    sprintf(data,"%s","DEFAULT.elas");
    in = fopen(data,"r");
    if ( !in )  return(1);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);
  
  /* read parameters */
  sol->nbcl = 0;
  while ( !feof(in) ) {
    ret = fscanf(in,"%s",data);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(data); i++) data[i] = tolower(data[i]);
    
    /* check for condition type */
    if ( !strcmp(data,"dirichlet")  || !strcmp(data,"load") ) {
      fscanf(in,"%d",&ncld);
      for (i=sol->nbcl; i<sol->nbcl+ncld; i++) {
        pcl = &sol->cl[i];
        if ( !strcmp(data,"load") )             pcl->typ = Load;
        else  if ( !strcmp(data,"dirichlet") )  pcl->typ = Dirichlet;
        else {
          fprintf(stdout,"  %%%% Unknown condition: %s\n",data);
          continue;
        }
        
        /* check for entity */
        fscanf(in,"%d %s ",&pcl->ref,buf);
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        fscanf(in,"%c",&pcl->att);
        pcl->att = tolower(pcl->att);
        if ( (pcl->typ == Dirichlet) && (pcl->att != 'v' && pcl->att != 'f') ) {
          fprintf(stdout,"  %%%% Wrong format: %s\n",buf);
          continue;
        }
        else if ( (pcl->typ == Load) && (pcl->att != 'v' && pcl->att != 'f' && pcl->att != 'n') ) {
          fprintf(stdout,"  %%%% Wrong format: %s\n",buf);
          continue;
        }
        if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") )          pcl->elt = LS_Ver;
        else if ( !strcmp(buf,"edges") || !strcmp(buf,"edge") )          pcl->elt = LS_Edg;
        else if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") )  pcl->elt = LS_Tri;
        
        if ( pcl->att != 'f' && pcl->att != 'n' ) {
          for (j=0; j<mesh->dim; j++) {
            fscanf(in,"%f ",&fp1);
            pcl->u[j] = fp1;
          }
        }
        else if ( pcl->att == 'n' ) {
          fscanf(in,"%f ",&fp1);
          pcl->u[0] = fp1;
        }
      }
      sol->nbcl += ncld;
    }
    /* gravity or body force */
    else if ( !strcmp(data,"gravity") ) {
      info.load |= (1 << 0);
      for (j=0; j<mesh->dim; j++) {
        fscanf(in,"%f ",&fp1);
        info.gr[j] = fp1;
      }
    }
    else if ( !strcmp(data,"lame") ) {
      fscanf(in,"%d",&ncld);
      assert(ncld <= LS_MAT);
      sol->nmat = ncld;
      for (i=0; i<ncld; i++) {
        pm = &sol->mat[i];
        fscanf(in,"%d %f %f\n",&pm->ref,&fp1,&fp2);
        pm->lambda = fp1;
        pm->mu     = fp2;
      }
    }
    else if ( !strcmp(data,"youngpoisson") ) {
      fscanf(in,"%d",&ncld);
      sol->nmat = ncld;
      for (i=0; i<ncld; i++) {
        pm = &sol->mat[i];
        fscanf(in,"%d %f %f\n",&pm->ref,&fp1,&fp2);
        pm->lambda = (fp1 * fp2) / ((1.0+fp2) * (1.0-2.0*fp2));
        pm->mu     = fp1 / (2.0*( 1.0+fp2));
      }
    }
  }
  fclose(in);
  
  for (i=0; i<sol->nbcl; i++) {
    pcl = &sol->cl[i];
    sol->cltyp |= pcl->elt;
  }
  return(1);
}

static void endcod() {
  char   stim[32];
  
  chrono(OFF,&info.ctim[0]);
  printim(info.ctim[0].gdif,stim);
  fprintf(stdout,"\n   ELAPSED TIME  %s\n",stim);
}

int main(int argc,char **argv) {
  Mesh     intmesh;
  Mesh     extmesh;
  Mesh     inmesh;
  Sol      sol;
  pBucket  bucket;
  int      ier,count,countmax,count1,count2,count3,stop;
  char     stim[32],ch;
  
  fprintf(stdout,"  -- WRAPPING, Release %s (%s) \n",LS_VER,LS_REL);
  fprintf(stdout,"     %s\n",LS_CPY);
  fprintf(stdout,"    %s\n",COMPIL);
  
  /* trap exceptions */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);
  signal(SIGBUS,excfun);
  atexit(endcod);
  
  tminit(info.ctim,TIMEMAX);
  chrono(ON,&info.ctim[0]);
  
  /* default values */
  memset(&intmesh,0,sizeof(Mesh));
  memset(&extmesh,0,sizeof(Mesh));
  memset(&inmesh,0,sizeof(Mesh));
  memset(&sol,0,sizeof(Sol));
  sol.cl  = (Cl*)calloc(LS_CL,sizeof(Cl));
  sol.mat = (Mat*)calloc(LS_MAT,sizeof(Mat));
  sol.nit = LS_MAXIT;
  
  info.imprim = 0;
  info.ddebug = 0;
  info.ncpu   = 1;
  info.zip    = 0;
  info.typ    = P1;
  info.scale  = 1; // scale extmesh by default
  info.nit    = MAXIT;
  /* default template */
  extmesh.name = "sphere";
  
  /* command line */
  if ( !parsar(argc,argv,&extmesh,&intmesh,&sol) )  return(1);
  
  /* load data */
  if ( info.imprim )   fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&info.ctim[1]);
  
  if ( !loadMesh(&extmesh) )  return(1);
  if ( !loadMesh(&intmesh) )  return(1);
  if(extmesh.ne==0) fprintf(stdout,"  -- Warning --- no tetrahedra found \n");
  assert(extmesh.ne>0);
  
  ier = loadSol(&sol);
  if ( !ier ) return(1);
  else if ( ier > 0 && sol.np != extmesh.np ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    free(sol.u0);
    sol.np = sol.ne = 0;
  }
  if ( !sol.dim )  sol.dim = extmesh.dim;
  if ( !parsop(&extmesh,&sol) )  return(1);
  chrono(OFF,&info.ctim[1]);
  printim(info.ctim[1].gdif,stim);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);
  
  chrono(ON,&info.ctim[2]);
  fprintf(stdout,"\n  %s\n   MODULE WARPING-ISCD : %s (%s)\n  %s\n",LS_STR,LS_VER,LS_REL,LS_STR);
  if ( info.imprim )   fprintf(stdout,"\n  -- INITIALISATION\n");
  
  bucket = newBucket_3d(&intmesh,BUCKSIZ);
  memset(sol.u,0,sol.dim*sol.np*sizeof(double));
  /* Initialisation of the algorithm */
  count = 1;
  count2 = 0;
  count3 = 0;
  if( info.scale) { fprintf(stdout,"\n  -- Using default spherical template ..  \n");
    countmax = initializationsphere(&extmesh,&intmesh);}
  else {
    fprintf(stdout,"\n  -- Using envelope template ..  \n");
    countmax = initialization(&extmesh,&intmesh);
  }
  if (!saveMesh(&extmesh)) return (1);
  //if( !copyMesh(&extmesh,&inmesh)) return (1); // copy and shrink extmesh
  intmesh.ref = 1;
  if ( info.imprim ) fprintf(stdout,"\n  -- PHASE 1\n");
  stop = 1;
  
  while( count <= 0.8*countmax && count3 <=info.nit ) {
    sol.err = LS_RES;
    if (!elasti1_3d(&extmesh,&sol) ) break;
    count1 = count;
    count =  distance(&extmesh,&intmesh,bucket,&sol,count);
    if (count1 < count ) count2 = 0;
    else count2++;
    if (!moveMesh(&extmesh,&sol)) return (1);
    /* reset sol at 0 */
    memset(sol.u,0,sol.dim*sol.np*sizeof(double));
    count3 ++;
   if ( info.imprim ) fprintf(stdout,"     - count = %d, countmax = %d, count2 = %d, count3 = %d, \n",count,countmax,count2,count3);
//    if (count3 >= 10 ) {
//      if (count3%3 == 0 ) {
//        fprintf(stdout,"  -- Continue? Y/n   \n");
//        scanf("%c",&ch);
//        if(ch=='n') break;
//      }
//    }

    
  }
  chrono(OFF,&info.ctim[2]);
  printim(info.ctim[2].gdif,stim);
  
  if ( info.imprim ) fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);
  
//  if ( info.imprim )   fprintf(stdout,"\n  -- PHASE 2\n");
//  
//  /* inflate a small mesh at the interior */
//  /* elasticity final */
//  
//  /* update references */
//  for (k=1; k<=extmesh.np; k++) {
//    ppt = &extmesh.point[k];
//    if  (ppt->ref == 25 )  {
//      ppti = &inmesh.point[k];
//      ppti->ref = 25;
//      /* update sol */
//      for (l=0; l<3; l++) {
//        sol.u[(extmesh.dim)*(k-1)+l] = (ppt->c[l]-ppti->c[l])/10;
//      }
//    }
//  }
//  
//  /* kill load conditions */
//  for (i=0; i<sol.nbcl; i++) {
//    pcl = &sol.cl[i];
//    if ( pcl->typ == Load )  pcl->typ = None;
//  }
//  for ( k=1; k<=10; k++) {
//    fprintf(stdout,"\n  -- PHASE 2 %d \n",k);
//    if (!elasti1_3d(&inmesh,&sol) )  return(1);
//    if (!moveMesh(&inmesh,&sol)) return (1);
//  }
//  if (!saveSurf(&inmesh,extmesh.name)) return (1);
 
  
  /* save file */
  if ( info.imprim ) fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",extmesh.name);
  if (!saveMesh(&extmesh)) return (1);
  chrono(ON,&info.ctim[1]);
  free(extmesh.tetra);
  if ( !info.zip ) free(extmesh.point);
  free(intmesh.tetra);
  if ( !info.zip ) free(intmesh.point);
  free(bucket->head);
  free(bucket->link);
  free(bucket);
  
  /* free mem */
  if ( sol.u )  free(sol.u);
  chrono(OFF,&info.ctim[0]);
  fprintf(stdout,"\n  %s\n   END OF MODULE WARPING. %s \n  %s\n",LS_STR,stim,LS_STR);
  return(0);
}
