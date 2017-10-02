#include "warping.h"
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
  fprintf(stdout,"-p                Print execution steps  \n");
  fprintf(stdout,"-debug            Print intermediary meshes \n");
  
  
  fprintf(stdout,"\n**  File specifications\n");
  fprintf(stdout,"-t template[.mesh]   input triangulation to be deformed (default %s) \n","sphere");
  fprintf(stdout,"-s source[.mesh]     input triangulation to be warped\n");
  fprintf(stdout,"-o output[.mesh]     output filename \n");
  
  fprintf(stdout,"\n** Parameters\n");
  fprintf(stdout,"-nit n            iterations (default %d)\n",MAXIT);
  fprintf(stdout,"-load  l          pressure magnitude (default %d)\n",LOAD);
  fprintf(stdout,"-lame  l          Lame coefficients  (default lambda = %e mu = %e)\n",LS_E * LS_NU /( (1 + LS_NU)*(1-2*LS_NU)),LS_E /( 2*(1 + LS_NU)));
  fprintf(stdout,"-yp  l            Young and Poisson coefficients  (default nu = %e E = %d)\n",LS_NU,LS_E);
  
  exit(1);
}

static int parsar(int argc,char *argv[],pMesh extmesh, pMesh intmesh,pSol sol) {
  int      i;
  char    *ptr;
  
  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
          
        case 'h':
          usage(argv[0]);
          
          break;
        case 'd':
          if ( !strcmp(argv[i],"-debug") ) {
              info.debug = 1;
           }
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
          if ( !strcmp(argv[i],"-s") ){
            ++i;
            intmesh->name = argv[i];
          }
          break;
          
        case 'p':
          if ( !strcmp(argv[i],"-p") ){
            info.imprim = 1;
          }
          break;
        case 'o':
          if ( !strcmp(argv[i],"-o") ){
            ++i;
            extmesh->nameout = argv[i];
          }
          break;
          
        case 'l':
          if ( !strcmp(argv[i],"-load") ) {
            ++i;
            if ( i < argc && isdigit(argv[i][0]) )
              info.load = strtod(argv[i],NULL);
            else
              --i;
          }
          else if ( !strcmp(argv[i],"-lame") ) {
            ++i;
            if ( i < argc && isdigit(argv[i][0]) )
              info.lambda = strtod(argv[i],NULL);
            else
              --i;
            ++i;
            if ( i < argc && isdigit(argv[i][0]) )
              info.mu = strtod(argv[i],NULL);
            else
              --i;
          }
          break;
          
        case 'y':
          if ( !strcmp(argv[i],"-yp") ) {
            ++i;
            if ( i < argc && isdigit(argv[i][0]) )
              info.nu = strtod(argv[i],NULL);
            else
              --i;
            ++i;
            if ( i < argc && isdigit(argv[i][0]) )
              info.e = strtod(argv[i],NULL);
            else
              --i;
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
      }
      else if (intmesh->name == NULL) {
        intmesh->name = argv[i];
      }
      else if (extmesh->ref==0) {
        extmesh->ref = atoi(argv[i]);
      }
      else {
        fprintf(stdout,"  Argument %s ignored\n",argv[i]);
        usage(argv[0]);
      }
    }
    i++;
  }
  
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
  
  if ( !extmesh->nameout ) {
    extmesh->nameout = (char *)calloc(128,sizeof(char));
    assert(extmesh->nameout);
    strcpy(extmesh->nameout,extmesh->name);
    ptr = strstr(extmesh->nameout,".mesh");
    if ( ptr ) *ptr = '\0';
    strcat(extmesh->nameout,".d.mesh");
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
  int      ier,count,countmax,count1,count2,count3,stop,it,ret;
  char     stim[32],str1[200],str2[200],*ptr;
  
  fprintf(stdout,"  -- WARPING, Release %s (%s) \n",LS_VER,LS_REL);
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
  sol.nit = LS_MAXIT;
  
  info.imprim = 0;
  info.debug  = 0;
  info.nit    = MAXIT;
  info.load   = LOAD;
  info.lambda = LS_E * LS_NU /( (1 + LS_NU)*(1-2*LS_NU));
  info.mu     = LS_E /( 2*(1 + LS_NU));
  info.nu     = LS_NU;
  info.e      = LS_E;
  /* default template */
  extmesh.name = "sphere";
  
  /* command line */
  if ( !parsar(argc,argv,&extmesh,&intmesh,&sol) )  return(1);
  
  if ( info.imprim )  { fprintf(stdout,"  \n -- Lame coefficients lambda = %e mu = %e \n",info.lambda,info.mu);
    fprintf(stdout,"  \n -- Young and Poisson coefficients nu = %e E = %e \n",info.nu,info.e);
    fprintf(stdout,"  \n -- Load = %e \n",info.load);
  }
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
  it = 0;
  strcpy(str1,extmesh.name);
  ptr = strstr(str1,".mesh");
  if (! ptr) {
  strcpy(str2,"sphere");
  ret = strcmp(str1, str2);
  }
  else { strcpy(str2,"sphere.mesh");
    ret = strcmp(str1, str2);
  }

  if( !ret ) { fprintf(stdout,"\n  -- Using default spherical template ..  \n");
    countmax = initializationsphere(&extmesh,&intmesh);
  }
  else {
    fprintf(stdout,"\n  -- Using envelope template ..  \n");
    countmax = initialization(&extmesh,&intmesh);
  }
  intmesh.ref = 1;
  
  if ( info.imprim ) {
    fprintf(stdout,"\n  -- Template internal reference ref = %d ..  \n", extmesh.ref);
    fprintf(stdout,"\n  -- PHASE 1\n");
  }
  stop = 1;
  
  if (info.debug) {
    it++;
    saveMesh(&extmesh,it);
  }

  
//  while( count <= 0.9*countmax && count3 <=info.nit ) {
  while( count3 <=info.nit ) {
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
    if ( info.imprim )  fprintf(stdout,"     - count = %d, countmax = %d, count2 = %d, count3 = %d, \n",count,countmax,count2,count3);
    if (info.debug) {
    it++;
      saveMesh(&extmesh,it);
    }
  }
  chrono(OFF,&info.ctim[2]);
  printim(info.ctim[2].gdif,stim);
  
  if ( info.imprim ) fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);
  
  
  /* save file */
  if ( info.imprim ) fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",extmesh.name);
  if(!info.debug) if (!saveMesh(&extmesh,0)) return (1);
  chrono(ON,&info.ctim[1]);
  free(extmesh.tetra);
  free(extmesh.point);
  free(intmesh.tetra);
  free(intmesh.point);
  free(bucket->head);
  free(bucket->link);
  free(bucket);
  
  /* free mem */
  if ( sol.u )  free(sol.u);
  chrono(OFF,&info.ctim[0]);
  fprintf(stdout,"\n  %s\n   END OF MODULE WARPING. %s \n  %s\n",LS_STR,stim,LS_STR);
  return(0);
}
