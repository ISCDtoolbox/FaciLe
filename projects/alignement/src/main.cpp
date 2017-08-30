#include "mesh.h"
#include "align.h"

char* source = "source.mesh";
char* target = "target.mesh";

//Parameters for Super4PCS
float overlap = 0.2;
float delta = 5.0;
int n_points = 1000;

//Parameters for ICP
float inlierDist = 1.0;
int maxIt = 100;

#include <cstring>

float* transformVert(const float* mat, float* vert){
  float *V = new float[3];
  V[0]=V[1]=V[2]=0;
  V[0] = mat[0]*vert[0] + mat[1]*vert[1] + mat[2]*vert[2] + mat[3];
  V[1] = mat[4]*vert[0] + mat[5]*vert[1] + mat[6]*vert[2] + mat[7];
  V[2] = mat[8]*vert[0] + mat[9]*vert[1] + mat[10]*vert[2] + mat[11];
  return V;
}

void getArgs(int argc, char **argv) {
  int i = 1;
  while (i < argc) {
    //Input
    if (!strcmp(argv[i], "-i")) {
      source = argv[++i];
      target = argv[++i];
    }
    //Super4PCS
    else if (!strcmp(argv[i], "-o")) {
      overlap = atof(argv[++i]);
    }
    else if (!strcmp(argv[i], "-d")) {
      delta = atof(argv[++i]);
    }
    else if (!strcmp(argv[i], "-n")) {
      n_points = atoi(argv[++i]);
    }
    //ICP
    else if (!strcmp(argv[i], "-m")) {
      maxIt = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-c")) {
      inlierDist = atof(argv[++i]);
    }
    //Usage
    else if (!strcmp(argv[i], "-h")) {
      std::cout << "Usage: " << argv[0] << " -i source.mesh target.mesh" << std::endl;
      std::cout << "\t[ -o overlap (" << overlap << ") ]" << std::endl;
      std::cout << "\t[ -d delta (" << delta << ") ]" << std::endl;
      std::cout << "\t[ -n n_points (" << n_points << ") ]" << std::endl;
      std::cout << "\t[ -c inlier (" << inlierDist << ") ]" << std::endl;
      std::cout << "\t[ -m max iterations (" << maxIt << ") ]" << std::endl;
      exit(0);
    }
    else if (argv[i][0] == '-') {
      std::cout << "Unknown flag" << std::endl;
      exit(-1);
    };
    i++;
  }
}

int main(int argc, char ** argv ) {

  getArgs(argc, argv);
  Mesh *sourceMesh = new Mesh(source);
  Mesh *targetMesh = new Mesh(target);

  std::cout << "vert = " << sourceMesh->vertices[0] << " " << sourceMesh->vertices[1] << " " << sourceMesh->vertices[2] << std::endl;

  //Compute the super4PCS registration
  const float *mat = super4PCS(*targetMesh, *sourceMesh, overlap, delta, n_points);

  //Transform the target object
  for (int i = 0 ; i < sourceMesh->vertices.size()/3 ; i++){
    float* newVert = transformVert(mat, &(sourceMesh->vertices[3*i]));
    sourceMesh->vertices[3*i+0] = newVert[0];
    sourceMesh->vertices[3*i+1] = newVert[1];
    sourceMesh->vertices[3*i+2] = newVert[2];
  }
  sourceMesh->write("output.mesh");

  std::cout << "vert = " << sourceMesh->vertices[0] << " " << sourceMesh->vertices[1] << " " << sourceMesh->vertices[2] << std::endl;


  const float *matICP = icp(*targetMesh, *sourceMesh, maxIt, inlierDist);
  //Transform the target object
  for (int i = 0 ; i < sourceMesh->vertices.size()/3 ; i++){
    float* newVert = transformVert(matICP, &(sourceMesh->vertices[3*i]));
    sourceMesh->vertices[3*i+0] = newVert[0];
    sourceMesh->vertices[3*i+1] = newVert[1];
    sourceMesh->vertices[3*i+2] = newVert[2];
  }

  sourceMesh->write("output2.mesh");

  return 0;
}