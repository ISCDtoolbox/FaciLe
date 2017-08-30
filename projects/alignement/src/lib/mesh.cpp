extern "C"
{
  #include "libmesh5.h"
}

#include "mesh.h"

Mesh::Mesh(char * mesh_path){
  std::string meshFile    = std::string(mesh_path);

  //Initialisation
  int nPts, nTri, nNor, nTet, nNorAtV;
  int ver, dim;
  double* tmp = new double[3];

  //READING .mesh
  int inm = GmfOpenMesh(mesh_path,GmfRead,&ver,&dim);
  if ( !inm ){
    std::cout << "Unable to open mesh file " << mesh_path << std::endl;
    exit(-1);
  }

  //GETTING SIZES
  nPts    = GmfStatKwd(inm, GmfVertices);
  nTri    = GmfStatKwd(inm, GmfTriangles);
  nTet    = GmfStatKwd(inm, GmfTetrahedra);
  nNor    = GmfStatKwd(inm, GmfNormals);
  nNorAtV = GmfStatKwd(inm, GmfNormalAtVertices);
  if ( !nPts || !nTri ){
    std::cout << "Missing data in mesh file" << mesh_path << std::endl;
    exit(-1);
  }
  vertices.resize(3 * nPts);
  refVert.resize(nPts);
  indTri.resize(3 * nTri);
  refTri.resize(nTri);
  indTet.resize(4 * nTet);
  refTet.resize(nTet);


  //VERTICES & INDICES
  GmfGotoKwd(inm,GmfVertices);
  for (int k = 0; k < nPts; k++){
    GmfGetLin(inm,GmfVertices,&tmp[0],&tmp[1],&tmp[2], &refVert[k]);
    vertices[3*k + 0] = tmp[0];
    vertices[3*k + 1] = tmp[1];
    vertices[3*k + 2] = tmp[2];
  }
  GmfGotoKwd(inm,GmfTriangles);
  for (int k = 0; k < nTri; k++){
    GmfGetLin(inm,GmfTriangles,&indTri[3*k],&indTri[3*k+1], &indTri[3*k+2], &refTri[k]);
    indTri[3*k]-=1;
    indTri[3*k+1]-=1;
    indTri[3*k+2]-=1;
  }
  nbVertices = vertices.size()/3;

  //TETRAHEDRON
  if(nTet){
    GmfGotoKwd(inm,GmfTetrahedra);
    for (int k = 0 ; k < nTet ; k++)
    GmfGetLin(inm,GmfTetrahedra,&indTet[4*k], &indTet[4*k+1], &indTet[4*k+2], &indTet[4*k+3],&refTet[k]);
  }

  //NORMALS
  std::vector<float> tmp_normals;
  std::vector<int> NormalAtVertices;
  if(nNor)
  tmp_normals.resize(3 * nNor);
  normals.resize(3 * nNor);
  if(nNorAtV)
  NormalAtVertices.resize(2 * nNorAtV + 2);

  if(nNor && nNorAtV){
    GmfGotoKwd(inm,GmfNormals);
    for (int k = 0 ; k < nNor ; k++) {
      GmfGetLin(inm,GmfNormals,&tmp[0],&tmp[1],&tmp[2]);
      double dd = 1.0 / std::sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
      tmp_normals[3*k + 0] = tmp[0] * dd;
      tmp_normals[3*k + 1] = tmp[1] * dd;
      tmp_normals[3*k + 2] = tmp[2] * dd;
    }
    GmfGotoKwd(inm,GmfNormalAtVertices);
    for (int k = 0 ; k < nNorAtV; k++)
    GmfGetLin(inm,GmfNormalAtVertices, &NormalAtVertices[2*k+0], &NormalAtVertices[2*k+1]);
    for(int i = 0 ; i < nNorAtV - 1 ; i++){
      int indV = NormalAtVertices[2*i+0] - 1;
      int indN = NormalAtVertices[2*i+1] - 1;
      normals[3 * indV + 0] =  tmp_normals[3*indN+0];
      normals[3 * indV + 1] =  tmp_normals[3*indN+1];
      normals[3 * indV + 2] =  tmp_normals[3*indN+2];
    }
  }

  std::cout << "Succesfully opened  " << mesh_path << std::endl;
  GmfCloseMesh(inm);


  //Lecture du .sol
  std::string solFile = meshFile.substr(0, meshFile.size()-5) + ".sol";
  int ver2, dim2;
  int inMeshSol = 0;
  inMeshSol = GmfOpenMesh((char*)(solFile.c_str()), GmfRead, &ver2, &dim2);
  if(inMeshSol){
    int type, offset, typtab[GmfMaxTyp];
    int nSol      = GmfStatKwd(inMeshSol, GmfSolAtVertices, &type, &offset, &typtab);
    std::vector<float> values(nSol);
    GmfGotoKwd(inMeshSol, GmfSolAtVertices);
    for(int i = 0 ; i< nSol ; i++){
      double val;
      if ( ver2 == GmfFloat )
      GmfGetLin(inMeshSol, GmfSolAtVertices, &values[i]);
      else{
        GmfGetLin(inMeshSol, GmfSolAtVertices, &val);
        values[i] = val;
      }
    }
    std::cout << "Succesfully opened  " << solFile << std::endl;
    GmfCloseMesh(inMeshSol);

    //From palette
    float mini = 1e9;
    float maxi = -1e9;
    for(int i = 0 ; i < values.size() ; i++){
      mini = std::min(mini, values[i]);
      maxi = std::max(maxi, values[i]);
    }
    /*
    std::vector<glm::vec4> P;
    P.push_back(glm::vec4(0.0f,   1.0f,   0.0f,   0.0f));
    P.push_back(glm::vec4(0.35, 0.75,   0.75, 0.0f));
    P.push_back(glm::vec4(0.5, 0.0f, 1.0f, 0.0f));
    P.push_back(glm::vec4(0.65, 0.0f, 0.75, 0.75));
    P.push_back(glm::vec4(1.0f,   0.0f, 0.0f, 1.0f));

    for(int k = 0 ; k < values.size() ; k++){
      float val = (maxi - values[k]) / (maxi - mini);
      for(int i = 0 ; i < P.size() ; i++){
        if(i==P.size()-1){
          float fac       = (val - P[i-1][0])/(P[i][0] - P[i-1][0]);
          glm::vec4 col   = (1-fac)*P[i-1] + fac*P[i];
          colors.push_back(col[1]);
          colors.push_back(col[2]);
          colors.push_back(col[3]);
        }
        else{
          if( (val>=P[i][0]) && (val<=P[i+1][0]) ){
            float fac       = (val - P[i][0])/(P[i+1][0] - P[i][0]);
            glm::vec4 col   = (1-fac)*P[i] + fac*P[i+1];
            colors.push_back(col[1]);
            colors.push_back(col[2]);
            colors.push_back(col[3]);
            break;
          }
        }
      }
    }
    */
  }
  else{
    std::cout << "Failed to open .sol " << solFile << std::endl;
  }
}
int Mesh::write(char* out_path){
  std::ofstream output(out_path);
  if (output.is_open())
  {
    output << "MeshVersionFormatted 2\nDimension 3\n";

    output << "Vertices\n";
    output << int(vertices.size()/3) << "\n";
    for(int i = 0; i < vertices.size()/3; i ++){
        output << vertices[3*i+0] << " " << vertices[3*i+1] << " " << vertices[3*i+2] << " 0\n";
    }

    output << "Triangles\n";
    output << int(indTri.size()/3) << "\n";
    for(int i = 0; i < indTri.size()/3; i ++){
        output << indTri[3*i+0]+1 << " " << indTri[3*i+1]+1 << " " << indTri[3*i+2]+1 << " 0\n";
    }

    output.close();
  }
  else std::cout << "Unable to open file";
  return 0;
}
