#ifndef MESHDEF
#define MESHDEF

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

class Mesh{
public:
  int nbVertices;

  std::vector<float>  vertices;
  std::vector<float>  normals;
  std::vector<float>  colors;
  std::vector<int>    refTri, refTet, refVert;
  std::vector<int>    indTri, indTet;
  std::vector<int>    adjacent;

  std::string         meshfile;
  std::string         name;

  Mesh(char * mesh_path);
  int write(char * out_path);

};

#endif
