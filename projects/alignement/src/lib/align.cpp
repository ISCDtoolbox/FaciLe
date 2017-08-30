#include "align.h"

#include <vector>

using namespace GlobalRegistration;

struct TransformVisitor {
    inline void operator() (
            float fraction,
            float best_LCP,
            Eigen::Ref<Match4PCSBase::MatrixType> /*transformation*/) {
        printf("done: %d%c best: %f                  \r",
               static_cast<int>(fraction * 100), '%', best_LCP);
        fflush(stdout);
    }
    constexpr bool needsGlobalTransformation() const { return false; }
};

const float* super4PCS(Mesh &sourceMesh, Mesh &targetMesh, float overlap, float delta, int n_points){
  /*
  Note: We have to invert source and target,
  as the resulting matrix is going from target to source
  */
  std::vector<Point3D> set1, set2;
  std::vector<tripple> tris1, tris2;

  for(int i = 0 ; i < sourceMesh.vertices.size()/3; i++){
    Point3D pt(
      sourceMesh.vertices[3*i+0],
      sourceMesh.vertices[3*i+1],
      sourceMesh.vertices[3*i+2]
    );
    set1.push_back(pt);
  }
  for(int i = 0 ; i < targetMesh.vertices.size()/3; i++){
    set2.push_back(
      Point3D(
        targetMesh.vertices[3*i+0],
        targetMesh.vertices[3*i+1],
        targetMesh.vertices[3*i+2]
      )
    );
  }

  Match4PCSOptions options;
  Match4PCSBase::MatrixType mat;
  bool overlapOk = options.configureOverlap(overlap);
  if(! overlapOk )  {
      std::cerr << "Invalid overlap configuration. ABORT" << std::endl;
      /// TODO Add proper error codes
      exit(-3);
  }
  options.sample_size = n_points;
  options.delta = delta;
  // Match and return the score (estimated overlap or the LCP).
  typename Point3D::Scalar score = 0;

  constexpr Utils::LogLevel loglvl = Utils::Verbose;
  using TrVisitorType = typename std::conditional <loglvl==Utils::NoLog,Match4PCSBase::DummyTransformVisitor,TransformVisitor>::type;
  Utils::Logger logger(loglvl);

  MatchSuper4PCS matcher(options, logger);
  logger.Log<Utils::Verbose>( "Use Super4PCS" );
  score = matcher.ComputeTransformation<TrVisitorType>(set1, &set2, mat);

  std::cout << mat.cast<double>() << std::endl;
  mat.transposeInPlace();

  return mat.data();
}

const float* icp(Mesh &sourceMesh, Mesh &targetMesh, int maxIt, float inlierDist){
  // define a 3 dim problem with 10000 model points
  // and 10000 template points:
  int sourceSampling = 10000;
  int targetSampling = 10000;

  // allocate model and template memory
  double* M = (double*)calloc(sourceMesh.vertices.size(),sizeof(double));
  double* T = (double*)calloc(targetMesh.vertices.size(),sizeof(double));

  //Create models
  for(int i = 0 ; i < sourceMesh.vertices.size() ; i++)
    M[i] = (double)(sourceMesh.vertices[i]);
  for(int i = 0 ; i < targetMesh.vertices.size() ; i++)
    T[i] = (double)(targetMesh.vertices[i]);

  // start with identity as initial transformation
  // in practice you might want to use some kind of prediction here
  Matrix R = Matrix::eye(3);
  Matrix t(3,1);

  // run point-to-plane ICP (-1 = no outlier threshold)
  std::cout << std::endl << "Running ICP (point-to-point, no outliers)" << std::endl;
  IcpPointToPoint icp(M,sourceSampling,3);
  icp.setMaxIterations(maxIt);
  double residual = icp.fit(T,targetSampling,R,t,inlierDist);

  // results
  std::cout << std::endl << "Transformation results:" << std::endl;
  std::cout << "R:" << std::endl << R << std::endl << std::endl;
  std::cout << "t:" << std::endl << t << std::endl << std::endl;
  std::cout << "Residual:"<<residual << std::endl;

  float* mat;
  mat = new float[16];
  for(int i = 0 ; i < 3 ; i++){
    for(int j = 0 ; j < 3 ; j++){
      mat[4*i+j] = R.val[i][j];
    }
  }
  mat[12]=mat[13]=mat[14]=0;
  mat[3]=t.val[0][0];
  mat[7]=t.val[1][0];
  mat[11]=t.val[2][0];
  mat[15]=1;

  // free memory
  free(M);
  free(T);
  return mat;
}
