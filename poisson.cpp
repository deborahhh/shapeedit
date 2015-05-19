#include "poisson.h"
#include <math.h>
#include <iostream>
#include "Eigen/SparseCore"
#include "Eigen/SparseCholesky"
#include "Eigen/SVD"

// access to the original (undeformed) mesh
//
extern Mesh g_originalMesh;

using namespace Eigen;

// you may have to tweak this to get smooth motion
//
int g_iterationsPerSecond = 8;


/*
 * GLOBAL VARIABLES
 */

Mesh g_mesh;

MatrixXd g_weightMatrix;
vector<vector<double> > g_weightVector;

vector<vector<Cvec2> > g_initialPositions;
vector<vector<Cvec2> > g_targetPositions;
vector<MatrixXd> g_rotations;

int iterations = 0;


/*
 * HELPER FUNCTIONS
 */

// calculated the cotangent of the angle <abc
double getCotangent(Cvec2 a, Cvec2 b, Cvec2 c) {
  Cvec2 adj = (c-b) * (dot(c-b,a-b) / norm2(c-b));
  Cvec2 opp = (a-b) - adj;
  return norm(adj)/norm(opp);
}

MatrixXd getWeights(Mesh& m) {

  int n = m.getNumVertices();
  MatrixXd weights = MatrixXd::Zero(n,n);

  int tf = m.getNumFaces();
  for (int i = 0; i < tf; i++) {

    Mesh::Face f = m.getFace(i);
    assert(f.getNumVertices() == 3);

    // get vertices 
    Mesh::Vertex a = f.getVertex(0), b = f.getVertex(1), c = f.getVertex(2);

    // get positions
    Cvec2 pa = a.getPosition(), pb = b.getPosition(), pc = c.getPosition();

    // get vertices indexes
    int ia = a.getIndex(), ib = b.getIndex(), ic = c.getIndex();

    // get cotagents
    double cota = getCotangent(pc,pa,pb), cotb = getCotangent(pa,pb,pc), cotc = getCotangent(pb,pc,pa);

    // add cotangents to matrix
    weights(ia,ia) += cotb + cotc;
    weights(ib,ib) += cota + cotc;
    weights(ic,ic) += cota + cotb;
    weights(ia,ib) -= cotc;
    weights(ib,ia) -= cotc;
    weights(ia,ic) -= cotb;
    weights(ic,ia) -= cotb;
    weights(ib,ic) -= cota;
    weights(ic,ib) -= cota;

  }

  return weights;
}


// gets a vector of positions from a mesh
vector<vector<Cvec2> > getPositionsFromMesh(Mesh& m) {

  int n = m.getNumVertices();
  vector<vector<Cvec2> > positions(n);

  int e = m.getNumEdges();

  for (int idx = 0; idx < e; idx++) {
    const Mesh::Edge ed = m.getEdge(idx);

    Mesh::Vertex vi = ed.getVertex(0), vj = ed.getVertex(1);
    int i = vi.getIndex(), j = vj.getIndex();
    double wij = -g_weightMatrix(i,j);

    positions.at(i).push_back(vi.getPosition()-vj.getPosition());
    positions.at(j).push_back(vj.getPosition()-vi.getPosition());
  }

  return positions;
}






// This function is called at the beginning of the program
// before any handles have been added.
//
void initPoisson(Mesh& m) {

  g_weightMatrix = getWeights(m);
  int e = m.getNumEdges();
  int n = m.getNumVertices();

  vector<vector<double> > weightVector(n);
  
  for (int idx = 0; idx < e; idx++) {
    const Mesh::Edge ed = m.getEdge(idx);

    Mesh::Vertex vi = ed.getVertex(0), vj = ed.getVertex(1);
    int i = vi.getIndex(), j = vj.getIndex();
    double wij = -g_weightMatrix(i,j);

    weightVector.at(i).push_back(wij);
    weightVector.at(j).push_back(wij);
  }

  g_weightVector = weightVector;

  g_initialPositions = getPositionsFromMesh(m);

  g_mesh = m;

}


// This function is called immediately after each movement
// with an updated set of handles and the most recent mesh.
//
void afterMove(Mesh& m, vector<handleType>& handles) {


  g_targetPositions = getPositionsFromMesh(m);

  iterations = 0;

}


/*
 * Helper functions for iterations
 */


// calculates the rotation matrices
void getRs() {

  int n = g_initialPositions.size();
  vector<MatrixXd> rs(n);

  for (int idx = 0; idx < n; idx++) {
    vector<Cvec2> vecP = g_initialPositions.at(idx);
    vector<Cvec2> vecP2 = g_targetPositions.at(idx);
    vector<double> vecW = g_weightVector.at(idx);
    int nbs = vecP.size();
    assert(nbs == vecW.size());
    assert(nbs == vecP2.size());

    MatrixXd matP(2,nbs);
    MatrixXd matP2(2,nbs);
    MatrixXd matW = MatrixXd::Zero(nbs,nbs);

    for (int j = 0; j < nbs; j++) {
      matP.col(j) << vecP.at(j)[0], vecP.at(j)[1];
      matP2.col(j) << vecP2.at(j)[0], vecP2.at(j)[1];
      matW(j,j) = vecW.at(j);
    }

    MatrixXd matS = matP * matW * matP2.transpose();

    JacobiSVD<MatrixXd> svd(matS, ComputeThinU | ComputeThinV);

    MatrixXd matR = svd.matrixV() * svd.matrixU().transpose();
    rs.at(idx) = matR;

  }

  g_rotations = rs;
}


// calculates the matrix B for the system
MatrixXd getMatrixB(MatrixXd& weights, vector<MatrixXd> rs) {

  Mesh m = g_mesh;

  int n = m.getNumVertices();
  MatrixXd matB = MatrixXd::Zero(n,2);

  int e = m.getNumEdges();
  for (int idx = 0; idx < e; idx++) {

    // get edge and its vertices
    Mesh::Edge ed = m.getEdge(idx);
    Mesh::Vertex vi = ed.getVertex(0), vj = ed.getVertex(1);
    int i = vi.getIndex(), j = vj.getIndex();

    // look up edge weight
    double wij = -weights(i,j);

    // get edge vector
    Vector2d eij;
    eij << (vi.getPosition() - vj.getPosition())[0], (vi.getPosition() - vj.getPosition())[1];

    // compute added vector and add to matB
    MatrixXd added = (wij/2) * (rs.at(i) + rs.at(j)) * eij;
    matB.row(i) << matB.row(i) + added.transpose();
    matB.row(j) << matB.row(j) - added.transpose();
  }

  return matB;
}

// compare function to order handles by vertex index
bool compare_handle(const handleType& a, const handleType& b) {
  return a.first < b.first;
}

// removes rows and cols of the system to adapt to fixed handles
pair<MatrixXd,MatrixXd> getSystem(MatrixXd matA, MatrixXd matB, vector<handleType>& handles) {

  /*
   * adapt system to take handles into account
   */

  // remove columns and pass to B side
  int n = matA.cols();
  assert(n == matA.rows());
  std::sort (handles.begin(), handles.end(), compare_handle);
  int hsz = handles.size();

  MatrixXd newA0(n,n-hsz);

  int counter = 0;
  for (int i = 0; i < n; i++) {
    // if current vertex is a fixed one, add to the B side
    if (counter < hsz && handles.at(counter).first == i) {
      Cvec2 pos = handles.at(counter).second;
      MatrixXd pt(1,2);
      pt << pos[0], pos[1];

      MatrixXd thiscol(n,1);
      thiscol = matA.col(i);
      matB -= thiscol * pt;
      counter++;
    }
    // else, just copy the column to the new matrix
    else {
      newA0.col(i-counter) = matA.col(i);
    }
  }


  // remove rows
  MatrixXd newA(n-hsz,n-hsz);
  MatrixXd newB(n-hsz,2);
  counter = 0;
  for (int i = 0; i < n; i++) {
    if (counter < hsz && handles.at(counter).first == i) {
      counter++;
    }
    else {
      for (int j = 0; j < n-hsz; j++) {
        newA(i-counter,j) = newA0(i,j);
      }
      newB(i-counter,0) = matB(i,0);
      newB(i-counter,1) = matB(i,1);
    }
  }

  return std::make_pair(newA,newB);
}

// solve system, given the matrices
VectorXd solveSystem(MatrixXd matA, MatrixXd matB) {

  int sz = matA.rows();
  assert(sz == matA.cols());

  // solve system
  MatrixXd matA2 = MatrixXd::Zero(2*sz,2*sz);
  for (int i = 0; i < sz; i++) {
    for (int j = 0; j < sz; j++) {
      matA2(i,j) = matA(i,j);
      matA2(i+sz,j+sz) = matA(i,j);
    }
  }

  VectorXd vecB2(2*sz);
  vecB2.block(0,0,sz,1) = matB.block(0,0,sz,1);
  vecB2.block(sz,0,sz,1) = matB.block(0,1,sz,1);


  SparseMatrix<double> newA = matA2.sparseView();
  SimplicialLLT<SparseMatrix<double> > solver;
  solver.compute(newA);
  VectorXd x = solver.solve(vecB2);

  return x;
}




/*
 * END OF STUFF
 */

// This function does a single iteration of the algorithm.
// It must return true if convergence has been achieved,
// and return false if it hasn't been achieved.
// It is called after each movement, and depending on whether
// convergence was reached, is scheduled to call again
// after 1/iterationsPerSecond seconds. If it is already scheduled
// to call when a movement is done, then it is not called,
// although its next call will come with updated handles.
// However in this case, afterMove does get called.
//
bool doIteration(Mesh& m, vector<handleType>& handles) {

  // get rotation matrices
  getRs();

  // get B matrix of the system
  MatrixXd matB = getMatrixB(g_weightMatrix, g_rotations);

  // remove rows to adapt to handles
  pair<MatrixXd,MatrixXd> system = getSystem(g_weightMatrix, matB, handles);

  // solve system
  VectorXd x = solveSystem(system.first,system.second);

  int n = m.getNumVertices();
  int hsz = handles.size();

  // get vertex positions
  vector<Cvec2> positions(n);
  int counter = 0;
  for (int i = 0; i < n; i++) {
    if (counter < hsz && handles.at(counter).first == i) {
      positions.at(i) = handles.at(counter).second;
      counter++;
    }
    else {
      positions.at(i) = Cvec2(x(i-counter),x(n-hsz+i-counter));
    }
  }



  // update vertex positions
  double error = 0;
  for (int i = 0; i < n; i++) {
    error += norm2(m.getVertex(i).getPosition() - positions.at(i));
    m.getVertex(i).setPosition(positions.at(i));
  }


  g_targetPositions = getPositionsFromMesh(m);
  error /= n; 

  // return true;
  cout << "error " << error << endl;

  return (error < 10e-4);
}

