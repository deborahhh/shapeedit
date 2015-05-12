#include "poisson.h"
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

// This function is called at the beginning of the program
// before any handles have been added.
//
void initPoisson(Mesh& m) {
// TODO
//
}

// This function is called immediately after each movement
// with an updated set of handles and the most recent mesh.
//
void afterMove(Mesh& m, vector<handleType>& handles) {
// TODO
//
}


/*
 *
 * STUFF 
 *
 */



/* helper functions */

bool compare_handle(const handleType& a, const handleType& b) {
  return a.first < b.first;
}


double getCotangent(Mesh::Vertex v, Mesh::Vertex w) {
  // TODO
  return 1.0;
}


vector<Matrix2d> getRs(Mesh& m) {

  int n = m.getNumVertices();
  vector<Matrix2d> rs(n);

  int e = m.getNumEdges();
  vector<vector<Cvec2> > positions(n);
  vector<vector<double> > weights(n);

  for (int idx = 0; idx < e; idx++) {
    const Mesh::Edge ed = m.getEdge(idx);

    Mesh::Vertex vi = ed.getVertex(0), vj = ed.getVertex(1);
    int i = vi.getIndex(), j = vj.getIndex();
    double wij = getCotangent(vi,vj);

    positions.at(i).push_back(vj.getPosition());
    positions.at(j).push_back(vi.getPosition());
    weights.at(i).push_back(wij);
    weights.at(j).push_back(wij);
  }

  for (int idx = 0; idx < n; idx++) {
    vector<Cvec2> vecP = positions.at(idx);
    vector<double> vecW = weights.at(idx);
    int nbs = vecP.size();
    assert(nbs == vecW.size());

    MatrixXd matP(2,nbs);
    MatrixXd matW = MatrixXd::Zero(nbs,nbs);

    for (int j = 0; j < nbs; j++) {
      matP.col(j) << vecP.at(j)[0], vecP.at(j)[1];
      matW(j,j) = vecW.at(j);
    }


    MatrixXd matS = matP * matW * matP.transpose();

    JacobiSVD<MatrixXd> svd(matS, ComputeThinU | ComputeThinV);

    Matrix2d matR = svd.matrixV() * svd.matrixU().transpose();
    rs.at(idx) = matR;
  }




/*

  int n = m.getNumVertices();
  vector<Matrix2d> rs(n);



  for (int i = 0; i < n; ++i) {
    // neighbor vertex
    const Mesh::Vertex v = m.getVertex(i);

    Mesh::VertexIterator it(v.getIterator()), it0(it);

    vector<Cvec2> ps;
    vector<double> ws;
    int nbs = 0;

    do {
      Mesh::Vertex w = it.getVertex();
      ps.push_back(w.getPosition());
      ws.push_back(getCotangent(v,w));
      nbs++;
    }
    while (++it != it0);      // go around once the 1ring

    MatrixXd matP(2,nbs);
    MatrixXd matW = MatrixXd::Zero(nbs,nbs);
    for (int j = 0; j < nbs; j++) {
      matP.col(j) << ps.at(j)[0], ps.at(j)[1];
      matW(j,j) = ws.at(j);
    }

    MatrixXd matS = matP * matW * matP.transpose();

    JacobiSVD<MatrixXd> svd(matS, ComputeThinU | ComputeThinV);

    Matrix2d matR = svd.matrixV() * svd.matrixU().transpose();
    rs.at(i) = matR;
  }*/

  return rs;
}

pair<SparseMatrix<double>,MatrixXd> getSystem(Mesh& m, vector<Matrix2d> rs) {

  int n = m.getNumVertices();
  SparseMatrix<double> weights(n,n);
  MatrixXd bs = MatrixXd::Zero(n,2);

  int e = m.getNumEdges();
  for (int idx = 0; idx < e; idx++) {
    const Mesh::Edge ed = m.getEdge(idx);

    Mesh::Vertex vi = ed.getVertex(0), vj = ed.getVertex(1);
    int i = vi.getIndex(), j = vj.getIndex();
    double wij = getCotangent(vi,vj);

    weights.insert(i,j) = -wij;
    weights.insert(j,i) = -wij;
    weights.coeffRef(i,i) += wij;
    weights.coeffRef(j,j) += wij;


    Vector2d eij;
    eij << (vi.getPosition() - vj.getPosition())[0], (vi.getPosition() - vj.getPosition())[1];
    cout << "rsi\n" << rs.at(i) << endl;
    cout << "rsj\n" << rs.at(j) << endl;
    cout << "eij" << eij << endl;
    MatrixXd added = (wij/2) * (rs.at(i) + rs.at(j)) * eij;
    cout << "added\n" << added << endl;
    int test1 = added.cols(), test2 = added.rows();
    cout << "bs\n" << bs << endl;
    bs.row(i) << bs.row(i) + added.transpose();
    cout << "bs1\n" << bs << endl;
    bs.row(j) << bs.row(j) - added.transpose();
    cout << "bs2\n" << bs << endl;
  }


  pair<SparseMatrix<double>,MatrixXd> rtn(weights,bs);
  return rtn;

/*
  int n = m.getNumVertices();
  SparseMatrix<double> weights(n,n);
  MatrixXd bs = MatrixXd::Zero(n,2);

  for (int i = 0; i < n; ++i) {
    const Mesh::Vertex v = m.getVertex(i);

    Mesh::VertexIterator it(v.getIterator()), it0(it);
    do
    {
      Mesh::Vertex w = it.getVertex();
      int j = w.getIndex();
      double wij = getCotangent(v,w);
      weights.insert(i,j) = -wij;
      weights.coeffRef(i,i) += wij;

      Vector2d eij;
      eij << (v.getPosition() - w.getPosition())[0],  (v.getPosition() - w.getPosition())[1];
      MatrixXd added = (wij/2) * (rs.at(i) + rs.at(j)) * eij;
      bs.row(i) << added;
      // can use here it.getVertex(), it.getFace()
    }
    while (++it != it0);      // go around once the 1ring
  }
  pair<SparseMatrix<double>,MatrixXd> rtn(weights,bs);
  return rtn;

  pair<SparseMatrix<double>,MatrixXd> rtn(weights,bs);
  return rtn;*/
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
// TODO
// Currently, this just moves every vertex directly.
//
  // for (int i = 0; i < m.getNumFaces(); ++i) {
  //   cout << "face" << i << endl;
  //   const Mesh::Face f = m.getFace(i);

  //   int total = f.getNumVertices();
  //   cout << "# vertices " << total << endl;
  //   for (int j = 0; j < total; j++) {
  //     cout << "vertex: " << f.getVertex(j).getIndex() << endl;
  //   }

  //                                   // go around once the 1ring
  // }
  // get Rs based on current mesh
  vector<Matrix2d> rs = getRs(m);


  // get system to update mesh
  pair<SparseMatrix<double>,MatrixXd> system = getSystem(m, rs);


  // adapt system to take handles into account
  int n = m.getNumVertices();
  std::sort (handles.begin(), handles.end(), compare_handle);
  int hsz = handles.size();

  SparseMatrix<double> newA0(n,n-hsz);
  SparseMatrix<double> matA = system.first;
  MatrixXd matB0 = system.second;

  int counter = 0;

  printf("lol");
  for (int i = 0; i < n; i++) {
    if (counter < hsz && handles.at(counter).first == i) {
      Cvec2 pos = handles.at(counter).second;
      MatrixXd pt(1,2);
      pt << pos[0], pos[1];

      MatrixXd thiscol(n,1);
      thiscol = matA.col(i);
      matB0 -= thiscol * pt;
      counter++;
    }
    else {
      newA0.col(i-counter) = matA.col(i);
    }
  }

  SparseMatrix<double> newA(n-hsz,n-hsz);
  MatrixXd matB(n-hsz,2);
  // remove rows
  counter = 0;
  for (int i = 0; i < n; i++) {
    if (counter < hsz && handles.at(counter).first == i) {
      counter++;
    }
    else {
      for (int j = 0; j < n-hsz; j++) {
        newA.coeffRef(i-counter,j) = newA0.coeffRef(i,j);
      }
      matB(i-counter,0) = matB0(i,0);
      matB(i-counter,1) = matB0(i,1);
    }
  }

  // for (int i = 0; i < hsz; i++) {
  //   int row = handles.at(hsz-i-1).first; // start from the last one
  //   matB.block(row,0,n-row,2) = matB.block(row+1,0,n-row,2);
  // }
  // matB.conservativeResize(n-hsz,2);

  printf("lol");
  // solve system
  SparseMatrix<double> matA2(2*(n-hsz),2*(n-hsz));
  // matA2.block(0,0,n,n-hsz) = newA.block(0,0,n,n-hsz);
  // matA2.block(n,n-hsz,n,n-hsz) = newA.block(0,0,n,n-hsz);
  for (int i = 0; i < n-hsz; i++) {
    for (int j = 0; j < n-hsz; j++) {
      matA2.insert(i,j) = newA.coeffRef(i,j);
      matA2.insert(i+n-hsz,j+n-hsz) = newA.coeffRef(i,j);
    }
  }

  printf("lol");
  VectorXd vecB2(2*(n-hsz));
  vecB2.block(0,0,n-hsz,1) = matB.block(0,0,n-hsz,1);
  vecB2.block(n-hsz,0,n-hsz,1) = matB.block(0,1,n-hsz,1);


  SimplicialLLT<SparseMatrix<double> > solver;
  solver.compute(matA2);
  VectorXd x = solver.solve(vecB2);


  printf("lol");
  // get vertex positions
  vector<Cvec2> positions(n);
  counter = 0;
  for (int i = 0; i < n; i++) {
    if (counter < hsz && handles.at(counter).first == i) {
      positions.at(i) = handles.at(counter).second;
      counter++;
    }
    else {
      positions.at(i) = Cvec2(x(i-counter),x(n-hsz+i-counter));
    }
  }
  printf("lol");


  // update vertex positions
  double error = 0;
  for (int i = 0; i < n; i++) {
    error += norm2(m.getVertex(i).getPosition() - positions.at(i));
    m.getVertex(i).setPosition(positions.at(i));
  }

  printf("lol");
 	// for (vector<handleType>::iterator i = handles.begin(); i != handles.end(); ++i)
 	// 	m.getVertex(i->first).setPosition(i->second);

  // return error < SOME_CONSTANT
 	return true;
}

