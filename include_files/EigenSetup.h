//#define EIGEN_USE_MKL_ALL
//#define EIGEN_DEFAULT_DENSE_INDEX_TYPE MKL_INT
#include <Eigen/Dense>
#include <Eigen/SparseCore>

// setup relevant vector/matrix datatypes from Eigen

using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXi;
using Eigen::SparseMatrix;
typedef Eigen::SparseMatrix<double,Eigen::RowMajor> Sparse_Matrix;
typedef Eigen::Triplet<double> triplet;
typedef Eigen::MatrixXd Full_matrix;

// matrix type for the unsigned int
typedef Matrix<unsigned int,Dynamic,Dynamic> MatrixUI;

