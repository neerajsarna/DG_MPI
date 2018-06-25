	#include <Eigen/Dense>
	#include <Eigen/SparseCore>
	typedef Eigen::SparseMatrix<double,Eigen::RowMajor> Sparse_Matrix;


	struct system_data
	{		
		// flux matrices
		Sparse_Matrix Ax;
		Sparse_Matrix Ay;
		Sparse_Matrix P;
		
		// boundary matrix
		std::vector<Sparse_Matrix> B; 
		std::vector<Sparse_Matrix> penalty;
		std::vector<Sparse_Matrix> penalty_B;

	};
