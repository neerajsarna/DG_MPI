#include "EigenSetup.h"


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

		// by default we have a false
		bool bc_inhomo_time = false;
		bool have_force = false;
	};
