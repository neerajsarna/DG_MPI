#include "EigenSetup.h"


	struct system_data
	{		
		// flux matrices
		Sparse_Matrix Ax;
		Sparse_Matrix Ay;
		Sparse_Matrix P;

		Sparse_Matrix Ax_mod;
		Sparse_Matrix Ay_mod;
		
		// boundary matrix
		std::vector<Sparse_Matrix> B; 
		std::vector<Sparse_Matrix> penalty;
		std::vector<Sparse_Matrix> penalty_B;
		std::vector<Sparse_Matrix> BC_Operator;	// the boundary operator in the target functional

		// by default we have a false
		bool bc_inhomo_time = false;
		bool have_force = false;
		double exact_target_value = 0; // exact value of the target
	};
