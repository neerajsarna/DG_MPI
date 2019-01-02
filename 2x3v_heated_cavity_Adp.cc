#include "run_problem.h"
#include "read_matrices.h"

using namespace dealii;


// given the number of equations and the number of boundary conditions, following routine develops system matrices for
// our ADJOINT moment system. 
void develop_system_adjoint(system_data &system_matrices,const int &M,
					const int &neqn_M,const int &nbc_M,const double &Kn);	// develop the adjiont system

// same as above but for the PRIMAL problem.
void develop_system(system_data &system_matrices,const int &M,
					const int &neqn_M,const int &nbc_M,const double &Kn);	// develop the primal system


// Same as above but takes in vectors of number of equation and rest of the parameters. Then 
// for every entry inside the vector, it calls the above routines. 
// Routine for PRIMAL problem.
std::vector<system_data>
develop_complete_system(const std::vector<int> &M,
						const std::vector<int> &neqn_M,
						const std::vector<int> &nbc_M,
						const double Kn);

// Routine for ADJOINT problem. 
std::vector<system_data>
develop_complete_system_adjoint(const std::vector<int> &M,
						const std::vector<int> &neqn_M,
						const std::vector<int> &nbc_M,
						const double Kn);

// class which contains the initial and boundary conditions
template<int dim>
class
ic_bc:public ic_bc_base<dim>
{
	public:
		ic_bc() {;};
		virtual double ic(const Point<dim> &p,const int &id);
		virtual void exact_solution(const Point<dim> &p,Vector<double> &value,const double &t);
				virtual void force(const Point<dim> &p,Vector<double> &value,const double &t);
		virtual void force(Vector<double> &value,
						   const Vector<double> &force_vec,
						   const Point<dim> &p,const double &t);
		virtual void bc_inhomo(const Sparse_Matrix &B,const unsigned int &bc_id,
								Vector<double> &value,const double &t);
};

// same as above but for the adjoint problem
template<int dim>
class
ic_bc_adjoint:public ic_bc_base<dim>
{
	public:
		ic_bc_adjoint() {;};
		virtual double ic(const Point<dim> &p,const int &id);
		virtual void exact_solution(const Point<dim> &p,Vector<double> &value,const double &t);
		virtual void force(const Point<dim> &p,Vector<double> &value,const double &t);
		virtual void force(Vector<double> &value,
						   const Vector<double> &force_vec,
						   const Point<dim> &p,const double &t);
		virtual void bc_inhomo(const Sparse_Matrix &B,const unsigned int &bc_id,
								Vector<double> &value,const double &t);
};

// compute the modulus of a matrix
Full_matrix compute_Amod(const Sparse_Matrix &A)
{
      EigenSolver<MatrixXd> ES(A);
      Full_matrix vecs = ES.pseudoEigenvectors();
      VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();

      Full_matrix Amod = vecs*vals.cwiseAbs().asDiagonal()*vecs.inverse();

      return(Amod);

}

// set boundary ids in a square
void 
set_square_bid(Triangulation<2> &triangulation)
{
	typename Triangulation<2>::active_cell_iterator
											 cell = triangulation.begin_active(),
											 endc = triangulation.end();

    const double left_edge = 0;
    const double right_edge = 1;

    for(; cell != endc ; cell++)
    	if(cell->is_locally_owned())
    	{
                for (unsigned int face = 0 ; face < GeometryInfo<2>::faces_per_cell ; face++)
              
                  if (cell->face(face)->at_boundary())
                  { 
                    double x_cord = cell->face(face)->center()(0);
                    double y_cord = cell->face(face)->center()(1);

                    // this boundary ID is coherent (-1) with the matlab code

                    // Following are the boundary ids:
                    // Right Wall = 0
                    // Top wall = 1
                    // left wall = 2
                    // bottom wall = 3
                    
                   // right
                    if (x_cord == right_edge)
                      cell->face(face)->set_boundary_id(0);

                    // top edge
                    if (y_cord == right_edge)
                      cell->face(face)->set_boundary_id(1);

                    // left edge
                    if (x_cord == left_edge)
                      cell->face(face)->set_boundary_id(2);

                    // bottom 
                    if (y_cord == left_edge)
                      cell->face(face)->set_boundary_id(3);
                   }
    	}
}



int main(int argc, char *argv[])
{
      using namespace dealii;
      
      const unsigned int num_threads = atoi(argv[1]);
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, num_threads);

      const int dim = 2;
      const int dim_problem = 2;
      const int poly_degree = 0;

     
     // store the number of equations for a given M (starts at M = 3)
     const std::vector<int> neqn_M = {13,22,34,50,70,95,125,161,203,252,308,372,444,525,615,715,825,946};
     const std::vector<int> nbc_M = {5,8,14,20,30,40,55,70,91,112,140,168,204,240,285,330,385,440};
     const double Kn = 0.1;

     const unsigned int num_systems = 4;
     std::vector<int> M(num_systems);
     std::vector<int> M_adjoint(num_systems);
     M[0] = 3;
     M[1] = 5;
     M[2] = 7;
     M[3] = 9;

     M_adjoint[0] = 3;
     M_adjoint[1] = 5;
     M_adjoint[2] = 7;
     M_adjoint[3] = 9;

     // develop system matrices for the primal problem.
	 std::vector<system_data> system_matrices = develop_complete_system(M,neqn_M,nbc_M,Kn);
	 // develop system matrices for the error computation. System matrices are same as the primal problem
	 // but are bigger for error computation.
	 std::vector<system_data> system_matrices_error = develop_complete_system(M_adjoint,neqn_M,nbc_M,Kn);
	 // system matrices for adjoint. 
     std::vector<system_data> system_matrices_adjoint = develop_complete_system_adjoint(M_adjoint,neqn_M,nbc_M,Kn);
     
      // create a rectangular mesh 
      Triangulation<dim> triangulation;
      // number of repetitions in every direction
      unsigned int repetitions = atoi(argv[2]);

      GridGenerator::subdivided_hyper_cube(triangulation,repetitions);
      triangulation.refine_global(1);
      set_square_bid(triangulation);

      ic_bc<dim> initial_boundary;	
      ic_bc_adjoint<dim> initial_boundary_adjoint;	

      // folder where all the output goes
      std::string foldername = "2x3v_heated_cavity_Adp/";

      // Maximum number of equations in all the different moment systems we loaded. Following computation assumes that
      // entries inside system_matrices have increasing value of $M$ (maximum tensorial degree in the Hermite expansion).
      const unsigned int max_neqn_primal = system_matrices[system_matrices.size()-1].Ax.rows();
      const unsigned int max_neqn_adjoint = system_matrices_adjoint[system_matrices_adjoint.size()-1].Ax.rows();


      // we give the exact value of the target functional from the DVM solution. Sufficient to put the value in the 
      // first element.
      // target functional here is the total heat flux in the y-direction i.e. \int_{(0,1)\times (0,1)}q_y dx dy.
      system_matrices[0].exact_target_value = 0.154185;


      run_problem<dim> Run_Problem(system_matrices,	  // system data
       								system_matrices_error,
				  			  		system_matrices_adjoint, // adjoint data
							  		triangulation, // triangulation
							  		poly_degree,
							  		&initial_boundary,
					          		&initial_boundary_adjoint,
					          		foldername,
					          		max_neqn_primal,
					          		max_neqn_adjoint,
					          		dim_problem);

}


void develop_system(system_data &system_matrices,const int &M,const int &neqn_M,
												 const int &nbc_M,const double &Kn)
{
	std::cout << "developing systems: " << std::endl;
	// we first initialise all the matrices
	system_matrices.Ax.resize(neqn_M,neqn_M);
	system_matrices.Ay.resize(neqn_M,neqn_M);
	system_matrices.P.resize(neqn_M,neqn_M);

	system_matrices.B.resize(4);
	system_matrices.penalty_B.resize(4);
	system_matrices.penalty.resize(4);
	system_matrices.BC_Operator.resize(4);

	for(unsigned int id = 0 ; id < 4 ; id ++)
	{
		system_matrices.B[id].resize(nbc_M,neqn_M);
		system_matrices.penalty_B[id].resize(neqn_M,neqn_M);
		system_matrices.penalty[id].resize(neqn_M,nbc_M);
		system_matrices.BC_Operator[id].resize(neqn_M,nbc_M);
	}

	std::vector<triplet> Row_Col_Value;
	std::string filename = "3v_Moments/Ax/Ax" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.Ax,Row_Col_Value);

	filename = "3v_Moments/P/P" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.P,Row_Col_Value);
	system_matrices.P = system_matrices.P/Kn;
	system_matrices.P.makeCompressed();

	filename = "3v_Moments/Bwall/Bwall" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.B[0],Row_Col_Value);

	filename = "3v_Moments/Bwall/penalty_odd_wall" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.penalty[0],Row_Col_Value);	
	
	std::vector<Sparse_Matrix> rotator(4);

	for(unsigned int id = 0 ; id < 4 ; id++ )
	{
		rotator[id].resize(neqn_M,neqn_M);
		filename = "3v_Moments/Rotator/rotator" + std::to_string(M)
					 + '_' + std::to_string(id+1) + ".txt";
					 
		build_triplet(Row_Col_Value,filename);
		build_matrix_from_triplet(rotator[id],Row_Col_Value);		
	}

	system_matrices.Ax.makeCompressed();
	system_matrices.Ay = rotator[1].transpose() * system_matrices.Ax * rotator[1];
	system_matrices.Ay.makeCompressed();

	for(unsigned int i = 0 ; i < 4 ; i++)
	{
		system_matrices.B[i] = system_matrices.B[0] * rotator[i];
		system_matrices.penalty[i] = rotator[i].transpose() * system_matrices.penalty[0];
		system_matrices.BC_Operator[i] =  system_matrices.penalty[i];
		system_matrices.penalty_B[i] = system_matrices.penalty[i] * system_matrices.B[i];

		system_matrices.B[i].makeCompressed();
		system_matrices.penalty[i].makeCompressed();
		system_matrices.BC_Operator[i].makeCompressed();
		system_matrices.penalty_B[i].makeCompressed();
	}

}


void develop_system_adjoint(system_data &system_matrices,const int &M,const int &neqn_M,
							const int &nbc_M,const double &Kn)
{
	// first we develop a temporary primal system. Then, by changing matrices in this system
	// we develop the adjoint system. 
	system_data temp;
	develop_system(temp,M,neqn_M,nbc_M,Kn);

	// allocate memory
	system_matrices.Ax.resize(neqn_M,neqn_M);
	system_matrices.Ay.resize(neqn_M,neqn_M);
	system_matrices.P.resize(neqn_M,neqn_M);

	system_matrices.B.resize(4);
	system_matrices.penalty_B.resize(4);
	system_matrices.penalty.resize(4);
	system_matrices.BC_Operator.resize(4);

	for(unsigned int id = 0 ; id < 4 ; id ++)
	{
		system_matrices.B[id].resize(nbc_M,neqn_M);
		system_matrices.penalty_B[id].resize(neqn_M,neqn_M);
		system_matrices.penalty[id].resize(neqn_M,nbc_M);
		system_matrices.BC_Operator[id].resize(neqn_M,nbc_M); // the boundary operator for computation of adjiont functional
	}

	// reverse the direction of advection
	system_matrices.Ax = - temp.Ax;
	system_matrices.Ay = - temp.Ay;

	system_matrices.Ax.makeCompressed();
	system_matrices.Ay.makeCompressed();

	// production remains the same
	system_matrices.P = temp.P;
	system_matrices.P.makeCompressed();

	// boundary conditions also flip signs
	std::vector<int> bc_id_primal(4); // the id of primal which is the id of adjoint
	bc_id_primal[0] = 2;	// adjoint boundary at x = 1 is the primal boundary at x = 0 (reversal in advection direction)
	bc_id_primal[1] = 3;
	bc_id_primal[2] = 0;
	bc_id_primal[3] = 1;

	for(unsigned int i = 0 ; i < 4 ; i++)
	{
		system_matrices.B[i] = temp.B[bc_id_primal[i]];
		system_matrices.penalty_B[i] = temp.penalty_B[bc_id_primal[i]];
		system_matrices.penalty[i] = temp.penalty[bc_id_primal[i]];
		system_matrices.BC_Operator[i] = temp.BC_Operator[bc_id_primal[i]];

		system_matrices.B[i].makeCompressed();
		system_matrices.penalty_B[i].makeCompressed();
		system_matrices.penalty[i].makeCompressed();
		system_matrices.BC_Operator[i].makeCompressed();
	}
}

std::vector<system_data>
develop_complete_system(const std::vector<int> &M,
						const std::vector<int> &neqn_M,
						const std::vector<int> &nbc_M,
						const double Kn)
{
     std::vector<system_data> system_matrices(M.size());
     Assert(*std::max_element(M.begin(),M.end())<=20,ExcNotImplemented());

     for(unsigned int i = 0 ; i < M.size(); i++)
     {
     	develop_system(system_matrices[i],M[i],neqn_M[M[i]-3],nbc_M[M[i]-3],Kn);
     	system_matrices[i].bc_inhomo_time = true;

     	system_matrices[i].Ax_mod = compute_Amod(system_matrices[i].Ax).sparseView();
      	system_matrices[i].Ay_mod = compute_Amod(system_matrices[i].Ay).sparseView();
      	system_matrices[i].Ax_mod.makeCompressed();
      	system_matrices[i].Ay_mod.makeCompressed();
     }

     return(system_matrices);
}

std::vector<system_data>
develop_complete_system_adjoint(const std::vector<int> &M,
						const std::vector<int> &neqn_M,
						const std::vector<int> &nbc_M,
						const double Kn)
{
     std::vector<system_data> system_matrices(M.size());
     Assert(*std::max_element(M.begin(),M.end())<=20,ExcNotImplemented());

     for(unsigned int i = 0 ; i < M.size(); i++)
     {
     	develop_system_adjoint(system_matrices[i],M[i],neqn_M[M[i]-3],nbc_M[M[i]-3],Kn);
     	system_matrices[i].bc_inhomo_time = false;
     	system_matrices[i].have_force = true;

     	system_matrices[i].Ax_mod = compute_Amod(system_matrices[i].Ax).sparseView();
      	system_matrices[i].Ay_mod = compute_Amod(system_matrices[i].Ay).sparseView();
      	system_matrices[i].Ax_mod.makeCompressed();
      	system_matrices[i].Ay_mod.makeCompressed();
     }

     return(system_matrices);
}


template<int dim>
double 
ic_bc<dim>::ic(const Point<dim> &p,const int &id)
{
	const double x = p[0];
	const double y = p[1];

	const double density = 0;	
	// const double density = exp(-pow((x-0.5),2)*100);	

	switch (id)
	{
		case 0:
		{
			return(density);
			break;
		}
		default:
		{
			return(0);
			break;
		}
	}
	
	return 0;
}

template<int dim>
void 
ic_bc<dim>::exact_solution(const Point<dim> &p,Vector<double> &value,const double &t)
{
	const double x = p[0];
	const double y = p[1];

	value(0) = 0;

}

template<int dim>
void 
ic_bc<dim>::force(const Point<dim> &p,Vector<double> &value,const double &t)
{
	const double x = p[0];
	const double y = p[1];

	Assert(value.size() != 0 ,ExcNotImplemented());
	value(0) = 0;

}

template<int dim>
void 
ic_bc<dim>::force(Vector<double> &value,
				  		 const Vector<double> &force_vec,
				  		 const Point<dim> &p,const double &t)
{
	Assert(value.size() != 0,ExcNotImplemented());
	value = 0;
}



template<int dim>
void 
ic_bc<dim>::bc_inhomo(const Sparse_Matrix &B,const unsigned int &bc_id,
								Vector<double> &value,const double &t)
{

	const int num_bc = B.rows();
	double thetaW;
	value.reinit(num_bc);
	value = 0;

	if (t <= 1)
        thetaW = exp(-1/(1-pow((t-1),2))) * exp(1);
    else
        thetaW = 1;

    
	switch (bc_id)
	{
		case 0:
		case 1:
		case 2:
		{
			break;
		}

		case 3:
		{

			double bc_normal = 1.0;
			for (unsigned int m = 0 ; m < B.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(B,m); n ; ++n)
                    	if (n.col() == 3 || n.col() == 5 || n.col() == 6)
                    		value(n.row()) += bc_normal * thetaW * n.value()/sqrt(2.0);
                    	
			break;
		}

		default:
		{
			Assert(1 == 0, ExcMessage("should not have reached"));
			break;
		}
	}

}

template<int dim>
double 
ic_bc_adjoint<dim>::ic(const Point<dim> &p,const int &id)
{	
	return 0;
}

template<int dim>
void 
ic_bc_adjoint<dim>::exact_solution(const Point<dim> &p,Vector<double> &value,const double &t)
{
	value = 0;

}



template<int dim>
void 
ic_bc_adjoint<dim>::bc_inhomo(const Sparse_Matrix &B,const unsigned int &bc_id,
								Vector<double> &value,const double &t)
{

	const int num_bc = B.rows();
	double thetaW;
	value.reinit(num_bc);
	value = 0;

}

template<int dim>
void 
ic_bc_adjoint<dim>::force(const Point<dim> &p,Vector<double> &value,const double &t)
{
	Assert(value.size() != 0,ExcNotImplemented());
	value = 0;

}


template<int dim>
void 
ic_bc_adjoint<dim>::force(Vector<double> &value,
				  		 const Vector<double> &force_vec,
				  		 const Point<dim> &p,
				  		 const double &t)
{
	Assert(value.size() != 0,ExcNotImplemented());
	Assert(force_vec.size() != 0,ExcNotImplemented());
	value = 0;

	// mean value for qy 
	value(11) = sqrt(3.0/2.0); // coefficient for (0,3,0)
	value(8) = sqrt(1.0/2.0); // coefficient for (2,1,0)
	value(12) = sqrt(1.0/2.0); // coefficient for (0,1,2)
}

