#include "run_problem.h"
#include "read_matrices.h"

using namespace dealii;

void develop_system(system_data &system_matrices);
void develop_system_adjoint(system_data &system_matrices);


std::vector<system_data>
develop_complete_system(const std::vector<int> &M,
						const std::vector<int> &neqn_M,
						const std::vector<int> &nbc_M,
						const double Kn);

std::vector<system_data>
develop_complete_system_adjoint(const std::vector<int> &M,
						const std::vector<int> &neqn_M,
						const std::vector<int> &nbc_M,
						const double Kn);

// now we specify the iniital and the boundary conditions
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

Full_matrix compute_Amod(const Sparse_Matrix &A)
{
      EigenSolver<MatrixXd> ES(A);
      Full_matrix vecs = ES.pseudoEigenvectors();
      VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();

      Full_matrix Amod = vecs*vals.cwiseAbs().asDiagonal()*vecs.inverse();

      return(Amod);

}

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
      const int dim_problem = 1;
      const int poly_degree = 0;
      std::string foldername = "1D_wave";

     std::vector<system_data> system_matrices(1);
     std::vector<system_data> system_matrices_adjoint(1);

     develop_system(system_matrices[0]);
     develop_system_adjoint(system_matrices_adjoint[0]);

     system_matrices[0].have_force = true;
     system_matrices[0].Ax_mod = compute_Amod(system_matrices[0].Ax).sparseView();
     system_matrices[0].Ay_mod = compute_Amod(system_matrices[0].Ay).sparseView();
 
     system_matrices_adjoint[0].have_force = true;
     system_matrices_adjoint[0].Ax_mod = compute_Amod(system_matrices_adjoint[0].Ax).sparseView();
     system_matrices_adjoint[0].Ay_mod = compute_Amod(system_matrices_adjoint[0].Ay).sparseView();

	  // develop mesh
     //for(int i = 1 ; i < 30 ; i++ ) 
     //{

      std::vector<unsigned int> repetitions(dim);

      const double left_edge = 0;
      const double right_edge = 1;
      Point<dim> p1;
      Point<dim> p2;

      // corners of the diagonal
      p1(0) = left_edge;
      p1(1) = left_edge;

      p2(0) = right_edge;
      p2(1) = right_edge;

      repetitions[0] = atoi(argv[2]);
      repetitions[1] = dim_problem == 1 ?  1 : repetitions[0];

      Triangulation<dim> triangulation;
      GridGenerator::subdivided_hyper_rectangle(triangulation,repetitions,p1,p2);

      typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

      for(; cell != endc ; cell++)
      		cell->set_refine_flag(RefinementCase<dim>::cut_axis(0));

      triangulation.execute_coarsening_and_refinement();
      set_square_bid(triangulation);

      ic_bc<dim> initial_boundary;
      ic_bc_adjoint<dim> initial_boundary_adjoint;

      const unsigned int max_neqn_primal = system_matrices[system_matrices.size()-1].Ax.rows();
      const unsigned int max_neqn_adjoint = system_matrices_adjoint[system_matrices_adjoint.size()-1].Ax.rows();

	  run_problem<dim> Run_Problem(system_matrices,	  			 // system data
       								system_matrices,
				  			  		system_matrices_adjoint, 	// adjoint data
							  		triangulation, 				// triangulation
							  		poly_degree,
							  		&initial_boundary,
					          		&initial_boundary_adjoint,
					          		foldername,
					          		max_neqn_primal,
					          		max_neqn_adjoint,
					          		dim_problem);
}


void develop_system(system_data &system_matrices)
{

	const unsigned int n_eqn = 2;

	// we can use the same sparsity pattern for all the test cases in this case
	system_matrices.Ax.resize(n_eqn,n_eqn);
	system_matrices.Ay.resize(n_eqn,n_eqn);
	system_matrices.P.resize(n_eqn,n_eqn);

	// for(unsigned int i = 0 ; i < n_eqn ; i++)
	// 	system_matrices.Ax.coeffRef(i,i) = 1;	

	system_matrices.Ax.coeffRef(0,1) = 1;
	system_matrices.Ax.coeffRef(1,0) = 1;

	// system_matrices.P.coeffRef(0,0) = -10;
	// system_matrices.P.coeffRef(1,1) = -10;

	system_matrices.Ax.makeCompressed();
	system_matrices.Ay.makeCompressed();
	system_matrices.P.makeCompressed();
	
	//allocation of boundary matrices
	system_matrices.B.resize(4);
	system_matrices.penalty.resize(4);
	system_matrices.penalty_B.resize(4);


	for (unsigned int i = 0 ; i < 4 ; i ++)
	{
		system_matrices.B[i].resize(1,n_eqn);
		system_matrices.penalty[i].resize(n_eqn,1);
		system_matrices.penalty_B[i].resize(n_eqn,n_eqn);
	}

	// boundary conditions for wave equation
	// right boundary x = 0
	unsigned int bc_id = 0;
	system_matrices.B[bc_id].coeffRef(0,0) = 1;
	system_matrices.B[bc_id].coeffRef(0,1) = -1;

	system_matrices.penalty[bc_id].coeffRef(0,0) = -0.5;
	system_matrices.penalty[bc_id].coeffRef(1,0) = 0.5;

	system_matrices.penalty_B[bc_id] = system_matrices.penalty[bc_id] * system_matrices.B[bc_id];

	bc_id = 2;
	system_matrices.B[bc_id].coeffRef(0,0) = 1;
	system_matrices.B[bc_id].coeffRef(0,1) = 1;

	system_matrices.penalty[bc_id].coeffRef(0,0) = -0.5;
	system_matrices.penalty[bc_id].coeffRef(1,0) = -0.5;

	system_matrices.penalty_B[bc_id] = system_matrices.penalty[bc_id] * system_matrices.B[bc_id];


	// boundary conditions for advection equation
	// unsigned int bc_id = 2;
	// for(unsigned int i = 0 ; i < n_eqn ; i++)
	// system_matrices.B[bc_id].coeffRef(i,i) = 1;

	// for(unsigned int i = 0 ; i < n_eqn ; i++)
	// system_matrices.penalty[bc_id].coeffRef(i,i) = -1;

	// for(unsigned int i = 0 ; i < n_eqn ; i++)
	// system_matrices.penalty_B[bc_id].coeffRef(i,i) = -1;

	// loop over all the boundaries
	for (unsigned int i = 0 ; i < 4 ; i ++)
	{
	system_matrices.B[i].makeCompressed();
	system_matrices.penalty[i].makeCompressed();
	system_matrices.penalty_B[i].makeCompressed();		
	}
}


void develop_system_adjoint(system_data &system_matrices)
{

	unsigned int n_eqn = 2;

	system_data temp;
	develop_system(temp);

	// allocate memory
	system_matrices.Ax.resize(n_eqn,n_eqn);
	system_matrices.Ay.resize(n_eqn,n_eqn);
	system_matrices.P.resize(n_eqn,n_eqn);

	system_matrices.Ax.makeCompressed();
	system_matrices.Ay.makeCompressed();
	system_matrices.P.makeCompressed();
	
	//allocation of boundary matrices
	system_matrices.B.resize(4);
	system_matrices.penalty.resize(4);
	system_matrices.penalty_B.resize(4);


	for (unsigned int i = 0 ; i < 4 ; i ++)
	{
		system_matrices.B[i].resize(1,n_eqn);
		system_matrices.penalty[i].resize(n_eqn,1);
		system_matrices.penalty_B[i].resize(n_eqn,n_eqn);
	}

	// reverse the direction of advection
	system_matrices.Ax = - temp.Ax;

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

		system_matrices.B[i].makeCompressed();
		system_matrices.penalty_B[i].makeCompressed();
		system_matrices.penalty[i].makeCompressed();
	}

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
ic_bc<dim>::force(const Point<dim> &p,Vector<double> &value,const double &t)
{
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
	const double x = p[0];
	
	value(0) = M_PI * cos(M_PI * x);//+ exp(-pow((x-0.5),2)*100) * 10;
	value(1) = -200 * exp(-100*pow((x-0.5),2)) * (x-0.5);// + sin(M_PI * x) * 10;

}



template<int dim>
void 
ic_bc<dim>::bc_inhomo(const Sparse_Matrix &B,const unsigned int &bc_id,
								Vector<double> &value,const double &t)
{

	const int num_bc = B.rows();
	double thetaW = 0;
	value.reinit(num_bc);
	value = 0;


}

template<int dim>
double 
ic_bc_adjoint<dim>::ic(const Point<dim> &p,const int &id)
{
	return(0);
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
	value = 0;
}


template<int dim>
void 
ic_bc_adjoint<dim>::force(Vector<double> &value,
				  		 const Vector<double> &force_vec,
				  		 const Point<dim> &p,
				  		 const double &t)
{
	Assert(value.size() != 0,ExcNotInitialized());
	Assert(force_vec.size() != 0,ExcNotInitialized());
	const double x = p[0];
	value = 0;

	// for(unsigned int i = 0 ; i < value.size() ; i++)
	// 	value(i) = 1;
	value(0) = exp(-pow((x-1),2)*100);
}


template<int dim>
void 
ic_bc<dim>::exact_solution(const Point<dim> &p,Vector<double> &value,const double &t)
{
	const double x = p[0]; // we need to shift the x coordinate for the reference solution
	const double y = p[1];


	value = 0;

   	
   	value(0) = exp(-pow((x-0.5),2)*100);
	value(1) = sin(M_PI * x);

}
