#include "run_problem.h"

using namespace dealii;
void develop_system(system_data &system_matrices);
void develop_system_adjoint(system_data &system_matrices);

Full_matrix compute_Amod(const Sparse_Matrix &A)
{
      EigenSolver<MatrixXd> ES(A);
      Full_matrix vecs = ES.pseudoEigenvectors();
      VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();

      Full_matrix Amod = vecs*vals.cwiseAbs().asDiagonal()*vecs.inverse();

      return(Amod);

}

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
      const int poly_degree = 0;
      std::string foldername = "2D_advection_gaussian";

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
      int repetitions = atoi(argv[2]);
      Triangulation<dim> triangulation;
      GridGenerator::subdivided_hyper_cube(triangulation,repetitions);
      triangulation.refine_global(1);		// we need to refine atleast once
      set_square_bid(triangulation);


      ic_bc<dim> initial_boundary;
      ic_bc_adjoint<dim> initial_boundary_adjoint;

	  run_problem<dim> Run_Problem(system_matrices,	  			 // system data
       								system_matrices,
				  			  		system_matrices_adjoint, 	// adjoint data
							  		triangulation, 				// triangulation
							  		poly_degree,
							  		&initial_boundary,
					          		&initial_boundary_adjoint,
					          		foldername);
     //}
}


void develop_system(system_data &system_matrices)
{

	// we can use the same sparsity pattern for all the test cases in this case
	system_matrices.Ax.resize(1,1);
	system_matrices.Ay.resize(1,1);
	system_matrices.P.resize(1,1);

	system_matrices.Ax.coeffRef(0,0) = 1;
	system_matrices.Ay.coeffRef(0,0) = 1;
	system_matrices.P.coeffRef(0,0) = 0;

	system_matrices.Ax.makeCompressed();
	system_matrices.Ay.makeCompressed();
	system_matrices.P.makeCompressed();
	
	//allocation of boundary matrices
	system_matrices.B.resize(4);
	system_matrices.penalty.resize(4);
	system_matrices.penalty_B.resize(4);


	for (unsigned int i = 0 ; i < 4 ; i ++)
	{
		system_matrices.B[i].resize(1,1);
		system_matrices.penalty[i].resize(1,1);
		system_matrices.penalty_B[i].resize(1,1);
	}
	//boundary at x = 1 
	unsigned int bc_id = 0;

	system_matrices.B[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty_B[bc_id].coeffRef(0,0) = 0;

	// top boundary  y = 1
	bc_id = 1;

	system_matrices.B[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty_B[bc_id].coeffRef(0,0) = 0;

	// left boundary x = 0
	bc_id = 2;

	system_matrices.B[bc_id].coeffRef(0,0) = 1;
	system_matrices.penalty[bc_id].coeffRef(0,0) = -1;
	system_matrices.penalty_B[bc_id].coeffRef(0,0) = -1;

	// boundary at y = 0, bottom boundary
	bc_id = 3;

	system_matrices.B[bc_id].coeffRef(0,0) = 1;
	system_matrices.penalty[bc_id].coeffRef(0,0) = -1;
	system_matrices.penalty_B[bc_id].coeffRef(0,0) = -1;

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

	// we can use the same sparsity pattern for all the test cases in this case
	system_matrices.Ax.resize(1,1);
	system_matrices.Ay.resize(1,1);
	system_matrices.P.resize(1,1);

	system_matrices.Ax.coeffRef(0,0) = -1;
	system_matrices.Ay.coeffRef(0,0) = -1;
	system_matrices.P.coeffRef(0,0) = 0;

	system_matrices.Ax.makeCompressed();
	system_matrices.Ay.makeCompressed();
	system_matrices.P.makeCompressed();
	
	//allocation of boundary matrices
	system_matrices.B.resize(4);
	system_matrices.penalty.resize(4);
	system_matrices.penalty_B.resize(4);


	for (unsigned int i = 0 ; i < 4 ; i ++)
	{
		system_matrices.B[i].resize(1,1);
		system_matrices.penalty[i].resize(1,1);
		system_matrices.penalty_B[i].resize(1,1);
	}
	//boundary at x = 1 
	unsigned int bc_id = 0;

	system_matrices.B[bc_id].coeffRef(0,0) = 1;
	system_matrices.penalty[bc_id].coeffRef(0,0) = -1;
	system_matrices.penalty_B[bc_id].coeffRef(0,0) = -1;

	// top boundary  y = 1
	bc_id = 1;

	system_matrices.B[bc_id].coeffRef(0,0) = 1;
	system_matrices.penalty[bc_id].coeffRef(0,0) = -1;
	system_matrices.penalty_B[bc_id].coeffRef(0,0) = -1;

	// left boundary x = 0
	bc_id = 2;

	system_matrices.B[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty_B[bc_id].coeffRef(0,0) = 0;

	// boundary at y = 0, bottom boundary
	bc_id = 3;

	system_matrices.B[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty_B[bc_id].coeffRef(0,0) = 0;

	// loop over all the boundaries
	for (unsigned int i = 0 ; i < 4 ; i ++)
	{
	system_matrices.B[i].makeCompressed();
	system_matrices.penalty[i].makeCompressed();
	system_matrices.penalty_B[i].makeCompressed();		
	}
}



template<int dim>
double 
ic_bc<dim>::ic(const Point<dim> &p,const int &id)
{
	const double x = p[0];
	const double y = p[1];

	return(0);
}

template<int dim>
void 
ic_bc<dim>::exact_solution(const Point<dim> &p,Vector<double> &value,const double &t)
{
	const double x = p[0];
	const double y = p[1];

	value(0) = exp(-pow((x-0.5),2)*100) * exp(-pow((y-0.5),2)*100);

}

template<int dim>
void 
ic_bc<dim>::force(const Point<dim> &p,Vector<double> &value,const double &t)
{
	value = 0;
}


template<int dim>
void
ic_bc<dim>::force(Vector<double> &value,
	  const Vector<double> &force_vec,
	  const Point<dim> &p,const double &t)
{
	const double x = p[0];
	const double y = p[1];

	value(0) = (-200.*(-1. + x + y))*exp(-100*(0.5 - x + pow(x,2) -y + pow(y,2)));
	//value(0) = M_PI * (cos(M_PI * x)*sin(M_PI*y)+cos(M_PI * y)*sin(M_PI*x));
}

template<int dim>
void 
ic_bc<dim>::bc_inhomo(const Sparse_Matrix &B,const unsigned int &bc_id,
								Vector<double> &value,const double &t)
{
	const int num_bc = B.rows();
	value.reinit(num_bc);
	value = 0;
}

template<int dim>
double 
ic_bc_adjoint<dim>::ic(const Point<dim> &p,const int &id)
{
	const double x = p[0];
	const double y = p[1];

	return(0);
}

template<int dim>
void 
ic_bc_adjoint<dim>::exact_solution(const Point<dim> &p,Vector<double> &value,const double &t)
{
	const double x = p[0];
	const double y = p[1];

	value(0) = 0;

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
	  const Point<dim> &p,const double &t)
{
	const double x = p[0];
	const double y = p[1];

	value(0) = (200.*(-1. + x + y))*exp(-100*(0.5 - x + pow(x,2) -y + pow(y,2)));
	//value(0) = -M_PI * (cos(M_PI * x)*sin(M_PI*y)+cos(M_PI * y)*sin(M_PI*x));
	
}

template<int dim>
void 
ic_bc_adjoint<dim>::bc_inhomo(const Sparse_Matrix &B,const unsigned int &bc_id,
								Vector<double> &value,const double &t)
{
	const int num_bc = B.rows();
	value.reinit(num_bc);
	value = 0;
}
