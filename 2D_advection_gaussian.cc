#include "run_problem.h"

using namespace dealii;
void develop_system(system_data &system_matrices);
void develop_system_adjoint(system_data &system_matrices);
//template<int dim> void set_boundary_id(Triangulation<dim> &triangulation);


// compute the derivative of exp(-(x-mu1)^2*var1-(y-mu1)^2*var1)
double 
dxdy_bi_modal(const double &mu,const double &var,
			  const double &x,const double &y)
{
	return(-2 * var * exp(-(pow(x-mu,2)+pow(y-mu,2))*var) * (x-mu + y-mu));
}	

Full_matrix compute_Amod(const Sparse_Matrix &A)
{
      EigenSolver<MatrixXd> ES(A);
      Full_matrix vecs = ES.pseudoEigenvectors();
      VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();

      Full_matrix Amod = vecs*vals.cwiseAbs().asDiagonal()*vecs.inverse();

      return(Amod);

}

// now we specify the inital and the boundary conditions
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
	typename Triangulation<2>::active_cell_iterator cell = triangulation.begin_active(),
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
      {
      	if(dim_problem == 1)
      	{
      		cell->set_refine_flag(RefinementCase<dim>::cut_axis(0));
      	}
      	else
      	{
      		cell->set_refine_flag();
      	}
      }
      triangulation.execute_coarsening_and_refinement();

      set_square_bid(triangulation);


      ic_bc<dim> initial_boundary;
      ic_bc_adjoint<dim> initial_boundary_adjoint;

      const unsigned int max_neqn_primal = system_matrices[system_matrices.size()-1].Ax.rows();
      const unsigned int max_neqn_adjoint = system_matrices_adjoint[system_matrices_adjoint.size()-1].Ax.rows();


      if(dim_problem == 1)
      {
      	system_matrices[0].Ay.coeffRef(0,0) = 0;
      	system_matrices[0].Ay.makeCompressed();

      	system_matrices[0].B[3].coeffRef(0,0) = 0;
      	system_matrices[0].B[3].makeCompressed();
		system_matrices[0].penalty[3].coeffRef(0,0) = 0;
		system_matrices[0].penalty[3].makeCompressed();
		system_matrices[0].penalty_B[3].coeffRef(0,0) = 0;
		system_matrices[0].penalty_B[3].makeCompressed();

      	system_matrices_adjoint[0].B[1].coeffRef(0,0) = 0;
      	system_matrices_adjoint[0].B[1].makeCompressed();
		system_matrices_adjoint[0].penalty[1].coeffRef(0,0) = 0;
		system_matrices_adjoint[0].penalty[1].makeCompressed();
		system_matrices_adjoint[0].penalty_B[1].coeffRef(0,0) = 0;
		system_matrices_adjoint[0].penalty_B[1].makeCompressed();
      }

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
	//system_matrices.Ay.coeffRef(0,0) = 0;
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
	//system_matrices.Ay.coeffRef(0,0) = 0;
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
	Vector<double> temp_p(2);

	for(unsigned int space = 0 ; space < dim; space++)
		temp_p(space) = p[space];

	const double x = temp_p[0];
	const double y = temp_p[1];

	return(0);
}

template<int dim>
void 
ic_bc<dim>::exact_solution(const Point<dim> &p,Vector<double> &value,const double &t)
{
	Vector<double> temp_p(2);

	for(unsigned int space = 0 ; space < dim; space++)
		temp_p(space) = p[space];

	const double x = temp_p[0];
	const double y = temp_p[1];

	const double mu1 = 0.5;
	const double mu2 = 0.8;
	const double var1 = 100;
	const double var2 = 50;

	//value(0) = sin(M_PI * x);
	value(0) = (exp(-(pow((x-mu1),2)+pow((y-mu1),2))*var1) 
				+ exp(-(pow((x-mu2),2)+pow((y-mu2),2))*var2)) * 0.5;
	//value(0) = exp(-pow((x-1),2)*100);
	//value(0) = sin(M_PI * x) * sin(M_PI * y);

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
	Vector<double> temp_p(2);

	for(unsigned int space = 0 ; space < dim; space++)
		temp_p(space) = p[space];

	const double x = temp_p[0];
	const double y = temp_p[1];

	const double mu1 = 0.5;
	const double mu2 = 0.8;
	const double var1 = 100;
	const double var2 = 50;


	//value(0) = (-200.*(-1. + x + y))*exp(-100*(0.5 - x + pow(x,2) -y + pow(y,2)));
	//value(0) = M_PI * (cos(M_PI * x)*sin(M_PI*y)+cos(M_PI * y)*sin(M_PI*x));
	//value(0) = -100 * exp(-100 * pow(x-0.5,2)) * (-1 + 2 * x);
	//value(0) = M_PI * cos(M_PI * x);
	//value(0) = -200 * exp(-100*pow((x-1),2)) * (x-1);
	value(0) = 0.5 * (dxdy_bi_modal(mu1,var1,x,y) + dxdy_bi_modal(mu2,var2,x,y));
}

template<int dim>
void 
ic_bc<dim>::bc_inhomo(const Sparse_Matrix &B,const unsigned int &bc_id,
								Vector<double> &value,const double &t)
{
	const int num_bc = B.rows();
	value.reinit(num_bc);
	value = 0;

	if(bc_id == 2)
		value(0) = 0;
}

template<int dim>
double 
ic_bc_adjoint<dim>::ic(const Point<dim> &p,const int &id)
{
	Vector<double> temp_p(2);

	for(unsigned int space = 0 ; space < dim; space++)
		temp_p(space) = p[space];

	const double x = temp_p[0];
	const double y = temp_p[1];

	return(0);
}

template<int dim>
void 
ic_bc_adjoint<dim>::exact_solution(const Point<dim> &p,Vector<double> &value,const double &t)
{
	Vector<double> temp_p(2);

	for(unsigned int space = 0 ; space < dim; space++)
		temp_p(space) = p[space];

	const double x = temp_p[0];
	const double y = temp_p[1];

	//value(0) = exp(-pow((x-0.5),2)*100) * exp(-pow((y-0.5),2)*100);
	//value(0) = exp(-pow((x-0.5),2)*100);
	//value(0) = sin(M_PI * x) * sin(M_PI * y);
	value(0) = sin(M_PI * x);

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
	Vector<double> temp_p(2);

	for(unsigned int space = 0 ; space < dim; space++)
		temp_p(space) = p[space];

	const double x = temp_p[0];
	const double y = temp_p[1];

	//value(0) = (200.*(-1. + x + y))*exp(-100*(0.5 - x + pow(x,2) -y + pow(y,2)));
	//value(0) = (200.*( x - y))*exp(-100*(0.5 - x + pow(x,2) -y + pow(y,2)));
	//value(0) = -M_PI * (cos(M_PI * x)*sin(M_PI*y)+cos(M_PI * y)*sin(M_PI*x));
	//value(0) = 100 * exp(-100 * pow(x-0.5,2)) * (-1 + 2 * x);
	//value(0) = -M_PI * cos(M_PI * x);

	//value(0) = exp(-pow((x-0.8),2)*100);
	value(0) = 1;
}

template<int dim>
void 
ic_bc_adjoint<dim>::bc_inhomo(const Sparse_Matrix &B,const unsigned int &bc_id,
								Vector<double> &value,const double &t)
{
	const int num_bc = B.rows();
	value.reinit(num_bc);
	value = 0;

	if(bc_id == 0)
		value(0) = 0;
}
