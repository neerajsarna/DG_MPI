#include "run_problem.h"
#include "read_matrices.h"

using namespace dealii;

void develop_system_adjoint(system_data &system_matrices,const int &M,
					const int &neqn_M,const int &nbc_M,const double &Kn);	// develop the adjiont system

void develop_system(system_data &system_matrices,const int &M,
					const int &neqn_M,const int &nbc_M,const double &Kn);	// develop the primal system


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

void count_bc_id(Triangulation<2> &triangulation)
{
	  typename Triangulation<2>::active_cell_iterator cell = triangulation.begin_active(),
	  													 endc = triangulation.end();


	  
	  int id0 = 0,id1 = 0, id2 = 0, id3 = 0; 

	  for(; cell != endc ; cell++)   
	  	if (cell->is_locally_owned())
	  {

	  	for(unsigned int face = 0 ; face < GeometryInfo<2>::faces_per_cell ; face++)
	  		if (cell->face(face)->at_boundary())
	  			switch (cell->face(face)->boundary_id())
	  			{
	  				case 0:
	  				{
	  					id0++;
	  					break;
	  				}

	  				case 1:
	  				{
	  					id1++;
	  					break;
	  				}

	  				case 2:
	  				{
	  					id2++;
	  					break;
	  				}

	  				case 3:
	  				{
	  					id3++;
	  					break;
	  				}
	  			}
	  		}
	  		
	  

	  	  	std::cout << "id0 " << id0 
	  			<< " id1 " << id1 
	  			<< " id2 " << id2 
	  			<< " id3 " << id3 << std::endl;

	  		fflush(stdout);


}


int main(int argc, char *argv[])
{
      using namespace dealii;
      
      const unsigned int num_threads = atoi(argv[1]);
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, num_threads);

      const int dim = 2;
      const int poly_degree = 0;

     
     // store the number of equations for a given M (starts at M = 3)
     std::vector<int> neqn_M;
     std::vector<int> nbc_M;
     const double Kn = 0.1;

     const unsigned int num_systems = 5;
     std::vector<int> M(num_systems);
     std::vector<int> M_adjoint(num_systems);
     M[0] = 6;
     M[1] = 14;
     M[2] = 16;
     M[3] = 18;
     M[4] = 20;

     M_adjoint[0] = 6;
     M_adjoint[1] = 8;
     M_adjoint[2] = 10;
     M_adjoint[3] = 11;
     M_adjoint[4] = 12;
	
     neqn_M.resize(18);
     nbc_M.resize(18);

     for(unsigned int i = 3 ; i <= 20; i++)
     {
     	neqn_M[i-3]=i;
     	nbc_M[i-3]=int(i/2);
     }

	 std::vector<system_data> system_matrices = develop_complete_system(M,neqn_M,nbc_M,Kn);
	 std::vector<system_data> system_matrices_error = develop_complete_system(M_adjoint,neqn_M,nbc_M,Kn);	// required for error computation
     std::vector<system_data> system_matrices_adjoint = develop_complete_system_adjoint(M_adjoint,neqn_M,nbc_M,Kn);
     
      // create a rectangular mesh 
      Triangulation<dim> triangulation;
      Point<dim> p1;
      Point<dim> p2;
      std::vector<unsigned int> repetitions(dim);

      const double left_edge = 0;
      const double right_edge = 1;

      // corners of the diagonal
      p1(0) = left_edge;
      p1(1) = left_edge;

      p2(0) = right_edge;
      p2(1) = right_edge;

     
      repetitions[0] = atoi(argv[2]);
      repetitions[1] = 1;

      //The diagonal of the rectangle is the line joining p1 and p2
      GridGenerator::subdivided_hyper_rectangle(triangulation,repetitions,p1,p2);
      set_square_bid(triangulation);

     
      triangulation.signals.post_refinement.connect(std::bind (&set_square_bid,
                      								std::ref(triangulation)));

      ic_bc<dim> initial_boundary;	
      ic_bc_adjoint<dim> initial_boundary_adjoint;	

      std::string foldername = "2x1v_moments_Inflow_Adp/temp";

       run_problem<dim> Run_Problem(system_matrices,	  // system data
       								system_matrices_error,
				  			  		system_matrices_adjoint, // adjoint data
							  		triangulation, // triangulation
							  		poly_degree,
							  		&initial_boundary,
					          		&initial_boundary_adjoint,
					          		foldername);

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

	for(unsigned int id = 0 ; id < 4 ; id ++)
	{
		system_matrices.B[id].resize(nbc_M,neqn_M);
		system_matrices.penalty_B[id].resize(neqn_M,neqn_M);
		system_matrices.penalty[id].resize(neqn_M,nbc_M);
	}

	std::vector<triplet> Row_Col_Value;
	std::string filename = "1v_Moments/Ax/Ax" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.Ax,Row_Col_Value);

	filename = "1v_Moments/P/P" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.P,Row_Col_Value);
	system_matrices.P = system_matrices.P/Kn;
	system_matrices.P.makeCompressed();

	filename = "1v_Moments/Bwall/Bwall" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.B[0],Row_Col_Value);

	filename = "1v_Moments/Bwall/penalty_wall" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.penalty[0],Row_Col_Value);	
	
	std::vector<Sparse_Matrix> rotator(4);

	for(unsigned int id = 0 ; id < 4 ; id=id+2)
	{
		rotator[id].resize(neqn_M,neqn_M);
		filename = "1v_Moments/Rotator/rotator" + std::to_string(M)
					 + '_' + std::to_string(id+1) + ".txt";
					 
		build_triplet(Row_Col_Value,filename);
		build_matrix_from_triplet(rotator[id],Row_Col_Value);		
	}

	system_matrices.Ax.makeCompressed();

	for(unsigned int i = 0 ; i < 4 ; i=i+2)
	{
		system_matrices.B[i] = system_matrices.B[0] * rotator[i];
		system_matrices.penalty[i] = rotator[i].transpose() * system_matrices.penalty[0];

		system_matrices.penalty_B[i] = system_matrices.penalty[i] * system_matrices.B[i];

		system_matrices.B[i].makeCompressed();
		system_matrices.penalty[i].makeCompressed();
		system_matrices.penalty_B[i].makeCompressed();
	}

}


void develop_system_adjoint(system_data &system_matrices,const int &M,const int &neqn_M,
							const int &nbc_M,const double &Kn)
{
	system_data temp;
	develop_system(temp,M,neqn_M,nbc_M,Kn);

	// allocate memory
	system_matrices.Ax.resize(neqn_M,neqn_M);
	system_matrices.Ay.resize(neqn_M,neqn_M);
	system_matrices.P.resize(neqn_M,neqn_M);

	system_matrices.B.resize(4);
	system_matrices.penalty_B.resize(4);
	system_matrices.penalty.resize(4);

	for(unsigned int id = 0 ; id < 4 ; id ++)
	{
		system_matrices.B[id].resize(nbc_M,neqn_M);
		system_matrices.penalty_B[id].resize(neqn_M,neqn_M);
		system_matrices.penalty[id].resize(neqn_M,nbc_M);
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

		system_matrices.B[i].makeCompressed();
		system_matrices.penalty_B[i].makeCompressed();
		system_matrices.penalty[i].makeCompressed();
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
     	system_matrices[i].have_force = true;

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
	
	value = 0;
	//value(2) = x;
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
		case 1:
		case 3:
		{
			value = 0;
			break;
		}
		case 0:
		{

			double bc_normal = 1.0;
			for (unsigned int m = 0 ; m < B.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(B,m); n ; ++n)
                    	if (n.col() == 2)
                    		value(n.row()) += bc_normal * thetaW * n.value()/sqrt(2.0);
                    	
			break;
		}

		case 2:
		{
			double bc_normal = -1.0;
			for (unsigned int m = 0 ; m < B.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(B,m); n ; ++n)
                    	if (n.col() == 2)
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
	Assert(value.size() != 0,ExcNotImplemented());
	Assert(force_vec.size() != 0,ExcNotImplemented());
	const double x = p[0];
	value = 0;

	value(2) = pow(x-0.5,1);
	// value(0) = -M_PI * cos(M_PI*x);
	// value(2) = -sqrt(2) * M_PI * cos(M_PI * x); 
}

