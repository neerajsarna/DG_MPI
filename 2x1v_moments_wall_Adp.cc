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


int main(int argc, char *argv[])
{
      using namespace dealii;
      
      const unsigned int num_threads = atoi(argv[1]);
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, num_threads);

      const int dim = 2;
      const int dim_problem = 1;
      const int poly_degree = 0;

     
     // store the number of equations for a given M (starts at M = 3)
     std::vector<int> neqn_M;
     std::vector<int> nbc_M;
     const double Kn = 0.1;

     const unsigned int num_systems = 12;
     std::vector<int> M(num_systems);
     std::vector<int> M_adjoint(num_systems);
     M[0] = 6;
     M[1] = 8;
     M[2] = 10;
     M[3] = 12;
     M[4] = 14;
     M[5] = 16;
     M[6] = 18;
     M[7] = 20;
     M[8] = 22;
     M[9] = 24;
     M[10] = 26;
     M[11] = 28;
/*     M[1] = 14;
     M[2] = 16;
     M[3] = 40;*/
     
     // does not really matters for the unifrm refinement
     M_adjoint[0] = 8;
     M_adjoint[1] = 10;
     M_adjoint[2] = 12;
     M_adjoint[3] = 14;
     M_adjoint[4] = 16;
     M_adjoint[5] = 18;
     M_adjoint[6] = 20;
     M_adjoint[7] = 22;
     M_adjoint[8] = 24;
     M_adjoint[9] = 26;
     M_adjoint[10] = 28;
     M_adjoint[11] = 30;

     AssertThrow(M.size() == num_systems,ExcNotInitialized());
     AssertThrow(M_adjoint.size() == num_systems,ExcNotInitialized());

     // M_adjoint[1] = 16;
     // M_adjoint[2] = 18;
     // M_adjoint[3] = 20;
     // M_adjoint[4] = 20;
/*     M_adjoint[1] = 8;
     M_adjoint[2] = 10;
     M_adjoint[3] = 11;
*/	
     neqn_M.resize(40);
     nbc_M.resize(40);

     for(unsigned int i = 3 ; i <= 42; i++)
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
      typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

      for(; cell != endc ; cell++)
      		cell->set_refine_flag(RefinementCase<dim>::cut_axis(0));

      triangulation.execute_coarsening_and_refinement();
      set_square_bid(triangulation);

     
      ic_bc<dim> initial_boundary;	
      ic_bc_adjoint<dim> initial_boundary_adjoint;	

      std::string foldername = "2x1v_moments_wall_Adp/";

      const unsigned int max_neqn_primal = system_matrices[system_matrices.size()-1].Ax.rows();
      const unsigned int max_neqn_adjoint = system_matrices_adjoint[system_matrices_adjoint.size()-1].Ax.rows();

       run_problem<dim> Run_Problem(system_matrices,	  // system data
       								system_matrices_error, // system data for computing error
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
	// we first initialise all the matrices
	system_matrices.Ax.resize(neqn_M,neqn_M);
	system_matrices.Ay.resize(neqn_M,neqn_M);
	system_matrices.P.resize(neqn_M,neqn_M);

	system_matrices.B.resize(4);
	system_matrices.penalty_B.resize(4);
	system_matrices.BC_Operator.resize(4);
	system_matrices.penalty.resize(4);

	for(unsigned int id = 0 ; id < 4 ; id ++)
	{
		system_matrices.B[id].resize(nbc_M,neqn_M);
		system_matrices.penalty_B[id].resize(neqn_M,neqn_M);
		system_matrices.penalty[id].resize(neqn_M,nbc_M);
		system_matrices.BC_Operator[id].resize(neqn_M,nbc_M);
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

	filename = "1v_Moments/Bwall/penalty_odd_wall" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.penalty[0],Row_Col_Value);	
	
	filename = "1v_Moments/Bwall/penalty_odd_wall" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.BC_Operator[0],Row_Col_Value);	

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
		system_matrices.BC_Operator[i] = rotator[i].transpose() * system_matrices.BC_Operator[0];

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
		system_matrices.BC_Operator[id].resize(neqn_M,nbc_M);
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
	// value(0) = M_PI * cos(M_PI * x);
	// value(2) = sqrt(2) * M_PI * cos(M_PI * x); 
	value(2) = pow(x-0.5,2);
	//value(2) = x;
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

	// if (t <= 1)
 //        thetaW = exp(-1/(1-pow((t-1),2))) * exp(1);
 //    else
 //        thetaW = 1;

    
	// switch (bc_id)
	// {
	// 	case 1:
	// 	case 3:
	// 	{
	// 		value = 0;
	// 		break;
	// 	}
	// 	case 0:
	// 	{

	// 		double bc_normal = 1.0;
	// 		for (unsigned int m = 0 ; m < B.outerSize() ; m++)
 //                    for (Sparse_Matrix::InnerIterator n(B,m); n ; ++n)
 //                    	if (n.col() == 2)
 //                    		value(n.row()) += bc_normal * thetaW * n.value()/sqrt(2.0);
                    	
	// 		break;
	// 	}

	// 	case 2:
	// 	{
	// 		double bc_normal = -1.0;
	// 		for (unsigned int m = 0 ; m < B.outerSize() ; m++)
 //                    for (Sparse_Matrix::InnerIterator n(B,m); n ; ++n)
 //                    	if (n.col() == 2)
 //                    		value(n.row()) += bc_normal * thetaW * n.value()/sqrt(2.0);
                    	
	// 		break;
	// 	}

	// 	default:
	// 	{
	// 		Assert(1 == 0, ExcMessage("should not have reached"));
	// 		break;
	// 	}
	// }

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

	switch(bc_id)
	{
		case 2:
		{
			value(1) = 1;
			break;
		}
		case 0:
		{
			value(1) = 1;
			break;
		}
	}

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

}


template<int dim>
void 
ic_bc<dim>::exact_solution(const Point<dim> &p,Vector<double> &value,const double &t)
{
	const double shift = 0.5;
	const double x = p[0]-shift; // we need to shift the x coordinate for the reference solution
	const double y = p[1];


	value = 0;

	// M = 6 solution, poisson heat conduction
/*	value(2) = 0.05593013555027514 - 0.00949275700536624*cosh(4.47213595499958*x) + 
    0.13333333333333333*pow(x,2) - 0.2777777777777778*pow(x,4);

    value(4) =  -0.01154700538379252 + 0.008220968718599857*cosh(4.47213595499958*x) - 
    			0.11547005383792518*pow(x,2);
*/

    // M= 40 solution, poisson heat
    value(2) = 0.05673963649111439 - 1.821773279341305e-10*cosh(1.2078840101317283*x) - 
    4.161103543210112e-9*cosh(1.3078453876891523*x) - 
    1.741563843261075e-7*cosh(1.4201039622641722*x) - 
    3.0912702759239734e-6*cosh(1.5481284098563823*x) - 
    0.00003432538090275853*cosh(1.6964567741916035*x) - 
    0.0002201604696294331*cosh(1.8716019970135933*x) - 
    0.0008290727126440684*cosh(2.084723525254245*x) - 
    0.001833451880689609*cosh(2.3577417474531877*x) - 
    0.002595390845546351*cosh(2.7318724706353126*x) - 
    0.0025586976787132197*cosh(3.282073516752524*x) - 
    0.0018323509064832313*cosh(4.162412100887137*x) - 
    0.0008248154566700484*cosh(5.7687839178768066*x) - 
    0.00012511255244483833*cosh(9.551054823962051*x) - 
    9.18904076230852e-9*cosh(28.55918996701305*x) + 0.13333333333333333*pow(x,2) - 
    0.2777777777777778*pow(x,4);

	
   value(4) = -0.01154700538379252 + 
    1.5777019398452543e-10*cosh(1.2078840101317283*x) + 
    3.603621376197396e-9*cosh(1.3078453876891523*x) + 
    1.508238530576551e-7*cosh(1.4201039622641722*x) + 
    2.6771185889138916e-6*cosh(1.5481284098563823*x) + 
    0.0000297266518563661*cosh(1.6964567741916035*x) + 
    0.00019066455960820147*cosh(1.8716019970135933*x) + 
    0.0007179980307342391*cosh(2.084723525254245*x) + 
    0.0015878159052935571*cosh(2.3577417474531877*x) + 
    0.002247674404992714*cosh(2.7318724706353126*x) + 
    0.0022158971903699217*cosh(3.282073516752524*x) + 
    0.0015868624336619229*cosh(4.162412100887137*x) + 
    0.0007143111389103247*cosh(5.7687839178768066*x) + 
    0.0001083506487495429*cosh(9.551054823962051*x) + 
    7.957942736569902e-9*cosh(28.55918996701305*x) - 0.11547005383792518*pow(x,2);
 }
