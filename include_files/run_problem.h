#include "solve_system_SS_adaptive.h"

template<int dim>
class
run_problem
{
	public:
		run_problem(std::vector<system_data> &system_mat,	  // system data
					std::vector<system_data> &system_mat_error,	  // system data
				    std::vector<system_data> &system_mat_adj, // adjoint data
					Triangulation<dim> &triangulation, // triangulation
					const int poly_degree,
					ic_bc_base<dim> *ic_bc_primal,
					ic_bc_base<dim> *ic_bc_adjoint,
					const std::string &foldername,
          const unsigned int &max_equations_primal,
          const unsigned int &max_equations_adjoint);


    DoFHandler<dim> dummy_dof_handler_grid;      // a dummy dof handler for grid refinement
    DoFHandler<dim> dummy_dof_handler_velocity; // a dummy dof handler for the velocity space refinement
    
    FESystem<dim> dummy_fe_grid;
    FESystem<dim> dummy_fe_velocity;

		Vector<double> error_per_cell_velocity;
    Vector<double> error_per_cell_grid;
    Vector<double> error_per_cell_grid_target;
    Vector<double> store_user_index;          // dummy vector to transfer user index during grid refinement

		 double t;		// time of the computation

		 void write_error(const std::string &filename,const Triangulation<dim> &triangulation,const Vector<double> &error_per_cell);
     void write_grid(const std::string &filename,const Triangulation<dim> &triangulation);

    ConvergenceTable convergence_table;
    void develop_convergence_table(const double &error_primal,
                                   const double &error_adjoint,
                                   const double &min_h,
                                   const unsigned int &num_dofs);
    void print_convergence_table();

    struct PerCellErrorScratch
    {
         PerCellErrorScratch(const FiniteElement<dim> &fe,
                               const QGauss<dim> &   quadrature_int,
                               const QGauss<dim-1> &   quadrature_face); // initialisation

         PerCellErrorScratch(const PerCellErrorScratch &scratch);  // copy constructor
          
         FEValues<dim> fe_v;      
         FEFaceValues<dim> fe_v_face;          
         FESubfaceValues<dim> fe_v_subface;          
                
    };

    struct PerCellError
    {
      
     double local_contri;
     unsigned int active_index;

    };

    void compute_error_velocity(const Vector<double> &primal_solution,
                                const DoFHandler<dim> &dof_handler_primal,
                                Vector<double> &adjoint_solution,
                                const DoFHandler<dim> &dof_handler_adjoint,
                                const std::vector<system_data> &system_matrices,
                                ic_bc_base<dim> *ic_bc_primal);

    void compute_error_grid(const Vector<double> &primal_solution,
                                const DoFHandler<dim> &dof_handler_primal,
                                Vector<double> &adjoint_solution,
                                const DoFHandler<dim> &dof_handler_adjoint,
                                const std::vector<system_data> &system_matrices,
                                ic_bc_base<dim> *ic_bc_primal,
                                const Triangulation<dim> &triangulation);

    typename DoFHandler<dim>::cell_iterator return_child_refined_neighbor(const typename DoFHandler<dim>::cell_iterator &neighbor,
                                                                         const typename DoFHandler<dim>::active_cell_iterator &cell);

    void compute_error_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                PerCellErrorScratch &scratch,
                                        PerCellError &data,
                                        const std::vector<Vector<double>> &g,
                                        const std::vector<system_data> &system_matrices,
                                        const Vector<double> &adjoint_solution,
                                        const Vector<double> &primal_solution,
                                        const std::vector<Vector<double>> &component_to_system);


    void error_face(const typename DoFHandler<dim>::cell_iterator &neighbor,  // iterator of the neighbor
                                const typename DoFHandler<dim>::cell_iterator &cell,  // iterator of the cell
                                const Vector<double> &primal_solution,       // primal solution
                                const Vector<double> &adjoint_solution,                         // the adjoint solution
                                const std::vector<system_data> &system_matrices,          // matrices
                                const std::vector<Vector<double>> &component_to_system,         // component to system index of the adjoint
                                const FEValuesBase<dim> &fe_v_face,                       
                                double &result);

     void get_value_at_quad(const std::vector<types::global_dof_index> &local_dof_indices,
                                    const std::vector<Vector<double>> &component_to_system,
                                    Vector<double> &adjoint_value,  // value to be filled
                                    const unsigned int &dofs_per_component,
                                    const Vector<double> &adjoint_solution,
                                    const FEValuesBase<dim> &fe_v,
                                    const unsigned int &q);


        // // x^T A y
    double xAy(const Vector<double> &x,
              const Sparse_Matrix &A,
              const Vector<double> &y);

    void assemble_to_global(const PerCellError &data,Vector<double> &input);

    void compute_error(const Vector<double> &primal_solution,
                                const Vector<double> &adjoint_solution,
                                const DoFHandler<dim> &dof_handler_adjoint,
                                const std::vector<system_data> &system_matrices,
                                ic_bc_base<dim> *ic_bc_primal,
                                Vector<double> &error_vector,
                                const unsigned int &quad_points);

    void fill_index_vector_from_user_index();
    void fill_user_index_from_index_vector();

    void update_index_vector(Triangulation<dim> &triangulation,
                             const unsigned int &refinement_type);

    void update_grid_refine_flags(Triangulation<dim> &triangulation,
                              const unsigned int &refinement_type);

    void perform_grid_refinement_and_sol_transfer(Triangulation<dim> &triangulation,
                                                  Solve_System_SS_adaptive<dim> &solve_primal,
                                                  Solve_System_SS_adaptive<dim> &solve_adjoint);


    void write_fe_index(const std::string &filename,
                                const DoFHandler<dim> &dof_handler,
                                const std::vector<unsigned int> &n_eqn);


    void print_fe_index(const DoFHandler<dim> &dof_handler,
                        const std::vector<unsigned int> &n_eqn);

    Vector<double> modulus_Vec(const Vector<double> &vecIn);

    void compute_error_in_target(const Triangulation<dim> &triangulation,
                                     ic_bc_base<dim> *ic_bc_primal,
                                     ic_bc_base<dim> *ic_bc_adjoint,
                                     const DoFHandler<dim> &dof_handler_primal,
                                     const Vector<double> &primal_solution,
                                     const std::vector<system_data> &system_matrices_adjoint,
                                     const double &t);

    double return_face_length(const typename DoFHandler<dim>::face_iterator &face_itr);
};
