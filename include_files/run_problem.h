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
					ic_bc_base<dim> *ic_bc_adjoint);


		void compute_error(const hp::DoFHandler<dim> &dof_handler_primal,
						   const std::vector<system_data> &system_matrices,
						   const std::vector<unsigned int> &n_eqn_primal,
							const std::vector<unsigned int> &n_eqn_adjoint,
						   ic_bc_base<dim> *ic_bc,
						   const std::vector<Vector<double>> &adjoint_solution,
						   const std::vector<Vector<double>> &primal_solution,
						   const std::vector<Point<dim>> &cell_index_primal,
                           const std::vector<Point<dim>> &cell_index_adjoint);

		Vector<double> error_per_cell;

		double t;		// time of the computation


		struct PerCellError
		{
			
			double local_contri;
			unsigned int index;

		};

		struct PerCellErrorScratch
		{
			    PerCellErrorScratch(const hp::FECollection<dim> &fe,
                	        		const hp::QCollection<dim-1> &   quadrature); // initialisation

    			PerCellErrorScratch(const PerCellErrorScratch &scratch);	// copy constructor

    			hp::FEFaceValues<dim> fe_v_face;    			
		};

		void compute_error_per_cell(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
                                      	PerCellErrorScratch &scratch,
                                      	PerCellError &data,
                                        const std::vector<std::vector<Vector<double>>> &g,
                                        const std::vector<system_data> &system_matrices,
                                        const std::vector<unsigned int> &n_eqn_primal,
                                        const std::vector<unsigned int> &n_eqn_adjoint,
                                        const std::vector<Vector<double>> &adjoint_solution,
						   				const std::vector<Vector<double>> &primal_solution,
						   				const std::vector<Point<dim>> &cell_index_primal,
                                        const std::vector<Point<dim>> &cell_index_adjoint);


		void assemble_to_global(const PerCellError &data);

		Sparse_Matrix construct_An_effective(const Sparse_Matrix &An_cell,const Sparse_Matrix &An_neighbor);

		void write_error(const std::string &filename,const Triangulation<dim> &triangulation);
};
