#include "run_problem.h"

template<int dim>
run_problem<dim>::run_problem(std::vector<system_data> &system_mat_primal,	  // system data
                              std::vector<system_data> &system_mat_error_comp,    // system data to compute error
				  			              std::vector<system_data> &system_mat_adjoint, // adjoint data
							  Triangulation<dim> &triangulation, // triangulation
							  const int poly_degree,
							  ic_bc_base<dim> *ic_bc_primal,
					          ic_bc_base<dim> *ic_bc_adjoint)
{
	   Solve_System_SS_adaptive<dim> solve_primal(system_mat_primal,
	  								 			  triangulation,
	  								 			  poly_degree,
	  								 			  ic_bc_primal);

	   Solve_System_SS_adaptive<dim> solve_adjoint(system_mat_adjoint,
	  								 			  triangulation,
	  								 			  poly_degree,
	  								 			  ic_bc_adjoint);


	   const unsigned int refine_cycles = 2;
	   const double tolerance = 5e-2;		// tolerance
	   t = 0;						// we solve for the steady state so set t only initially
	   std::vector<std::vector<Vector<double>>> component_to_system = solve_primal.return_component_to_system(); 
		 error_per_cell.reinit(triangulation.n_active_cells());
	   std::vector<Vector<double>> temp;

	   for(unsigned int cycle = 0 ; cycle < refine_cycles ; cycle++)
	   {
	   		if(cycle == 1)
	   			break;

	   		std::cout << "solving primal: " << std::endl;
	   		solve_primal.run_time_loop(triangulation,cycle,refine_cycles,t,temp);


	   		if(cycle != refine_cycles-1)
	   		{
	   			std::cout << "solving adjoint: " << std::endl;
	   			solve_adjoint.run_time_loop(triangulation,cycle,refine_cycles,t,solve_primal.cellwise_sol);

	   			compute_error(solve_adjoint.dof_handler,
	   						      system_mat_error_comp,
	   					  	  solve_primal.n_eqn,
	   					  	  solve_adjoint.n_eqn,
	   					  	  solve_primal.initial_boundary,
	   					  	  solve_adjoint.cellwise_sol,
                    solve_primal.cellwise_sol,
                    solve_primal.cell_index_center,
                    solve_adjoint.cell_index_center);
	   		}


	   		// solve_primal.allocate_fe_index(cycle + 1);
      //   	solve_primal.distribute_dofs();

      //   	// interpolate onto the current solution
      //   	solve_primal.interpolate_higher_fe_index(solve_primal.cellwise_sol,solve_primal.cell_fe_index,
      //                                 				solve_primal.locally_owned_solution,component_to_system);

      //  		solve_primal.locally_relevant_solution = solve_primal.locally_owned_solution;
	   }


	 	

	std::string filename = "2x3v_moments_HC_Adp/M_4Adj" + std::string("/result0")
                                + "_Kn_" + "0p1" + std::string(".txt");

    solve_adjoint.create_output(filename);

    filename = "2x3v_moments_HC_Adp/error_M3.txt";

    write_error(filename,triangulation);


}

template<int dim>
void 
run_problem<dim>::compute_error(const hp::DoFHandler<dim> &dof_handler,
										 const std::vector<system_data> &system_matrices,
										 const std::vector<unsigned int> &n_eqn_primal,
										 const std::vector<unsigned int> &n_eqn_adjoint,
										 ic_bc_base<dim> *ic_bc,
                     const std::vector<Vector<double>> &adjoint_solution,
                     const std::vector<Vector<double>> &primal_solution,
                     const std::vector<Point<dim>> &cell_index_primal,
                     const std::vector<Point<dim>> &cell_index_adjoint)
{
	const typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
													 	                               endc = dof_handler.end();												   
	error_per_cell = 0;

	QGauss<dim-1> face_quadrature_basic(1);
	hp::QCollection<dim-1> face_quadrature;

	for (unsigned int i = 0 ; i < dof_handler.get_fe().size() ; i++)
		face_quadrature.push_back(face_quadrature_basic);

	PerCellError per_cell_error;
	PerCellErrorScratch per_cell_error_scratch(dof_handler.get_fe(),face_quadrature);

    std::vector<std::vector<Vector<double>>> g(n_eqn_adjoint.size());

    for (unsigned int i = 0 ; i < n_eqn_adjoint.size() ; i ++)
    {
          g[i].resize(4);
          for (unsigned int id = 0 ; id < 4 ; id++)
            ic_bc->bc_inhomo(system_matrices[i].B[id],id,g[i][id],t);
     }

	WorkStream::run(cell,
		endc,
		std::bind(&run_problem<dim>::compute_error_per_cell,
			this,
			std::placeholders::_1,
			std::placeholders::_2,
			std::placeholders::_3,
			std::cref(g),
			std::cref(system_matrices),
			std::cref(n_eqn_primal),
			std::cref(n_eqn_adjoint),
			std::cref(adjoint_solution),
			std::cref(primal_solution),
      std::cref(cell_index_primal),
      std::cref(cell_index_adjoint)),
		std::bind(&run_problem<dim>::assemble_to_global,
			this,
			std::placeholders::_1),
		per_cell_error_scratch,
		per_cell_error);
}


template<int dim>
void 
run_problem<dim>::compute_error_per_cell(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
                                      	PerCellErrorScratch &scratch,
                                      	PerCellError &data,
                                        const std::vector<std::vector<Vector<double>>> &g,
                                        const std::vector<system_data> &system_matrices,
                                        const std::vector<unsigned int> &n_eqn_primal,
										                    const std::vector<unsigned int> &n_eqn_adjoint,
                                        const std::vector<Vector<double>> &adjoint_solution,
                                        const std::vector<Vector<double>> &primal_solution,
                                        const std::vector<Point<dim>> &cell_index_primal,
                                        const std::vector<Point<dim>> &cell_index_adjoint)
 {

              // operations to avoid data races, computations are in steady state so we do not need to
 			        // change the value of temp_g.
              Assert(cell_index_adjoint[cell->index()].distance(cell_index_primal[cell->index()])<1e-15,ExcMessage("indexing changed"));
              std::vector<std::vector<Vector<double>>> temp_g = g;
              const unsigned int this_fe_index = cell->active_fe_index();

              // no hanging nodes check
              Assert(!cell->has_children(), ExcInternalError());
              const unsigned int index = cell->index();
              const double volume = cell->measure();

              data.local_contri = 0;
              data.index = index;


              // contribution from collisions P
              for (unsigned int m = 0 ; m < system_matrices[this_fe_index].P.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices[this_fe_index].P,m); n ; ++n)
                    	if(n.col() < n_eqn_primal[this_fe_index])
                  {
              // 0 because of finite volume
                    const double adjoint_value = adjoint_solution[index](n.row());
                    const double solution_value = primal_solution[index](n.col());

              // explicit euler update, 
                    data.local_contri +=  n.value() * adjoint_value * solution_value * volume;

                  }


              // loop over the faces, we assume no hanging nodes 
              for(unsigned int face  = 0; face < GeometryInfo<dim>::faces_per_cell; face++ )
              {
                scratch.fe_v_face.reinit(cell,face);
                const FEFaceValues<dim> &fe_v_face_temp = scratch.fe_v_face.get_present_fe_values();
                
                // normal to the face assuming cartesian grid
                Tensor<1,dim> normal_vec = fe_v_face_temp.normal_vector(0);

                const double nx = normal_vec[0];
                const double ny = normal_vec[1];

                const typename Triangulation<dim>::face_iterator face_itr = cell->face(face);
                const double face_length = face_itr->measure();

                  // construct An of the current cell
                  Sparse_Matrix An_cell = system_matrices[this_fe_index].Ax * nx
                                          + system_matrices[this_fe_index].Ay * ny;

                if (face_itr->at_boundary())
                {

                  const unsigned int bc_id = face_itr->boundary_id();

                  for (unsigned int m = 0 ; m < An_cell.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An_cell,m); n ; ++n)
                    	if(n.col() < n_eqn_primal[this_fe_index])
                  {
                      // 0 because of finite volume
              		 // 0 because of finite volume
                    const double adjoint_value = adjoint_solution[index](n.row());
                    const double solution_value = primal_solution[index](n.col());

              	    // explicit euler update, 
                    data.local_contri -=  n.value() * adjoint_value * solution_value * face_length;

                  }

                  // contribution from penalty_B
                  for (unsigned int m = 0 ; m < system_matrices[this_fe_index].penalty_B[bc_id].outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices[this_fe_index].penalty_B[bc_id],m); n ; ++n)
                    	if(n.col() < n_eqn_primal[this_fe_index])
                  {
                    // 0 because of finite volume
                    const double adjoint_value = adjoint_solution[index](n.row());
                    const double solution_value = primal_solution[index](n.col());

              		// explicit euler update
                    data.local_contri +=  n.value() * adjoint_value * solution_value * face_length;

                  }

                  // contribution from penalty * g
                  for (unsigned int m = 0 ; m < system_matrices[this_fe_index].penalty[bc_id].outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices[this_fe_index].penalty[bc_id],m); n ; ++n)
                  {
                    // 0 because of finite volume
                    const double adjoint_value = adjoint_solution[index](n.row());

                    data.local_contri -=  n.value()
                                          * temp_g[this_fe_index][bc_id](n.col()) 
                                          * adjoint_value * face_length;

                  }
                }
                else
                {

                   typename hp::DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face);
                   const unsigned int neighbor_fe_index = neighbor->active_fe_index();
                   const unsigned int index_neighbor = neighbor->index();

                   Sparse_Matrix An_neighbor = system_matrices[neighbor_fe_index].Ax * nx
                                              + system_matrices[neighbor_fe_index].Ay * ny;

                   Sparse_Matrix An_effective = construct_An_effective(An_cell,An_neighbor);

                   Sparse_Matrix Amod;

                   if (fabs(ny) < 1e-16)
                        Amod = system_matrices[std::max(this_fe_index,neighbor_fe_index)].Ax_mod;
                   else
                        Amod = system_matrices[std::max(this_fe_index,neighbor_fe_index)].Ay_mod;

                   Assert(!neighbor->has_children(), ExcInternalError());

                   // contribution from the present cell
                  for (unsigned int m = 0 ; m < An_cell.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An_cell,m); n ; ++n)
                    	if(n.col() < n_eqn_primal[this_fe_index])
                  {
                                  // 0 because of finite volume
                    const double adjoint_value = adjoint_solution[index](n.row());
                    const double solution_value = primal_solution[index](n.col());

                    data.local_contri -=  n.value() * adjoint_value * solution_value
                                                   * face_length/(2);

                  }

                  // we add the diffusion now
                  for (unsigned int m = 0 ; m < Amod.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(Amod,m); n ; ++n)
                      if(n.row() < n_eqn_adjoint[this_fe_index] && n.col() < n_eqn_primal[this_fe_index])
                  {
                    const double adjoint_value = adjoint_solution[index](n.row());
                    const double solution_value = primal_solution[index](n.col());

                    data.local_contri -=  n.value() * adjoint_value * solution_value
                                                   * face_length/(2);

                  }

                  // contribution from the neighboring cell
                  for (unsigned int m = 0 ; m < An_effective.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An_effective,m); n ; ++n)
                    	if(n.col() < n_eqn_primal[neighbor_fe_index])
                  {
                    const double adjoint_value = adjoint_solution[index](n.row());
                    const double neighbor_value = primal_solution[index_neighbor](n.col());

                    data.local_contri -=  n.value() * adjoint_value * neighbor_value
                                                   * face_length/(2);

                  }

                  // diffusion with the neighbouring cell
                  for (unsigned int m = 0 ; m < Amod.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(Amod,m); n ; ++n)
                      if(n.row() < n_eqn_adjoint[this_fe_index] && n.col() < n_eqn_primal[neighbor_fe_index])
                  {
                    const double adjoint_value = adjoint_solution[index](n.row());
                    const double neighbor_value = primal_solution[index_neighbor](n.col());
                    data.local_contri -=  n.value() * (-neighbor_value) * adjoint_value
                                                   * face_length/(2);

                  }

                } //end of else


              } 
              //end of loop over the faces


 }


template<int dim>
void
run_problem<dim>::assemble_to_global(const PerCellError &data)
 {
        error_per_cell(data.index) = data.local_contri;
 }

 template<int dim>
run_problem<dim>::PerCellErrorScratch::PerCellErrorScratch(const hp::FECollection<dim> &fe,
                                                                 const hp::QCollection<dim-1> &quadrature)
:
fe_v_face(fe,quadrature,update_normal_vectors)
{;}

template<int dim>
run_problem<dim>::PerCellErrorScratch::PerCellErrorScratch(const PerCellErrorScratch &scratch)
:
fe_v_face(scratch.fe_v_face.get_fe_collection(),
          scratch.fe_v_face.get_quadrature_collection(),
          update_normal_vectors)
{;}

template<int dim>
Sparse_Matrix
run_problem<dim>::construct_An_effective(const Sparse_Matrix &An_cell,const Sparse_Matrix &An_neighbor)
{
  Assert(An_cell.rows() != 0 ,ExcNotInitialized());
  Assert(An_cell.rows() == An_cell.cols(),ExcNotInitialized());

  Assert(An_neighbor.rows() == An_neighbor.cols(),ExcNotInitialized());
  Assert(An_neighbor.rows() != 0 ,ExcNotInitialized());

  const unsigned int rows_cell = An_cell.rows();
  const unsigned int rows_neighbor = An_neighbor.rows();

  // first we copy the biggest one
  Sparse_Matrix result = rows_cell >= rows_neighbor ? An_cell.block(0,0,rows_cell,rows_neighbor)
                                                    : An_neighbor.block(0,0,rows_cell,rows_neighbor);

  return(result);

}

template<int dim>
void 
run_problem<dim>::write_error(const std::string &filename,const Triangulation<dim> &triangulation)
{
	FILE *fp;
	fp = fopen(filename.c_str(),"w+");

	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

	for(; cell != endc ; cell++)
	{
		 for (unsigned int space = 0 ; space < dim ; space ++)
        		fprintf(fp, "%f\t",cell->center()(space));

         fprintf(fp, "%0.16f\n",error_per_cell(cell->index()));
	}

	fclose(fp);
}


template class run_problem<2>;
