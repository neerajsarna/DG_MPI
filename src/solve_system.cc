
#include "solve_system.h"


template<int dim>
Solve_System<dim>::Solve_System(system_data &system_mat,
								parallel::distributed::Triangulation<dim> &triangulation,
								const int poly_degree,
								ic_bc_base<dim> *ic_bc,
                std::string &foldername,
                const double min_h)
:
mpi_comm(MPI_COMM_WORLD),
dof_handler(triangulation),
fe_basic(poly_degree),
fe(fe_basic,system_mat.Ax.rows()),
locally_owned_cells(triangulation.n_active_cells()),
initial_boundary(ic_bc),
n_eqn(fe.n_components()),
system_matrices(system_mat),
output_foldername(foldername),
this_mpi_process(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
pout(std::cout,this_mpi_process==0)
// ,
// computing_timer(MPI_COMM_WORLD,
//                     pout,
//                     TimerOutput::summary,
//                     TimerOutput::wall_times)
{
	
	  // we store data of the system ( flux matrices and boundary matrices)

	  dof_handler.distribute_dofs(fe);
	  
	  locally_owned_dofs = dof_handler.locally_owned_dofs();
      DoFTools::extract_locally_relevant_dofs (dof_handler,
                          locally_relevant_dofs) ;
    

      locally_owned_solution.reinit(locally_owned_dofs,mpi_comm);
      locally_relevant_solution.reinit(locally_owned_dofs,locally_relevant_dofs,mpi_comm);


      create_IndexSet_triangulation(); // initialise the locally_owned_cells
      error_per_cell.reinit(locally_owned_cells,mpi_comm); //allocate memory for error per cell

      prescribe_initial_conditions();
      locally_relevant_solution = locally_owned_solution;

      // we need to fix this in case of moments
      max_speed = 1;

      // an approximation to delta_t
      dt = CFL * min_h/max_speed;

}


// prescribe the initial conditions
template<int dim>
void
Solve_System<dim>::prescribe_initial_conditions()
{

	typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
    Vector<double> component(dofs_per_cell);
    Vector<double> value(dofs_per_cell);


    for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
      component[i] = fe.system_to_component_index(i).first;

    for (; cell != endc; ++cell)
    	if (cell->is_locally_owned())
      	{
      		const Point<dim> location = cell->center();

      		for (unsigned int dof = 0 ; dof < dofs_per_cell ; dof++)
      		{
      			value[dof] = initial_boundary->ic(location,component[dof]);
      		}

      		cell->distribute_local_to_global(value,locally_owned_solution);

      	}

  
    locally_relevant_solution = locally_owned_solution;
}


template<int dim>
void 
Solve_System<dim>::run_time_loop()
{
    //TimerOutput::Scope timer_section(computing_timer,"Solving");
    const QGauss<dim-1> face_quadrature(1);

		typename DoFHandler<dim>::cell_iterator neighbor;
		typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();

    const UpdateFlags face_update_flags =  update_normal_vectors;

    FEFaceValues<dim> fe_v_face(fe, face_quadrature, face_update_flags);

		const unsigned int dofs_per_cell = fe.dofs_per_cell;
		const unsigned int dofs_per_component = dofs_per_cell/n_eqn;

		// only implemented a finite volume scheme
		Assert(dofs_per_component == 1,ExcNotImplemented());
		Assert(dofs_per_cell == n_eqn,ExcNotImplemented());
  		
  		Vector<double> component(dofs_per_cell);
    
    	for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
      		component[i] = fe.system_to_component_index(i).first;

      	std::vector<Vector<double>> component_to_system(n_eqn,
                                                          Vector<double>(dofs_per_component));

      for (unsigned int i = 0 ; i < n_eqn ; i ++)
        for (unsigned int j = 0 ; j < dofs_per_component ; j ++)
          	component_to_system[i](j) = fe.component_to_system_index(i,j); 

        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      	std::vector<types::global_dof_index> local_dof_indices_neighbor(dofs_per_cell);

      	Vector<double> neighbor_values(dofs_per_cell);
        std::vector<Vector<double>> g(4);

        for (unsigned int id = 0 ; id < 4 ; id++)
          initial_boundary->bc_inhomo(system_matrices.B[id],id,g[id],0);

      	double t = 0;
      	while (t < t_end)
      	{
      		if (t+dt > t_end)
      			dt = t_end-t;

          // we need to initialise the iterator again and again
          typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();

      		for (; cell != endc; ++cell)
      			if (cell->is_locally_owned())
      			{
            // no hanging nodes check
      				Assert(!cell->has_children(), ExcInternalError());
      				
      				cell->get_dof_indices(local_dof_indices);
      				const double volume = cell->measure();

              Vector<double> flux_contri(dofs_per_cell);

              flux_contri = 0;

          // loop over the faces, we assume no hanging nodes 
      				for(unsigned int face  = 0; face< GeometryInfo<dim>::faces_per_cell; face++ )
      				{
                fe_v_face.reinit(cell,face);
              
      		   		// normal to the face assuming cartesian grid
      				 	Tensor<1,dim> normal_vec = fe_v_face.normal_vector(0);

      					const double nx = normal_vec[0];
      					const double ny = normal_vec[1];

					     // construct An
      				 	Sparse_Matrix An = system_matrices.Ax * nx + system_matrices.Ay * ny;

      					const typename DoFHandler<dim>::face_iterator face_itr = cell->face(face);

                const double face_length = face_itr->measure();

      					if (face_itr->at_boundary())
      					{
                  const unsigned int bc_id = face_itr->boundary_id();
                  // only compute for times greater than zero, already computed for t= 0 before
                    if (system_matrices.bc_inhomo_time && t > 1e-16)
                                initial_boundary->bc_inhomo(system_matrices.B[bc_id],bc_id,
                                                          g[bc_id],t);

                  for (unsigned int m = 0 ; m < An.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An,m); n ; ++n)
                  {
              // 0 because of finite volume
                    unsigned int dof_sol = component_to_system[n.row()](0);

              // the solution id which meets An
                    unsigned int dof_sol_col = local_dof_indices[component_to_system[n.col()](0)];

              // explicit euler update
                    flux_contri(dof_sol) -=  n.value() * dt
                                             * locally_relevant_solution(dof_sol_col) * face_length/volume;

                  }

                  // contribution from penalty_B
                  for (unsigned int m = 0 ; m < system_matrices.penalty_B[bc_id].outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices.penalty_B[bc_id],m); n ; ++n)
                  {
                    unsigned int dof_sol = component_to_system[n.row()](0);
                    unsigned int dof_sol_col = local_dof_indices[component_to_system[n.col()](0)];
                    flux_contri(dof_sol) +=  n.value() * dt
                                             * locally_relevant_solution(dof_sol_col) * face_length/volume;

                  }
                  // contribution from penalty * g
                  for (unsigned int m = 0 ; m < system_matrices.penalty[bc_id].outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices.penalty[bc_id],m); n ; ++n)
                  {
              
                    unsigned int dof_sol = component_to_system[n.row()](0);
                    unsigned int dof_sol_col = local_dof_indices[component_to_system[n.col()](0)];
                    flux_contri(dof_sol) -=  n.value() * dt
                                             * g[bc_id](n.col()) * face_length/volume;

                  }
      					}
      					else
      					{
                  neighbor = cell->neighbor(face);

      						 neighbor->get_dof_indices(local_dof_indices_neighbor);

							     Assert(!neighbor->has_children(), ExcInternalError());

                  for (unsigned int m = 0 ; m < An.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An,m); n ; ++n)
                  {
								// 0 because of finite volume
      							unsigned int dof_sol = component_to_system[n.row()](0);

								// the solution part which meets An
      							unsigned int dof_sol_col = local_dof_indices[component_to_system[n.col()](0)];
      							unsigned int dof_neighbor_col = local_dof_indices_neighbor[component_to_system[n.col()](0)];

					 		// explicit euler update, the two comes from the flux 
      							flux_contri(dof_sol) -=  dt * n.value() * (locally_relevant_solution(dof_sol_col)
      								                       + locally_relevant_solution(dof_neighbor_col)) * face_length/(2 * volume);
              // now add the diffusion
                    flux_contri(dof_sol) -= max_speed * dt * (locally_relevant_solution(dof_sol_col)
                                          - locally_relevant_solution(dof_neighbor_col)) * face_length/(2 * volume);

      						}
       					}


       				} // end of else 

              // update the solution with the fluxes
              for (unsigned int i = 0 ; i < dofs_per_cell ; i ++)
                locally_owned_solution(local_dof_indices[i]) = locally_relevant_solution(local_dof_indices[i])
                                                               + flux_contri(i);

      			} // end of loop over cells

      			locally_relevant_solution = locally_owned_solution;
      			t += dt;
//      			std::cout << "time: " << t << std::endl;

      		}	// end of loop over time

          create_output();

          //computing_timer.print_summary();
          compute_error();

}

template<int dim>
void 
Solve_System<dim>::create_output()
{
  std::string filename = output_foldername + "/result" + std::to_string(this_mpi_process) +  ".txt";

  typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end(), cell = dof_handler.begin_active();

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  Vector<double> solution_value(n_eqn);

  FILE *fp;
  fp = fopen(filename.c_str(),"w+");

  AssertThrow(fp!=NULL,ExcFileNotOpen(filename));
  for(; cell != endc ; cell++)
    if (cell->is_locally_owned())
    {
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int space = 0 ; space < dim ; space ++)
        fprintf(fp, "%f\t",cell->center()(space));

      VectorTools::point_value(dof_handler,locally_owned_solution, cell->center(),solution_value); 

      for (unsigned int i = 0 ; i < n_eqn; i++)
        fprintf(fp, "%f\t",solution_value(i));

      fprintf(fp, "\n");
    }

  fclose(fp);

}

template<int dim>
void 
Solve_System<dim>::compute_error_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                          PerCellErrorScratch &scratch,
                                          PerCellError &data)
{

    if (cell->is_locally_owned())
    {
    data.error_value = 0;
    data.solution_value = 0;

    VectorTools::point_value(dof_handler,locally_owned_solution, 
                             cell->center(),data.solution_value); 
    initial_boundary->exact_solution(cell->center(),data.exact_solution,t_end);

    // overwrite the solution value with the error in this cell
    data.solution_value.sadd(1,-1,data.exact_solution);
    
    for (unsigned int i = 0 ; i < data.solution_value.size(); i++) 
      data.error_value = pow(data.solution_value(i),2) +data.error_value;

    data.cell_index = cell->index();

    data.volume = cell->measure();
    
  }
}

template<int dim>
void
Solve_System<dim>::copy_error_to_global(const PerCellError &data)
{
  error_per_cell(data.cell_index) = sqrt(data.error_value * data.volume);
}


template<int dim>
void
Solve_System<dim>::compute_error()
{
  PerCellError per_cell_error(n_eqn);
  PerCellErrorScratch per_cell_error_scratch;

  const typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                       endc = dof_handler.end();

  WorkStream::run (cell,
                   endc,
                   *this,
                   &Solve_System<dim>::compute_error_per_cell,
                   &Solve_System<dim>::copy_error_to_global,
                   per_cell_error_scratch,
                   per_cell_error);
  
  double errorTot = error_per_cell.l2_norm();
  pout << dof_handler.n_dofs() << "\t" << errorTot << std::endl;

}


template<int dim>
void
Solve_System<dim>::create_IndexSet_triangulation()
{
  typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end(), 
                                                 cell = dof_handler.begin_active();


  for(; cell!=endc;cell++)  
    if(cell->is_locally_owned())
      locally_owned_cells.add_index(cell->index());
    
}
template<int dim>
Solve_System<dim>::~Solve_System()
{
	dof_handler.clear();
}

// explicit initiation to avoid linker error
template class Solve_System<2>;


