
#include "solve_system_SS.h"


template<int dim>
Solve_System_SS<dim>::Solve_System_SS(system_data &system_mat,
                parallel::distributed::Triangulation<dim> &triangulation,
								const int poly_degree,
								ic_bc_base<dim> *ic_bc,
                std::string &foldername)
:
mpi_comm(MPI_COMM_WORLD),
dof_handler(triangulation),
fe_basic(poly_degree),
fe(fe_basic,system_mat.Ax.rows()),
initial_boundary(ic_bc),
system_matrices(system_mat),
output_foldername(foldername),
this_mpi_process(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
pout(std::cout,this_mpi_process==0),
computing_timer(MPI_COMM_WORLD,
                    pout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
{
      n_eqn = fe.n_components();
      distribute_dofs();

      prescribe_initial_conditions();
      locally_owned_solution.compress(VectorOperation::insert);
      locally_relevant_solution = locally_owned_solution;

      
      max_speed = compute_max_speed();

}

template<int dim>
void
Solve_System_SS<dim>::run_time_loop(parallel::distributed::Triangulation<dim> &triangulation)
{
      computing_timer.enter_subsection("solving");
          
      // do everything here 
      // start of refinement 
      const int refine_cycles = 1;

      const int times_refine = 1;

      int total_steps = 0;

      double t = 0;

      for (int cycle = 0 ; cycle < refine_cycles ; cycle++)
      {

        // an approximation to delta_t
        double min_length = min_h(triangulation);
        Utilities::MPI::min(min_length,mpi_comm);

        dt = CFL * min_length/(dim * max_speed);

        pout << "refinement cycle: " << cycle << std::endl;
        pout << "total cells: " << triangulation.n_global_active_cells() << std::endl;
        pout << "total dofs: " << dof_handler.n_dofs() << std::endl;

        typedef FilteredIterator<typename DoFHandler<dim>::active_cell_iterator> CellFilter;

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


          QGauss<dim-1> face_quadrature(1);
          PerCellAssemble per_cell_assemble;
          PerCellAssembleScratch per_cell_assemble_scratch(fe,face_quadrature);

          int step_count = 0;
          residual_ss = 100;

          std::vector<Vector<double>> g(4);

          for (unsigned int id = 0 ; id < 4 ; id++)
            initial_boundary->bc_inhomo(system_matrices.B[id],id,g[id],t);

    while (step_count < 100 || residual_ss > 1e-5 ) // we atleast run till t_end || residual_ss > 1e-8
    {
      WorkStream::run ( CellFilter(IteratorFilters::LocallyOwnedCell(),
        dof_handler.begin_active()),
      CellFilter(IteratorFilters::LocallyOwnedCell(), 
        dof_handler.end()),
      std::bind(&Solve_System_SS<dim>::assemble_per_cell,
        this,
        std::placeholders::_1,
        std::placeholders::_2,
        std::placeholders::_3,
        std::cref(component),
        std::cref(component_to_system),
        std::cref(t),
        std::cref(g)),
      std::bind(&Solve_System_SS<dim>::assemble_to_global,
        this,
        std::placeholders::_1,
        std::cref(component)),
      per_cell_assemble_scratch,
      per_cell_assemble);

      locally_owned_solution.compress(VectorOperation::insert); // synchornising 
      locally_owned_residual.compress(VectorOperation::insert);

      residual_ss = locally_owned_residual.l2_norm();
      locally_relevant_solution = locally_owned_solution;

      t += dt;
      step_count++;

      if (step_count%500 == 0)
        pout<< "residual: " << residual_ss << "time: " << t << std::endl;

    }

    total_steps += step_count;

  
    if (cycle != refine_cycles-1)  // no need for the last computation
    {

    //we now need to interpolate onto the new mesh
    locally_relevant_solution_temp = locally_owned_solution; // initialise the temporary

    for (unsigned int i = 0 ; i < times_refine ; i++)
    {
     parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> soltrans(dof_handler);


     typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end(), 
     cell = dof_handler.begin_active();  

     for(; cell != endc ; cell++)
      if (cell->is_locally_owned())
        cell->set_refine_flag();

      soltrans.prepare_for_coarsening_and_refinement(locally_relevant_solution_temp);
      triangulation.refine_global();
      distribute_dofs();        // reinitialise all the vectors
      soltrans.interpolate(locally_owned_solution); // interpolate to the new vector

      locally_owned_solution.compress(VectorOperation::insert);
      locally_relevant_solution = locally_owned_solution; // update the locally relevant solution
    }

  }
}

  pout << "Steps taken: " << total_steps << "final residual: "<< residual_ss << std::endl;


}

template<int dim>
void 
Solve_System_SS<dim>::distribute_dofs()
{
      dof_handler.distribute_dofs(fe);
    
      locally_owned_dofs = dof_handler.locally_owned_dofs();
      DoFTools::extract_locally_relevant_dofs (dof_handler,
                          locally_relevant_dofs) ;
    

      locally_owned_solution.reinit(locally_owned_dofs,mpi_comm);
      locally_owned_residual.reinit(locally_owned_dofs,mpi_comm);
      locally_relevant_solution.reinit(locally_owned_dofs,locally_relevant_dofs,mpi_comm);

      // allocate the memory error 
      IndexSet locally_owned_cells(dof_handler.get_triangulation().n_active_cells()); // helpful while computing error
      create_IndexSet_triangulation(locally_owned_cells); 
      error_per_cell.reinit(locally_owned_cells,mpi_comm);
}


// prescribe the initial conditions
template<int dim>
void
Solve_System_SS<dim>::prescribe_initial_conditions()
{

  typedef FilteredIterator<typename DoFHandler<dim>::active_cell_iterator> CellFilter;

  PerCellIC percellIC;
  PerCellICScratch percellICScratch;

  Vector<double> component(fe.dofs_per_cell);
  
  for (unsigned int i = 0 ; i < fe.dofs_per_cell ; i++)
      component[i] = fe.system_to_component_index(i).first;

  WorkStream::run (CellFilter(IteratorFilters::LocallyOwnedCell(),
                    dof_handler.begin_active()),
                   CellFilter(IteratorFilters::LocallyOwnedCell(), dof_handler.end()),
                   std::bind(&Solve_System_SS<dim>::compute_ic_per_cell,
                            this,
                            std::placeholders::_1,
                            std::placeholders::_2,
                            std::placeholders::_3,
                            std::cref(component)),
                   std::bind(&Solve_System_SS<dim>::copy_ic_to_global,
                            this,
                            std::placeholders::_1),
                   percellICScratch,
                   percellIC);


}

template<int dim>
void
Solve_System_SS<dim>::compute_ic_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                        PerCellICScratch &scratch,
                                        PerCellIC &data,
                                        const Vector<double> &component)
{
  
          data.dofs_per_cell = fe.dofs_per_cell;
          data.local_dof_indices.resize(data.dofs_per_cell);
          data.local_contri.reinit(data.dofs_per_cell);

          const Point<dim> location = cell->center();

          for (unsigned int dof = 0 ; dof < data.dofs_per_cell ; dof++)
            data.local_contri[dof] = initial_boundary->ic(location,component[dof]);
          
          cell->get_dof_indices(data.local_dof_indices);

  
}

template<int dim>
void
Solve_System_SS<dim>::copy_ic_to_global(const PerCellIC &data)
{
    for(unsigned int dof = 0 ; dof< data.dofs_per_cell ; dof++)
      locally_owned_solution(data.local_dof_indices[dof]) = data.local_contri(dof);
}


template<int dim>
void 
Solve_System_SS<dim>::assemble_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                      PerCellAssembleScratch &scratch,
                                      PerCellAssemble &data,
                                      const Vector<double> &component,
                                      const std::vector<Vector<double>> &component_to_system,
                                      const double &t,
                                      const std::vector<Vector<double>> &g)
 {

              // operations to avoid data races
              std::vector<Vector<double>> temp_g(4);
              for (unsigned int id = 0 ; id < 4 ; id ++)
                temp_g[id] = g[id];

              // no hanging nodes check
              Assert(!cell->has_children(), ExcInternalError());
              
              data.dofs_per_cell = fe.dofs_per_cell;
              data.local_dof_indices.resize(data.dofs_per_cell);
              cell->get_dof_indices(data.local_dof_indices);
              const double volume = cell->measure();

              data.local_contri.reinit(data.dofs_per_cell);
              data.local_contri = 0;

              // contribution from collisions P
              for (unsigned int m = 0 ; m < system_matrices.P.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices.P,m); n ; ++n)
                  {
              // 0 because of finite volume
                    unsigned int dof_sol = component_to_system[n.row()](0);

              // the solution id which meets An
                    unsigned int dof_sol_col = data.local_dof_indices[component_to_system[n.col()](0)];

              // explicit euler update, 
                    data.local_contri(dof_sol) +=  n.value() * dt
                                             * locally_relevant_solution(dof_sol_col);

                  }

              // loop over the faces, we assume no hanging nodes 
              for(unsigned int face  = 0; face< GeometryInfo<dim>::faces_per_cell; face++ )
              {
                scratch.fe_v_face.reinit(cell,face);
              
                // normal to the face assuming cartesian grid
                Tensor<1,dim> normal_vec = scratch.fe_v_face.normal_vector(0);

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
                                                          temp_g[bc_id],t);

                  for (unsigned int m = 0 ; m < An.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An,m); n ; ++n)
                  {
                      // 0 because of finite volume
                      unsigned int dof_sol = component_to_system[n.row()](0);

                     // the solution id which meets An
                      unsigned int dof_sol_col = data.local_dof_indices[component_to_system[n.col()](0)];

                    // explicit euler update
                      data.local_contri(dof_sol) -=  n.value() * dt
                                             * locally_relevant_solution(dof_sol_col) * face_length/volume;

                  }

                  // contribution from penalty_B
                  for (unsigned int m = 0 ; m < system_matrices.penalty_B[bc_id].outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices.penalty_B[bc_id],m); n ; ++n)
                  {
                    unsigned int dof_sol = component_to_system[n.row()](0);
                    unsigned int dof_sol_col = data.local_dof_indices[component_to_system[n.col()](0)];
                    data.local_contri(dof_sol) +=  n.value() * dt
                                             * locally_relevant_solution(dof_sol_col) * face_length/volume;

                  }
                  // contribution from penalty * g
                  for (unsigned int m = 0 ; m < system_matrices.penalty[bc_id].outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices.penalty[bc_id],m); n ; ++n)
                  {
              
                    unsigned int dof_sol = component_to_system[n.row()](0);
                    unsigned int dof_sol_col = data.local_dof_indices[component_to_system[n.col()](0)];
                    data.local_contri(dof_sol) -=  n.value() * dt
                                             * temp_g[bc_id](n.col()) * face_length/volume;

                  }
                }
                else
                {
                   typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face);
                   std::vector<types::global_dof_index> local_dof_indices_neighbor(data.dofs_per_cell);

                   neighbor->get_dof_indices(local_dof_indices_neighbor);

                   Assert(!neighbor->has_children(), ExcInternalError());

                  for (unsigned int m = 0 ; m < An.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An,m); n ; ++n)
                  {
                // 0 because of finite volume
                    unsigned int dof_sol = component_to_system[n.row()](0);

                // the solution part which meets An
                    unsigned int dof_sol_col = data.local_dof_indices[component_to_system[n.col()](0)];
                    unsigned int dof_neighbor_col = local_dof_indices_neighbor[component_to_system[n.col()](0)];

              // explicit euler update, the two comes from the flux 
                    data.local_contri(dof_sol) -=  dt * n.value() * (locally_relevant_solution(dof_sol_col)
                                             + locally_relevant_solution(dof_neighbor_col)) * face_length/(2 * volume);

                  }

                  // we add the diffusion now
                   for (unsigned int id = 0 ; id < n_eqn ; id++)
                   {
                       unsigned int dof_sol = component_to_system[id](0);                    
                       unsigned int dof_sol_col = data.local_dof_indices[component_to_system[id](0)];
                       unsigned int dof_neighbor_col = local_dof_indices_neighbor[component_to_system[id](0)];

                       data.local_contri(dof_sol) -= max_speed * dt * (locally_relevant_solution(dof_sol_col)
                                           - locally_relevant_solution(dof_neighbor_col)) * face_length/(2 * volume);
                   }
                } //end of else


              } 
              //end of loop over the faces


 }

template<int dim>
void
Solve_System_SS<dim>::assemble_to_global(const PerCellAssemble &data,const Vector<double> &component)
 {
        for (unsigned int i = 0 ; i < data.dofs_per_cell ; i ++)
        {
              locally_owned_solution(data.local_dof_indices[i]) = locally_relevant_solution(data.local_dof_indices[i])
                                                               + data.local_contri(i);

              if (component[i] <= 7) // only count till tensor degree 2
                locally_owned_residual(data.local_dof_indices[i]) = fabs(data.local_contri(i))/dt;          
        }

 }

template<int dim>
void 
Solve_System_SS<dim>::create_output(const std::string &filename)
{

  typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end(), cell = dof_handler.begin_active();

  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
  const unsigned int dofs_per_component = dofs_per_cell/n_eqn;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  
  std::vector<Vector<double>> component_to_system(n_eqn,Vector<double>(dofs_per_component));

  for (unsigned int i = 0 ; i < n_eqn ; i ++)
      for (unsigned int j = 0 ; j < dofs_per_component ; j ++)
            component_to_system[i](j) = fe.component_to_system_index(i,j); 

  const unsigned int to_print = 13;

  FILE *fp;
  fp = fopen(filename.c_str(),"w+");

  AssertThrow(fp!=NULL,ExcFileNotOpen(filename));
  for(; cell != endc ; cell++)
    if (cell->is_locally_owned())
    {
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int space = 0 ; space < dim ; space ++)
        fprintf(fp, "%f\t",cell->center()(space));

      for (unsigned int i = 0 ; i < to_print ; i++) 
      {
        const double sol_value =  locally_owned_solution(local_dof_indices[component_to_system[i](0)]);
        fprintf(fp, "%f\t",sol_value);
      }

      fprintf(fp, "\n");
    }

  fclose(fp);

}

template<int dim>
void 
Solve_System_SS<dim>::compute_error_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                             PerCellErrorScratch &scratch,
                                             PerCellError &data)
{

    data.error_value = 0;
    data.solution_value = 0;

    VectorTools::point_value(cell->get_dof_handler(),locally_owned_solution, 
                             cell->center(),data.solution_value); 
    initial_boundary->exact_solution(cell->center(),data.exact_solution,t_end);

    // overwrite the solution value with the error in this cell
    data.solution_value.sadd(1,-1,data.exact_solution);
    
    for (unsigned int i = 0 ; i < data.solution_value.size(); i++) 
      data.error_value = pow(data.solution_value(i),2) +data.error_value;

    data.cell_index = cell->index();

    data.volume = cell->measure();
    
  
}

template<int dim>
void
Solve_System_SS<dim>::copy_error_to_global(const PerCellError &data)
{
  error_per_cell(data.cell_index) = sqrt(data.error_value * data.volume);
}


template<int dim>
void
Solve_System_SS<dim>::compute_error()
{
  PerCellError per_cell_error(n_eqn);
  PerCellErrorScratch per_cell_error_scratch;

  typedef FilteredIterator<typename DoFHandler<dim>::active_cell_iterator> CellFilter;

  WorkStream::run ( CellFilter(IteratorFilters::LocallyOwnedCell(),
                    dof_handler.begin_active()),
                    CellFilter(IteratorFilters::LocallyOwnedCell(), 
                    dof_handler.end()),
                    std::bind(&Solve_System_SS<dim>::compute_error_per_cell,
                            this,
                            std::placeholders::_1,
                            std::placeholders::_2,
                            std::placeholders::_3),
                    std::bind(&Solve_System_SS<dim>::copy_error_to_global,
                            this,
                            std::placeholders::_1),
                   per_cell_error_scratch,
                   per_cell_error);
  
  double errorTot = error_per_cell.l2_norm();
  pout << dof_handler.n_dofs() << "\t" << errorTot << std::endl;

}


template<int dim>
void
Solve_System_SS<dim>::create_IndexSet_triangulation(IndexSet &locally_owned_cells)
{
  typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end(), 
                                                 cell = dof_handler.begin_active();

  for(; cell!=endc;cell++)  
    if(cell->is_locally_owned())
      locally_owned_cells.add_index(cell->index());
    
}


template<int dim>
double
Solve_System_SS<dim>::compute_max_speed()
{
    EigenSolver<MatrixXd> ES(system_matrices.Ax);
    VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal().cwiseAbs();

    return(vals.maxCoeff());
}


template<int dim>
double
Solve_System_SS<dim>::min_h(parallel::distributed::Triangulation<dim> &triangulation)
{

  typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(),
                                                                           endc = triangulation.end();


  double min_length = 10000;

  for(; cell != endc ; cell++)
    if (cell->is_locally_owned())
      for(unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
        if (cell->face(face)->measure() < min_length)
          min_length = cell->face(face)->measure();


  return(min_length);
}

template<int dim>
Solve_System_SS<dim>::PerCellAssembleScratch::PerCellAssembleScratch(const FiniteElement<dim> &fe,
  const Quadrature<dim-1> &   quadrature)
:
fe_v_face(fe,quadrature,update_normal_vectors)
{;}

template<int dim>
Solve_System_SS<dim>::PerCellAssembleScratch::PerCellAssembleScratch(const PerCellAssembleScratch &scratch)
:
fe_v_face(scratch.fe_v_face.get_fe(),
  scratch.fe_v_face.get_quadrature(),
  update_normal_vectors)
{;}

template<int dim>
Solve_System_SS<dim>::~Solve_System_SS()
{
  dof_handler.clear();
}
// explicit initiation to avoid linker error
template class Solve_System_SS<2>;


