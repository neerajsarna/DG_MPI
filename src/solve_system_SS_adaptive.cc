#include "solve_system_SS_adaptive.h"


template<int dim>
Solve_System_SS_adaptive<dim>::Solve_System_SS_adaptive(std::vector<system_data> &system_mat,
                                                        Triangulation<dim> &triangulation,
							                                         	const int poly_degree,
								                                        ic_bc_base<dim> *ic_bc)
:
mpi_comm(MPI_COMM_WORLD),
dof_handler(triangulation),
fe_basic(poly_degree),
initial_boundary(ic_bc),
system_matrices(system_mat),
this_mpi_process(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
pout(std::cout,this_mpi_process==0)
{

      develop_neqn();
      construct_fe_collection();
            
      max_speed = compute_max_speed();
      Assert(n_eqn.size() == fe.size(), ExcNotImplemented());

      cellwise_sol.resize(triangulation.n_active_cells()); // the size of triangulation does not change
      store_cell_index_center();
}


template<int dim>
void 
Solve_System_SS_adaptive<dim>::distribute_dofs()
{
      dof_handler.distribute_dofs(fe);
    
      locally_owned_dofs = dof_handler.locally_owned_dofs();
      DoFTools::extract_locally_relevant_dofs (dof_handler,
                          locally_relevant_dofs);
    
      locally_owned_solution.reinit(locally_owned_dofs,mpi_comm);
      locally_owned_residual.reinit(locally_owned_dofs,mpi_comm);
      locally_relevant_solution.reinit(locally_owned_dofs,locally_relevant_dofs,mpi_comm);

      // allocate the memory error 
      IndexSet locally_owned_cells(dof_handler.get_triangulation().n_active_cells()); // helpful while computing error
      create_IndexSet_triangulation(locally_owned_cells); 
    
}

template<int dim>
int
Solve_System_SS_adaptive<dim>::solve_steady_state(Triangulation<dim> &triangulation,
                                                  double &t,const std::vector<Vector<double>> &force_vector) 
{
        // an approximation to delta_t
        double min_length = min_h(triangulation);
        Utilities::MPI::min(min_length,mpi_comm);

        Assert(max_speed.size() != 0,ExcNotInitialized());
        current_max_speed = max_speed[current_max_fe_index()]; // speed corresponding to the current maximum system we have

        dt = CFL * min_length/(dim * current_max_speed);

        pout << "total cells: " << triangulation.n_global_active_cells() << std::endl;
        pout << "total dofs: " << dof_handler.n_dofs() << std::endl;

        typedef FilteredIterator<typename hp::DoFHandler<dim>::active_cell_iterator> CellFilter;

        std::vector<Vector<double>> component = return_component();

        std::vector<std::vector<Vector<double>>> component_to_system = return_component_to_system();

        Assert(component.size() !=0 , ExcNotInitialized());
        Assert(component_to_system.size() !=0 , ExcNotInitialized());

        QGauss<dim-1> face_quadrature_basic(1);
        hp::QCollection<dim-1> face_quadrature;

        for (unsigned int i = 0 ; i < fe.size() ; i++)
          face_quadrature.push_back(face_quadrature_basic);

        PerCellAssemble per_cell_assemble;
        PerCellAssembleScratch per_cell_assemble_scratch(fe,face_quadrature);

        int step_count = 0;
        double residual_ss = 100;

        std::vector<std::vector<Vector<double>>> g(n_eqn.size());

        for (unsigned int i = 0 ; i < n_eqn.size() ; i ++)
        {
          g[i].resize(4);
          for (unsigned int id = 0 ; id < 4 ; id++)
            initial_boundary->bc_inhomo(system_matrices[i].B[id],id,g[i][id],t);
        }


        while (step_count < 100 || residual_ss > 1e-7 ) // we atleast run till t_end || residual_ss > 1e-8
        {
                WorkStream::run ( CellFilter(IteratorFilters::LocallyOwnedCell(),
                  dof_handler.begin_active()),
                CellFilter(IteratorFilters::LocallyOwnedCell(), 
                  dof_handler.end()),
                std::bind(&Solve_System_SS_adaptive<dim>::assemble_per_cell,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3,
                  std::cref(component_to_system),
                  std::cref(t),
                  std::cref(g),
                  std::cref(force_vector)),
                std::bind(&Solve_System_SS_adaptive<dim>::assemble_to_global,
                  this,
                  std::placeholders::_1,
                  std::cref(component)),
                per_cell_assemble_scratch,
                per_cell_assemble);

                residual_ss = locally_owned_residual.l2_norm();
                locally_relevant_solution = locally_owned_solution;

                t += dt;
                step_count++;

      if (step_count%500 == 0)
        pout<< "residual: " << residual_ss << "time: " << t << std::endl;

    }

    return(step_count);

}

template<int dim>
void
Solve_System_SS_adaptive<dim>::run_time_loop(Triangulation<dim> &triangulation,
                                             const unsigned int &cycle,
                                             const unsigned int &refine_cycles,
                                             double &t,
                                             const std::vector<Vector<double>> &force_vector)
{      
      std::vector<std::vector<Vector<double>>> component_to_system = return_component_to_system();   

      locally_relevant_solution = locally_owned_solution; // the value set to locally owned solution is set here
    

        int total_steps = 0;
        std::cout << "cycle: " << cycle << std::endl;
        fflush(stdout);
      
        total_steps = solve_steady_state(triangulation,t,force_vector);

        pout << "Steps taken: " << total_steps << std::endl;

        if(cycle != refine_cycles-1)  // if we are not in the last stage
        {
          create_cellwise_solution(locally_owned_solution,cellwise_sol,component_to_system);
          cell_fe_index = create_cell_fe_index(); // fe index of every cell
        }


}


// prescribe the initial conditions
template<int dim>
void
Solve_System_SS_adaptive<dim>::prescribe_initial_conditions()
{

  typedef FilteredIterator<typename hp::DoFHandler<dim>::active_cell_iterator> CellFilter;

  PerCellIC percellIC;
  PerCellICScratch percellICScratch;

  std::vector<Vector<double>> component = return_component();
  

  WorkStream::run (CellFilter(IteratorFilters::LocallyOwnedCell(),
                    dof_handler.begin_active()),
                   CellFilter(IteratorFilters::LocallyOwnedCell(), dof_handler.end()),
                   std::bind(&Solve_System_SS_adaptive<dim>::compute_ic_per_cell,
                            this,
                            std::placeholders::_1,
                            std::placeholders::_2,
                            std::placeholders::_3,
                            std::cref(component)),
                   std::bind(&Solve_System_SS_adaptive<dim>::copy_ic_to_global,
                            this,
                            std::placeholders::_1),
                   percellICScratch,
                   percellIC);


}

template<int dim>
void
Solve_System_SS_adaptive<dim>::compute_ic_per_cell(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
                                        PerCellICScratch &scratch,
                                        PerCellIC &data,
                                        const std::vector<Vector<double>> &component)
{
  
          data.dofs_per_cell = fe[cell->active_fe_index()].dofs_per_cell;
          data.local_dof_indices.resize(data.dofs_per_cell);
          data.local_contri.reinit(data.dofs_per_cell);

          const Point<dim> location = cell->center();

          for (unsigned int dof = 0 ; dof < data.dofs_per_cell ; dof++)
            data.local_contri[dof] = initial_boundary->ic(location,component[cell->active_fe_index()](dof));
          
          cell->get_dof_indices(data.local_dof_indices);

  
}

template<int dim>
void
Solve_System_SS_adaptive<dim>::copy_ic_to_global(const PerCellIC &data)
{
    for(unsigned int dof = 0 ; dof< data.dofs_per_cell ; dof++)
      locally_owned_solution(data.local_dof_indices[dof]) = data.local_contri(dof);
}


template<int dim>
void 
Solve_System_SS_adaptive<dim>::assemble_per_cell(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
                                      PerCellAssembleScratch &scratch,
                                      PerCellAssemble &data,
                                      const std::vector<std::vector<Vector<double>>> &component_to_system,
                                      const double &t,
                                      const std::vector<std::vector<Vector<double>>> &g,
                                      const std::vector<Vector<double>> &force_vector)
 {

              // operations to avoid data races
              std::vector<std::vector<Vector<double>>> temp_g = g;

              // no hanging nodes check
              Assert(!cell->has_children(), ExcInternalError());
              
              const unsigned int this_fe_index = cell->active_fe_index();

              data.this_fe_index = this_fe_index;
              data.dofs_per_cell = fe[this_fe_index].dofs_per_cell;
              data.local_dof_indices.resize(data.dofs_per_cell);
              cell->get_dof_indices(data.local_dof_indices);
              const double volume = cell->measure();

              data.local_contri.reinit(data.dofs_per_cell);
              data.local_contri = 0;


              // contribution from collisions P
              for (unsigned int m = 0 ; m < system_matrices[this_fe_index].P.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices[this_fe_index].P,m); n ; ++n)
                  {
              // 0 because of finite volume
                    unsigned int dof_sol = component_to_system[this_fe_index][n.row()](0);

              // the solution id which meets An
                    unsigned int dof_sol_col = data.local_dof_indices[component_to_system[this_fe_index][n.col()](0)];

              // explicit euler update, 
                    data.local_contri(dof_sol) +=  n.value() * dt
                                             * locally_relevant_solution(dof_sol_col);

                  }

              if(system_matrices[this_fe_index].have_force)
              {
                Vector<double> force_value(n_eqn[this_fe_index]);
                initial_boundary->force(force_value,force_vector[cell->index()],
                                        cell->center(),t);

                for(unsigned int i = 0 ; i < n_eqn[this_fe_index] ; i++)
                    data.local_contri(component_to_system[this_fe_index][i](0)) += dt * force_value(i);                
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

                const typename hp::DoFHandler<dim>::face_iterator face_itr = cell->face(face);
                const double face_length = face_itr->measure();

                  // construct An of the current cell
                  Sparse_Matrix An_cell = system_matrices[this_fe_index].Ax * nx
                                          + system_matrices[this_fe_index].Ay * ny;

                if (face_itr->at_boundary())
                {

                  const unsigned int bc_id = face_itr->boundary_id();
                  // only compute for times greater than zero, already computed for t= 0 before
                    if (system_matrices[this_fe_index].bc_inhomo_time && t > 1e-16)
                                initial_boundary->bc_inhomo(system_matrices[this_fe_index].B[bc_id],bc_id,
                                                          temp_g[this_fe_index][bc_id],t);


                  for (unsigned int m = 0 ; m < An_cell.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An_cell,m); n ; ++n)
                  {
                      // 0 because of finite volume
                      unsigned int dof_sol = component_to_system[this_fe_index][n.row()](0);

                     // the solution id which meets An
                      unsigned int dof_sol_col = data.local_dof_indices[component_to_system[this_fe_index][n.col()](0)];

                    // explicit euler update
                      data.local_contri(dof_sol) -=  n.value() * dt
                                             * locally_relevant_solution(dof_sol_col) * face_length/volume;

                  }

                  // contribution from penalty_B
                  for (unsigned int m = 0 ; m < system_matrices[this_fe_index].penalty_B[bc_id].outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices[this_fe_index].penalty_B[bc_id],m); n ; ++n)
                  {
                    unsigned int dof_sol = component_to_system[this_fe_index][n.row()](0);
                    unsigned int dof_sol_col = data.local_dof_indices[component_to_system[this_fe_index][n.col()](0)];
                    data.local_contri(dof_sol) +=  n.value() * dt
                                             * locally_relevant_solution(dof_sol_col) * face_length/volume;

                  }
                  // contribution from penalty * g
                  for (unsigned int m = 0 ; m < system_matrices[this_fe_index].penalty[bc_id].outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices[this_fe_index].penalty[bc_id],m); n ; ++n)
                  {
              
                    unsigned int dof_sol = component_to_system[this_fe_index][n.row()](0);
                    unsigned int dof_sol_col = data.local_dof_indices[component_to_system[this_fe_index][n.col()](0)];
                    data.local_contri(dof_sol) -=  n.value() * dt
                                                  * temp_g[this_fe_index][bc_id](n.col()) * face_length/volume;

                  }
                }
                else
                {

                   typename hp::DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face);
                   const unsigned int neighbor_fe_index = neighbor->active_fe_index();
                   const unsigned int dofs_per_cell_neighbor = fe[neighbor_fe_index].dofs_per_cell;
                   std::vector<types::global_dof_index> local_dof_indices_neighbor(dofs_per_cell_neighbor);


                   neighbor->get_dof_indices(local_dof_indices_neighbor);

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
                  {
                    unsigned int dof_sol = component_to_system[this_fe_index][n.row()](0);
                    unsigned int dof_sol_col = data.local_dof_indices[component_to_system[this_fe_index][n.col()](0)];
                    data.local_contri(dof_sol) -=  dt * n.value() * (locally_relevant_solution(dof_sol_col))
                                                   * face_length/(2 * volume);

                  }

                  // we add the diffusion now
                  for (unsigned int m = 0 ; m < Amod.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(Amod,m); n ; ++n)
                      if(n.row() < n_eqn[this_fe_index] && n.col() < n_eqn[this_fe_index])
                  {
                    unsigned int dof_sol = component_to_system[this_fe_index][n.row()](0);
                    unsigned int dof_sol_col = data.local_dof_indices[component_to_system[this_fe_index][n.col()](0)];
                    data.local_contri(dof_sol) -=  dt * n.value() * (locally_relevant_solution(dof_sol_col))
                                                   * face_length/(2 * volume);

                  }

                  // contribution from the neighboring cell
                  for (unsigned int m = 0 ; m < An_effective.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An_effective,m); n ; ++n)
                  {
                    unsigned int dof_sol = component_to_system[this_fe_index][n.row()](0);
                    unsigned int dof_neighbor_col = local_dof_indices_neighbor[component_to_system[neighbor_fe_index][n.col()](0)];
                    data.local_contri(dof_sol) -=  dt * n.value() * (locally_relevant_solution(dof_neighbor_col))
                                                   * face_length/(2 * volume);

                  }

                  // diffusion with the neighbouring cell
                  for (unsigned int m = 0 ; m < Amod.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(Amod,m); n ; ++n)
                      if(n.row() < n_eqn[this_fe_index] && n.col() < n_eqn[neighbor_fe_index])
                  {
                    unsigned int dof_sol = component_to_system[this_fe_index][n.row()](0);
                    unsigned int dof_neighbor_col = local_dof_indices_neighbor[component_to_system[neighbor_fe_index][n.col()](0)];
                    data.local_contri(dof_sol) -=  dt * n.value() * (-locally_relevant_solution(dof_neighbor_col))
                                                   * face_length/(2 * volume);

                  }

                } //end of else


              } 
              //end of loop over the faces


 }


template<int dim>
void
Solve_System_SS_adaptive<dim>::assemble_to_global(const PerCellAssemble &data,const std::vector<Vector<double>> &component)
 {
        for (unsigned int i = 0 ; i < data.dofs_per_cell ; i ++)
        {
              locally_owned_solution(data.local_dof_indices[i]) = locally_relevant_solution(data.local_dof_indices[i])
                                                               + data.local_contri(i);

              if (component[data.this_fe_index][i] <= 7) // only count till tensor degree 2
                locally_owned_residual(data.local_dof_indices[i]) = fabs(data.local_contri(i))/dt;          
        }

 }

template<int dim>
void 
Solve_System_SS_adaptive<dim>::create_output(const std::string &filename)
{

  typename hp::DoFHandler<dim>::active_cell_iterator endc = dof_handler.end(), cell = dof_handler.begin_active();
  std::vector<std::vector<Vector<double>>> component_to_system = return_component_to_system();
  const unsigned int to_print = 13;   // print till the third order tensor

  FILE *fp;
  fp = fopen(filename.c_str(),"w+");

  AssertThrow(fp!=NULL,ExcFileNotOpen(filename));
  for(; cell != endc ; cell++)
    {
      const unsigned int this_fe_index = cell->active_fe_index();
      std::vector<types::global_dof_index> local_dof_indices(fe[this_fe_index].dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int space = 0 ; space < dim ; space ++)
        fprintf(fp, "%0.16f\t",cell->center()(space));

      for (unsigned int i = 0 ; i < n_eqn[this_fe_index] ; i++) 
        if(i < to_print)
      {
        const double sol_value =  locally_owned_solution(local_dof_indices[component_to_system[this_fe_index][i](0)]);
        fprintf(fp, "%0.16f\t",sol_value);
      }

      fprintf(fp, "\n");
    }

  fclose(fp);

}

// template<int dim>
// void 
// Solve_System_SS_adaptive<dim>::compute_error_per_cell(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
//                                              PerCellErrorScratch &scratch,
//                                              PerCellError &data)
// {

//     data.error_value = 0;
//     data.solution_value = 0;

//     VectorTools::point_value(cell->get_dof_handler(),locally_owned_solution, 
//                              cell->center(),data.solution_value); 
//     initial_boundary->exact_solution(cell->center(),data.exact_solution,t_end);

//     // overwrite the solution value with the error in this cell
//     data.solution_value.sadd(1,-1,data.exact_solution);
    
//     for (unsigned int i = 0 ; i < data.solution_value.size(); i++) 
//       data.error_value = pow(data.solution_value(i),2) +data.error_value;

//     data.cell_index = cell->index();

//     data.volume = cell->measure();
    
  
// }

// template<int dim>
// void
// Solve_System_SS<dim>::copy_error_to_global(const PerCellError &data)
// {
//   error_per_cell(data.cell_index) = sqrt(data.error_value * data.volume);
// }


// template<int dim>
// void
// Solve_System_SS<dim>::compute_error()
// {
//   PerCellError per_cell_error(n_eqn);
//   PerCellErrorScratch per_cell_error_scratch;

//   typedef FilteredIterator<typename DoFHandler<dim>::active_cell_iterator> CellFilter;

//   WorkStream::run ( CellFilter(IteratorFilters::LocallyOwnedCell(),
//                     dof_handler.begin_active()),
//                     CellFilter(IteratorFilters::LocallyOwnedCell(), 
//                     dof_handler.end()),
//                     std::bind(&Solve_System_SS<dim>::compute_error_per_cell,
//                             this,
//                             std::placeholders::_1,
//                             std::placeholders::_2,
//                             std::placeholders::_3),
//                     std::bind(&Solve_System_SS<dim>::copy_error_to_global,
//                             this,
//                             std::placeholders::_1),
//                    per_cell_error_scratch,
//                    per_cell_error);
  
//   double errorTot = error_per_cell.l2_norm();
//   pout << dof_handler.n_dofs() << "\t" << errorTot << std::endl;

// }


template<int dim>
void
Solve_System_SS_adaptive<dim>::create_IndexSet_triangulation(IndexSet &locally_owned_cells)
{
  typename hp::DoFHandler<dim>::active_cell_iterator endc = dof_handler.end(), 
                                                 cell = dof_handler.begin_active();

  for(; cell!=endc;cell++)  
      locally_owned_cells.add_index(cell->index());
    
}


template<int dim>
std::vector<double>
Solve_System_SS_adaptive<dim>::compute_max_speed()
{
    std::vector<double> max_speed(n_eqn.size());

    for(unsigned int i = 0 ; i < n_eqn.size(); i++)
    {
      EigenSolver<MatrixXd> ES(system_matrices[i].Ax);
      VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal().cwiseAbs();      
      max_speed[i] = vals.maxCoeff();
    }


    return(max_speed);
}


template<int dim>
double
Solve_System_SS_adaptive<dim>::min_h(Triangulation<dim> &triangulation)
{

  typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(),
                                                                           endc = triangulation.end();


  double min_length = 10000;

  for(; cell != endc ; cell++)
      for(unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
        if (cell->face(face)->measure() < min_length)
          min_length = cell->face(face)->measure();


  return(min_length);
}

template<int dim>
Solve_System_SS_adaptive<dim>::PerCellAssembleScratch::PerCellAssembleScratch(const hp::FECollection<dim> &fe,
                                                                              const hp::QCollection<dim-1> &quadrature)
:
fe_v_face(fe,quadrature,update_normal_vectors)
{;}

template<int dim>
Solve_System_SS_adaptive<dim>::PerCellAssembleScratch::PerCellAssembleScratch(const PerCellAssembleScratch &scratch)
:
fe_v_face(scratch.fe_v_face.get_fe_collection(),
          scratch.fe_v_face.get_quadrature_collection(),
          update_normal_vectors)
{;}

template<int dim>
void 
Solve_System_SS_adaptive<dim>::develop_neqn()
{
  n_eqn.resize(system_matrices.size());

  for(unsigned int i = 0 ; i < system_matrices.size() ; i++)
    n_eqn[i] = system_matrices[i].Ax.rows();

}

template<int dim>
void 
Solve_System_SS_adaptive<dim>::construct_block_structure(std::vector<int> &block_structure,
                                                         const std::vector<unsigned int> &n_eqn)
{
    AssertThrow(n_eqn.size() != 0, ExcNotInitialized());
    AssertThrow(std::is_sorted(std::begin(n_eqn),std::end(n_eqn)),
                ExcMessage("number of equations not sorted"));

    // the very first entry should be the number of equations in the first system
    block_structure.push_back(n_eqn[0]);

    for (unsigned long int i = 1 ; i < n_eqn.size() ; i++)
      block_structure.push_back(n_eqn[i]-n_eqn[i-1]);

    AssertDimension(block_structure.size(),n_eqn.size());
}

template<int dim>
void
Solve_System_SS_adaptive<dim>::construct_fe_collection()
{
    std::vector<int> block_structure;

    Assert(n_eqn.size() != 0,ExcNotInitialized());

    // first create the block structure for the finite element object 
    // which will then be used to construct the fe system

    construct_block_structure(block_structure,n_eqn);

    std::vector< std::vector< const FiniteElement<dim>* > > fe_list;
    std::vector<unsigned int> multiplicities;

    fe_list.resize(block_structure.size());
    
    for(unsigned int i = 0 ; i < block_structure.size() ; i++)
      multiplicities.push_back(block_structure[i]);

    switch (block_structure.size())
    {
      case 1:
      {
        fe_list[0].push_back(&fe_basic);
        break;
      }

      case 2:
      {
        fe_list[0].push_back(&fe_basic);
        fe_list[0].push_back(new FE_Nothing<dim>());

        fe_list[1].push_back(&fe_basic);
        fe_list[1].push_back(&fe_basic);

        break;
      }

      case 3:
      {
        fe_list[0].push_back(&fe_basic);
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());

        fe_list[1].push_back(&fe_basic);
        fe_list[1].push_back(&fe_basic);
        fe_list[1].push_back(new FE_Nothing<dim>());

        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(&fe_basic);
        break;
      }

      case 4:
      {
        fe_list[0].push_back(&fe_basic);
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());

        fe_list[1].push_back(&fe_basic);
        fe_list[1].push_back(&fe_basic);
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());

        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(new FE_Nothing<dim>());

        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        break;
      }

      case 5:
      {
        fe_list[0].push_back(&fe_basic);
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());

        fe_list[1].push_back(&fe_basic);
        fe_list[1].push_back(&fe_basic);
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());

        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(new FE_Nothing<dim>());
        fe_list[2].push_back(new FE_Nothing<dim>());

        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(new FE_Nothing<dim>());


        fe_list[4].push_back(&fe_basic);
        fe_list[4].push_back(&fe_basic);
        fe_list[4].push_back(&fe_basic);
        fe_list[4].push_back(&fe_basic);
        fe_list[4].push_back(&fe_basic);
        break;
      }

      case 6:
      {
        fe_list[0].push_back(&fe_basic);
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
	fe_list[0].push_back(new FE_Nothing<dim>());

        fe_list[1].push_back(&fe_basic);
        fe_list[1].push_back(&fe_basic);
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());
	fe_list[1].push_back(new FE_Nothing<dim>());

        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(new FE_Nothing<dim>());
        fe_list[2].push_back(new FE_Nothing<dim>());
	fe_list[2].push_back(new FE_Nothing<dim>());

        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(new FE_Nothing<dim>());
	fe_list[3].push_back(new FE_Nothing<dim>());


        fe_list[4].push_back(&fe_basic);
        fe_list[4].push_back(&fe_basic);
        fe_list[4].push_back(&fe_basic);
        fe_list[4].push_back(&fe_basic);
        fe_list[4].push_back(&fe_basic);
	fe_list[4].push_back(new FE_Nothing<dim>());
 

        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
	fe_list[5].push_back(&fe_basic);

       break;
      }

      case 7:
      {
        fe_list[0].push_back(&fe_basic);
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());

        fe_list[1].push_back(&fe_basic);
        fe_list[1].push_back(&fe_basic);
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());

        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(new FE_Nothing<dim>());
        fe_list[2].push_back(new FE_Nothing<dim>());
        fe_list[2].push_back(new FE_Nothing<dim>());
        fe_list[2].push_back(new FE_Nothing<dim>());

        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(new FE_Nothing<dim>());
        fe_list[3].push_back(new FE_Nothing<dim>());
        fe_list[3].push_back(new FE_Nothing<dim>());


        fe_list[4].push_back(&fe_basic);
        fe_list[4].push_back(&fe_basic);
        fe_list[4].push_back(&fe_basic);
        fe_list[4].push_back(&fe_basic);
        fe_list[4].push_back(&fe_basic);
        fe_list[4].push_back(new FE_Nothing<dim>());
        fe_list[4].push_back(new FE_Nothing<dim>());
 

        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(new FE_Nothing<dim>());

        fe_list[6].push_back(&fe_basic);
        fe_list[6].push_back(&fe_basic);
        fe_list[6].push_back(&fe_basic);
        fe_list[6].push_back(&fe_basic);
        fe_list[6].push_back(&fe_basic);
        fe_list[6].push_back(&fe_basic);
        fe_list[6].push_back(&fe_basic);

       break;
      }

      case 8:
      {
        fe_list[0].push_back(&fe_basic);
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());
        fe_list[0].push_back(new FE_Nothing<dim>());

        fe_list[1].push_back(&fe_basic);
        fe_list[1].push_back(&fe_basic);
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());
        fe_list[1].push_back(new FE_Nothing<dim>());

        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(&fe_basic);
        fe_list[2].push_back(new FE_Nothing<dim>());
        fe_list[2].push_back(new FE_Nothing<dim>());
        fe_list[2].push_back(new FE_Nothing<dim>());
        fe_list[2].push_back(new FE_Nothing<dim>());
        fe_list[2].push_back(new FE_Nothing<dim>());

        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(&fe_basic);
        fe_list[3].push_back(new FE_Nothing<dim>());
        fe_list[3].push_back(new FE_Nothing<dim>());
        fe_list[3].push_back(new FE_Nothing<dim>());
        fe_list[3].push_back(new FE_Nothing<dim>());


        fe_list[4].push_back(&fe_basic);//1
        fe_list[4].push_back(&fe_basic);//2
        fe_list[4].push_back(&fe_basic);//3
        fe_list[4].push_back(&fe_basic);//4
        fe_list[4].push_back(&fe_basic);//5
        fe_list[4].push_back(new FE_Nothing<dim>());//6
        fe_list[4].push_back(new FE_Nothing<dim>());//7
        fe_list[4].push_back(new FE_Nothing<dim>());//8
 

        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(&fe_basic);
        fe_list[5].push_back(new FE_Nothing<dim>());
        fe_list[5].push_back(new FE_Nothing<dim>());

        fe_list[6].push_back(&fe_basic);
        fe_list[6].push_back(&fe_basic);
        fe_list[6].push_back(&fe_basic);
        fe_list[6].push_back(&fe_basic);
        fe_list[6].push_back(&fe_basic);
        fe_list[6].push_back(&fe_basic);
        fe_list[6].push_back(&fe_basic);
        fe_list[6].push_back(new FE_Nothing<dim>());


        fe_list[7].push_back(&fe_basic);
        fe_list[7].push_back(&fe_basic);
        fe_list[7].push_back(&fe_basic);
        fe_list[7].push_back(&fe_basic);
        fe_list[7].push_back(&fe_basic);
        fe_list[7].push_back(&fe_basic);
        fe_list[7].push_back(&fe_basic);
        fe_list[7].push_back(&fe_basic);

       break;
      }

      default:
      {
        AssertThrow(1 == 0, ExcMessage("Should not have reached here"));
        break;
      }
    }

    for(unsigned int i = 0 ; i < fe_list.size() ;i++)
      fe.push_back(FESystem<dim>(fe_list[i],multiplicities));
}

template<int dim>
unsigned int 
Solve_System_SS_adaptive<dim>::current_max_fe_index()
{
    typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();

    unsigned int fe_index_max = 0;

    for(; cell != endc ; cell++)
        if(cell->active_fe_index() > fe_index_max)
          fe_index_max = cell->active_fe_index();


    return(fe_index_max);

}

// set the fe index as the present cycle
template<int dim>
void
Solve_System_SS_adaptive<dim>::allocate_fe_index(const unsigned int present_cycle,
                                                 const Vector<double> &error_per_cell,
                                                 Triangulation<dim> &triangulation)
{
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  const double frac_refine = 0.3;
  const double frac_coarse = 0;
  bool found_higher = false;

  //if(present_cycle == 0)
    for(; cell != endc ; cell++)
            cell->set_active_fe_index(present_cycle);
/*  else
  {
    // we use the triangulation object to mark cells for refinement
    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                            error_per_cell,
                                                            frac_refine,
                                                            frac_coarse);  
    for(; cell != endc ; cell++)
      {
        const unsigned int current_index = cell->active_fe_index();
        AssertThrow(current_index <=5 ,ExcNotImplemented());


        if(current_index < fe.size() && cell->refine_flag_set()) // else do nothing
        {
          cell->set_active_fe_index(current_index + 1);           // increase the current fe index
          cell->clear_refine_flag();                              // dont need the refinement flag anymore
        }

      }
  }
*/
            
}

template<int dim>
std::vector<Vector<double>> 
Solve_System_SS_adaptive<dim>::return_component()
{
    std::vector<Vector<double>> component(n_eqn.size());
  
  for(unsigned int i = 0 ; i < component.size(); i++)
  {
    component[i].reinit(fe[i].dofs_per_cell);
    for (unsigned int j = 0 ; j < fe[i].dofs_per_cell ; j++)
        component[i](j) = fe[i].system_to_component_index(j).first;
  }

  return(component);
}

template<int dim>
std::vector<std::vector<Vector<double>>>
Solve_System_SS_adaptive<dim>::return_component_to_system()
{
        std::vector<std::vector<Vector<double>>> component_to_system(n_eqn.size());
        const unsigned int dofs_per_comp = 1;     // considering a finite volume scheme

      // loop over all the different systems available in the system
      for (unsigned long int i = 0 ; i < n_eqn.size(); i++)
        // loop over all the equations of the particular system
        for (unsigned int j = 0 ; j < n_eqn[i]; j++)
          component_to_system[i].push_back(Vector<double>(dofs_per_comp));

      // now we allocate the values for component_to_system
      for (unsigned long int i = 0 ; i < n_eqn.size(); i++)
          for (unsigned int k = 0 ; k < n_eqn[i] ;k ++)
              for (unsigned int j = 0 ; j < dofs_per_comp ; j++)
                component_to_system[i][k](j) = fe[i].component_to_system_index(k,j);


      return(component_to_system);
}

// template<int dim>
// void 
// Solve_System_SS_adaptive<dim>::refine_and_interpolate(parallel::distributed::Triangulation<dim> &triangulation)
// {
//        //we now need to interpolate onto the new mesh
//      locally_relevant_solution_temp = locally_owned_solution; // initialise the temporary

//      parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> soltrans(dof_handler);


//      typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end(), 
//      cell = dof_handler.begin_active();  

//      for(; cell != endc ; cell++)
//       if (cell->is_locally_owned())
//         cell->set_refine_flag();

//       soltrans.prepare_for_coarsening_and_refinement(locally_relevant_solution_temp);
//       triangulation.refine_global();
//       distribute_dofs();        // reinitialise all the vectors
//       soltrans.interpolate(locally_owned_solution); // interpolate to the new vector

//       locally_owned_solution.compress(VectorOperation::insert);

// }


// returns the effective An of the neighbor
template<int dim>
Sparse_Matrix
Solve_System_SS_adaptive<dim>::construct_An_effective(const Sparse_Matrix &An_cell,const Sparse_Matrix &An_neighbor)
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

// store the center corresponding to the index
template<int dim>
void 
Solve_System_SS_adaptive<dim>::store_cell_index_center()
{
      typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin(), endc = dof_handler.end();

      for(; cell != endc ; cell++)
        cell_index_center.push_back(cell->center());
}

template<int dim>
std::vector<unsigned int>
Solve_System_SS_adaptive<dim>::create_cell_fe_index()
{
  std::vector<unsigned int> cell_fe_index(dof_handler.get_triangulation().n_active_cells());

  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin(), endc = dof_handler.end();  

  for(; cell != endc ; cell++)
    cell_fe_index[cell->index()] = cell->active_fe_index();

  return(cell_fe_index);

}
// takes in a solution based upon the dof and returns the cellwise solution
template<int dim>
void 
Solve_System_SS_adaptive<dim>::create_cellwise_solution(const LA::MPI::Vector &dofwise_sol,
                                                        std::vector<Vector<double>> &cellwise_sol,
                                                        const std::vector<std::vector<Vector<double>>> &component_to_system)
{

    typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin(), endc = dof_handler.end();  

    Assert(cellwise_sol.size() != 0,ExcNotInitialized());
    Assert(dofwise_sol.size() == dof_handler.n_dofs(),ExcNotInitialized());

    for(; cell != endc ; cell++)
    {
      std::vector<types::global_dof_index> local_dof_indices(cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
      const unsigned int this_fe_index = cell->active_fe_index();
      cellwise_sol[cell->index()].reinit(n_eqn[this_fe_index]);

      for(unsigned int i = 0 ; i < n_eqn[this_fe_index]; i++)
        cellwise_sol[cell->index()](i) = dofwise_sol(local_dof_indices[component_to_system[this_fe_index][i](0)]);

    }
}

template<int dim>
void 
Solve_System_SS_adaptive<dim>::interpolate_higher_fe_index(const std::vector<Vector<double>> &cellwise_sol,
                                                           const std::vector<unsigned int> &cell_fe_index,
                                                           LA::MPI::Vector &new_vec,
                                                           const std::vector<std::vector<Vector<double>>> &component_to_system)
{
    Assert(cellwise_sol.size() != 0,ExcNotInitialized());
    typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin(), endc = dof_handler.end();
    Assert(dof_handler.n_dofs() == new_vec.size(),ExcNotInitialized());
    Assert(cell_index_center.size() != 0 ,ExcNotInitialized());

    for(;cell != endc ; cell++)
    {
      // check whether the indexing has been preserved
      AssertThrow(cell_index_center[cell->index()].distance(cell->center()) < 1e-16,
                  ExcOrderingChanged(cell_index_center[cell->index()].distance(cell->center())));

      // or not during polynomial increase.
      std::vector<types::global_dof_index> local_dof_indices(cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
      const unsigned int this_fe_index = cell->active_fe_index();

      // loop over the old fe index
      for(unsigned int i = 0 ; i < n_eqn[cell_fe_index[cell->index()]] ; i ++)
        new_vec(local_dof_indices[component_to_system[this_fe_index][i](0)]) = cellwise_sol[cell->index()](i);


    }


}
template<int dim>
Solve_System_SS_adaptive<dim>::~Solve_System_SS_adaptive()
{
  dof_handler.clear();
}

// explicit initiation to avoid linker error
template class Solve_System_SS_adaptive<2>;


