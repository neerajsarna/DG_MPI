/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2016 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2009, 2010
 *         Timo Heister, University of Goettingen, 2009, 2010
 */
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/generic_linear_algebra.h>
namespace LA
{
  using namespace dealii::LinearAlgebraTrilinos;
// #if defined(DEAL_II_WITH_PETSC) && !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
//   using namespace dealii::LinearAlgebraPETSc;
// #  define USE_PETSC_LA
// #elif defined(DEAL_II_WITH_TRILINOS)
//   using namespace dealii::LinearAlgebraTrilinos;
// #else
// #  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
// #endif
}
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <fstream>
#include <iostream>


int main(int argc, char *argv[])
{
      using namespace dealii;
      
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  
      const int dim = 2;
      const int p = 0;

      MPI_Comm                                  mpi_communicator(MPI_COMM_WORLD);
      parallel::distributed::Triangulation<dim> triangulation(mpi_communicator);
      DoFHandler<dim>                           dof_handler(triangulation);
      FE_DGQ<dim>                                 fe(p);
      IndexSet                                  locally_owned_dofs;
      IndexSet                                  locally_relevant_dofs;
      LA::MPI::SparseMatrix                     system_matrix;
      const int this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator));

      GridGenerator::subdivided_hyper_cube(triangulation,3);

      dof_handler.distribute_dofs(fe);
      locally_owned_dofs = dof_handler.locally_owned_dofs ();
      DoFTools::extract_locally_relevant_dofs (dof_handler,
                          locally_relevant_dofs) ;
    
    
 
      LA::MPI::Vector locally_owned;
      LA::MPI::Vector locally_relevant;

      locally_owned.reinit(locally_owned_dofs,mpi_communicator);
      locally_relevant.reinit(locally_owned_dofs,locally_relevant_dofs,mpi_communicator);

 
      IndexSet::ElementIterator begin_owned = locally_owned_dofs.begin(), end_owned = locally_owned_dofs.end();
      IndexSet::ElementIterator b_relevant = locally_relevant_dofs.begin(), e_relevant = locally_relevant_dofs.end();

      for (; begin_owned != end_owned; ++begin_owned)
	{
		std::cout << "i am: " << this_mpi_process << std::endl;
		locally_owned((*begin_owned)) = this_mpi_process;
		std::cout << "i have dof: " << (*begin_owned) << "i have value: " << locally_owned((*begin_owned)) << std::endl; 
	}

 	locally_relevant = locally_owned;
     for (; b_relevant != e_relevant; ++b_relevant)
	{
		std::cout << "i am: " << this_mpi_process << std::endl;
		std::cout << "i have dof:" << (*b_relevant) << "i have value: " << locally_relevant((*b_relevant)) << std::endl; 
        }

       typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();


   	
  return 0;
}
