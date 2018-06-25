#include "ic_bc_base.h"

template<int dim>
ic_bc_base<dim>::ic_bc_base()
{;}

// explicit initiation to avoid linker error
template class ic_bc_base<2>;

