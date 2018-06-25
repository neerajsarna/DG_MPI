#include "include_dealii.h"

using namespace dealii;
template<int dim>
class
ic_bc_base
{
	public:
		ic_bc_base();

		virtual double ic(const Point<dim> &p,const int &id) = 0;

};