#include "include_dealii.h"
#include "EigenSetup.h"

using namespace dealii;


void 
build_triplet(std::vector<triplet> &Row_Col_Value,const std::string filename)
	{

		DeclException1 (ExcFileNotOpen, std::string,
                        << "Could not open " << arg1 );
		
		// as a precautionary measure, remove all the old value
		Row_Col_Value.clear();		

		// number of nonzeros in the system
		unsigned int nz;			

		// create an fstrem to read the file
		std::fstream in(filename.c_str());

		// string which will read the file line by line
		std::string line;

		// check whether we can open the file or not
		AssertThrow(in.is_open(),ExcFileNotOpen(filename));

		// get the first line in "in"
		std::getline(in,line);
		std::stringstream ss(line);
		ss >> nz;

		// shout out if the matrix is empty
		Assert(nz !=0 ,ExcMessage("Total number of non zeros should be > 0"));
		Row_Col_Value.reserve(nz);

		unsigned int counter = 0;
		while (getline(in, line))
		{
			std::stringstream ss2(line);

			unsigned int row;
			unsigned int col;
			double value;

			ss2 >> row;
			ss2 >> col;
			ss2 >> value;
			Row_Col_Value.push_back(triplet(row,col,value));

			counter ++;
		}

		// we check whether the correct number of values have been read or not
		AssertDimension(counter,nz);

	}

	void build_matrix_from_triplet(Sparse_Matrix &matrix,std::vector<triplet> &Row_Col_Value)
	{
		// first we check whether Row_Col_Value has some size or not
		Assert(Row_Col_Value.size() != 0,ExcMessage("Triplet not constructed"));

		// Now we check whether the sparse matrix has been resized or not
		Assert(matrix.rows() != 0 && matrix.cols() !=0,ExcNotInitialized());

		matrix.setFromTriplets(Row_Col_Value.begin(),Row_Col_Value.end());

		// compress the system matrix
		matrix.makeCompressed();

		// clear the memory consumed by the triplet
		Row_Col_Value.clear();
	}