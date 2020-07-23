#ifndef DISTANCES_H
#define DISTANCES_H

#include <functional>

class dist_solver { public:
	double* matrix;
	int ndims;
	int ncells;
	
	dist_solver(int ncomponents, int numcells);
	void fit_transform(function<double(int,int)> dist_func);
};

class cluster_solver { public:
	int max_cells;
	
	cluster_solver(int newmax_cells);
	void fit_transform(vector<bcell_double>& cells, vector<bcell_double>& outcells, vector<string,vector<string> >& id_mapping);
};

#endif
