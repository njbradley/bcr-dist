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
	vector<bcell*> cells;
	vector<bcell*> clusters;
	unordered_map<string,vector<string> > mapping;
	
	cluster_solver(int newmax_cells, vector<bcell*> cells);
	void fit_transform();
};

#endif
