#ifndef DISTANCES_CC
#define DISTANCES_CC

dist_solver::dist_solver(int ncomponents, int numcells): ndims(ncomponents), ncells(numcells) {
	matrix = new double[ndims*numcells];
}

dist_solver::fit_transform(function<double(int,int)> dist_func) {
	for (int i = 0; i < numcells; i ++) {
		double pos[ndims];
		double dists[i];
		for (int k = 0; k < ndims; k ++) {
			pos = 0;
		}
		for (int j = 0; j < i; j ++) {
			dists[j] = dist_func(i,j);
		}
	}
}
		
			




cluster_solver::cluster_solver(int newmax_cells): max_cells(newmax_cells) {
	
}

void cluster_solver::fit_transform(vector<bcell_double>& cells, vector<bcell_double*>& outcells, vector<string,vector<string> >& id_mapping) {
	int av_size = cells.size()/max_cells;
	double av_close_dist = 0;
	for (bcell_double& cell : cells) {
		double min = cell.distance(&cells[0]);
		for (int i = 1; i < cells.size(); i ++) {
			double dist = cell.distance(&cells[i]);
			if (dist < min) {
				min = dist;
			}
		}
		av_close_dist += min;
	}
	av_close_dist /= cells.size();
	
	vector<bcell_double*> cells_left;
	
	for (bcell_double& cell : cells) {
		cells_left.push_back(&cell);
	}
	
	while (cells_left.size() + outcells.size() > max_cells) {
		vector<bcell_double*> cluster;
		vector<bcell_double*> lastpoints {rand()%cells_left.size()};
		vector<bcell_double*> newpoints;
		
		// for (bcell_double* cell : lastpoints) {
		// 	for (
	}
}
	

#endif
