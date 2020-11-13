#ifndef DISTANCES_CC
#define DISTANCES_CC

#include <LBFGS.h>

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
		
			




cluster_solver::cluster_solver(int newmax_cells, vector<bcell*>& newcells): max_cells(newmax_cells), cells(newcells) {
	
}

void cluster_solver::fit_transform() {
	int av_size = cells.size()/max_cells;
	double av_close_dist = 0;
	for (bcell* cell : cells) {
		double min = cell->distance(cells[0]);
		for (int i = 1; i < cells.size(); i ++) {
			double dist = cell->distance(cells[i]);
			if (dist < min) {
				min = dist;
			}
		}
		av_close_dist += min;
	}
	av_close_dist /= cells.size();
	
	vector<bcell*> cells_left(cells);
	
	while (cells_left.size() + outcells.size() > max_cells) {
		vector<bcell*> cluster;
		int start_index = rand()%cells_left.size();
		vector<bcell*> lastpoints {cells_left[start_index]};
		cells_left.erase(cells_left.begin() + start_index);
		vector<bcell*> newpoints;
		while (lastpoints.size() > 0) {
			for (bcell* lastcell : lastpoints) {
				for (int i = cells_left.size()-1; i >= 0; i --) {
					double dist = lastcell->distance(cells_left[i]);
					if (dist < av_close_dist) {
						newpoints.push_back(dist);
						cells_left.erase(cells_left.begin() + i);
					}
				}
				cluster.push_back(lastcell);
			}
			lastpoints.clear();
			newpoints.swap(lastpoints);
		}
		
		bcell* center = nullptr;
		double dist;
		for (bcell* cell : cluster) {
			double total;
			for (bcell* other : cluster) {
				total += cell->distance(other);
			}
			if (center == nullptr or total < dist) {
				center = cell;
				dist = total;
			}
		}
		
		clusters.push_back(center);
		for (bcell* cell : cluster) {
			if (cell != center) {
				mapping[center->id].push_back(cell->id);
			}
		}
	}
}
	

#endif
