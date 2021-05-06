#ifndef FILEIO_H
#define FILEIO_H

#include "classes.h"
#include "table.h"

bool load_bd_data(string path_heavy, string path_light, vector<dsbcell>& cells);
bool load_bd_data(string path, vector<dsbcell>& cells);

bool load_10x_data(string path, vector<dsbcell>& cells);

void load_dekosky_data(string path, vector<dsbcell>& cells);

template <typename celltype>
void save_dist_matrix(string path, vector<celltype>& cells);

void read_dist_matrix(string path, float* matrix, int dim);








template <typename celltype>
void save_dist_matrix(string path, vector<celltype>& cells) {
	ofstream ofile(path);
	ofile << "cell-id";
	for (int i = 0; i < cells.size(); i ++) {
		ofile << '\t' << cells[i].id;
	}
	ofile << endl;
	for (int i = 0; i < cells.size(); i ++) {
		ofile << cells[i].id;
		for (int j = 0; j < i; j ++) {
			ofile << '\t' << cells[i].distance(&cells[j]);
		}
		ofile << endl;
	}
}

#endif
