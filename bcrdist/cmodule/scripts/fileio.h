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











template <typename celltype>
void save_dist_matrix(string path, vector<celltype>& cells) {
	vector<string> ids {"cell-id"};
	for (celltype& cell : cells) {
		ids.push_back(cell.id);
	}
	otablestream otable(path, &ids);
	for (int i = cells.size()-1; i >= 0; i --) {
		tablerow row;
		row.add("cell-id", cells[i].id);
		for (int j = 0; j < i; j ++) {
			double dist = cells[i].distance(&cells[j]);
			row.add(cells[j].id, dist);
		}
		otable.writeline(&row);
	}
	cout << "sucessfully saved " << cells.size() << " cells into the distance matrix " << path << endl;
}

#endif
