#ifndef FILEIO_C
#define FILEIO_C

#include "fileio.h"

void load_bd_data(string path, vector<bcell>& cells) {
	itablestream itable(path);
	tablerow row(&itable);
	while (!row.eof) {
		cells.emplace_back(row.get("Cell Label"), row.get("V"), row.get("AA CDR3"));
		row = tablerow(&itable);
	}
	cout << "sucessfully loaded " << cells.size() << " cell rows from file " << path << endl;
}

void save_dist_matrix(string path, vector<bcell>& cells) {
	vector<string> ids {"cell-id"};
	for (bcell& cell : cells) {
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
