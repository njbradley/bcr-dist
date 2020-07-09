#ifndef FILEIO_C
#define FILEIO_C

#include "fileio.h"
#include "cell.h"
#include "table.h"

void load_bd_data(string path_heavy, string path_light, vector<bcell_double>& cells) {
	itablestream itable_heavy(path_heavy);
	itablestream itable_light(path_light);
	tablerow hrow(&itable_heavy);
	tablerow lrow(&itable_light);
	
	while (!hrow.eof and !lrow.eof) {
		string id = hrow.get("Cell Label");
		if (lrow.get("Cell Label") != id) {
			int hid = atoi(id.c_str());
			int lid = atoi(lrow.get("Cell Label").c_str());
			if (hid > lid) {
				lrow = tablerow(&itable_light);
				continue;
			} else {
				hrow = tablerow(&itable_heavy);
				continue;
			}
		}
		string hcdr3 = hrow.get("AA CDR3");
		string lcdr3 = lrow.get("AA CDR3");
		if (hcdr3 != "[CDR3_not_canonical]" and lcdr3 != "[CDR3_not_canonical]") {
			cells.emplace_back(id, bcell_chain(hrow.get("V"), hcdr3), bcell_chain(lrow.get("V"), lcdr3));
		}
		hrow = tablerow(&itable_heavy);
		lrow = tablerow(&itable_light);
	}
	cout << "sucessfully loaded " << cells.size() << " cell rows from files " << path_heavy << ' ' << path_light << endl;
}



void save_dist_matrix(string path, vector<bcell_double>& cells) {
	vector<string> ids {"cell-id"};
	for (bcell_double& cell : cells) {
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
