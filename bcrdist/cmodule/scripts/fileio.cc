#ifndef FILEIO_C
#define FILEIO_C

#include "fileio.h"
#include "cell.h"
#include "table.h"
#include "data.h"

void load_bd_data(string path_heavy, string path_light, vector<dsbcell>& cells) {
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
	cout << "sucessfully loaded " << cells.size() << " double stranded cells from files " << path_heavy << ' ' << path_light << endl;
}

void load_bd_data(string path, vector<ssbcell>& cells) {
	itablestream itable(path);
	tablerow row(&itable);
	
	while (!row.eof) {
		if (row.get("AA CDR3") != "[CDR3_not_canonical]") {
			cells.emplace_back(row.get("Cell Label"), bcell_chain(row.get("V"), row.get("AA CDR3")));
		}
		row = tablerow(&itable);
	}
	cout << "sucessfully loaded " << cells.size() << " single stranded cells from file " << path << endl;
}


void load_10x_data(string path, vector<dsbcell>& cells) {
	itablestream itable(path);
	tablerow row(&itable);
	while (!row.eof) {
		string id = row.get("barcode");
		int match = -1;
		for (int i = 0; i < cells.size(); i ++) {
			if (cells[i].id == id) {
				match = i;
			}
		}
		string chain = row.get("chain");
		if (match == -1) {
			if (chain == "IGL" or chain == "IGK") {
				cells.emplace_back(id, bcell_chain("","",""), bcell_chain(row.get("v_gene"), row.get("cdr3")));
			} else if (chain == "IGH") {
				cells.emplace_back(id, bcell_chain(row.get("v_gene"), row.get("cdr3")), bcell_chain("","",""));
			}
		} else {
			if (chain == "IGL" or chain == "IGK") {
				if (!cells[match].light.valid) {
					cells[match].light = bcell_chain(row.get("v_gene"), row.get("cdr3"));
				}
			} else if (chain == "IGH") {
				if (!cells[match].heavy.valid) {
					cells[match].heavy = bcell_chain(row.get("v_gene"), row.get("cdr3"));
				}
			}
		}
		row = tablerow(&itable);
	}
	
	for (int i = cells.size()-1; i >= 0; i --) {
		if (!cells[i].heavy.valid or !cells[i].light.valid) {
			cells.erase(cells.begin() + i);
		}
	}
	cout << "sucessfully loaded in " << cells.size() << " cells from file " << path << endl;
}


void load_dekosky_data(string path, vector<dsbcell>& cells) {
	itablestream itable(path);
	tablerow row(&itable);
	int id = 0;
	while (!row.eof) {
		char heavy_light[] = {'H','L'};
		string aaseqs[2];
		int i = 0;
		for (char hl : heavy_light) {
			string nucseq = row.get(hl + string("3 Junction"));
			for (int j = 0; j < nucseq.length(); j += 3) {
				codon c {toupper(nucseq[j]), toupper(nucseq[j+1]), toupper(nucseq[j+2])};
				char aa = nuc_to_aa[c];
				aaseqs[i].push_back(aa);
			}
			i++;
		}
		string vhgene = row.get("VH Gene");
		string vlgene = row.get("VL Gene");
		int reads = atoi(row.get("Read Count").c_str());
		if (reads > 2) {
			cells.emplace_back(std::to_string(id), bcell_chain(vhgene, aaseqs[0]), bcell_chain(vlgene, aaseqs[1]));
		}
		row = tablerow(&itable);
		id ++;
	}
	cout << "sucessfully loaded " << cells.size() << " cells from " << path << endl;
}
			



void save_dist_matrix(string path, vector<dsbcell>& cells) {
	vector<string> ids {"cell-id"};
	for (dsbcell& cell : cells) {
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


void save_dist_matrix(string path, vector<ssbcell>& cells) {
	vector<string> ids {"cell-id"};
	for (ssbcell& cell : cells) {
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
