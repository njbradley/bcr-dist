#ifndef FILEIO_C
#define FILEIO_C

#include "fileio.h"
#include "cell.h"
#include "table.h"
#include "data.h"

bool load_bd_data(string path_heavy, string path_light, vector<dsbcell>& cells) {
	itablestream itable_heavy(path_heavy);
	itablestream itable_light(path_light);
	
	if (!itable_heavy.good() or !itable_light.good()) {
		return false;
	}
	
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
	return true;
}

bool load_bd_data(string path, vector<dsbcell>& cells) {
	itablestream itable(path);
	
	if (!itable.good()) {
		return false;
	}
	
	tablerow row(&itable);
	
	while (!row.eof) {
		bcell_chain heavy(row.get("BCR_Heavy_V_gene_Dominant"), row.get("BCR_Heavy_CDR3_Translation_Dominant"));
		bcell_chain light(row.get("BCR_Light_V_gene_Dominant"), row.get("BCR_Light_CDR3_Translation_Dominant"));
		if (heavy.valid and light.valid) {
			cells.emplace_back(row.get("Cell_Index"), heavy, light);
		} else {
			// for (pair<string,string> val : row.items) {
			// 	cout << val.first << ',' << val.second << ' ';
			// }
			//cout << row.get("Cell_Index") << endl;
			//cout << row.get("BCR_Heavy_V_gene_Dominant") << ' ' << row.get("BCR_Light_V_gene_Dominant") << endl;
			//cout << row.get("BCR_Heavy_CDR3_Translation_Dominant") << ' ' << row.get("BCR_Light_CDR3_Translation_Dominant") << endl;
		}
		row = tablerow(&itable);
	}
	cout << "sucessfully loaded " << cells.size() << " double stranded cells from file " << path << endl;
	return true;
}

bcell_chain get_10x_chain(tablerow* row) {
	if (row->get("cdr1") != "") {
		return bcell_chain(row->get("cdr1"), row->get("cdr2"), row->get("cdr3"));
	} else {
		return bcell_chain(row->get("v_gene"), row->get("cdr3"));
	}
}

bool load_10x_data(string path, vector<dsbcell>& cells) {
	itablestream itable(path);
	
	if (!itable.good()) {
		return false;
	}
	
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
				cells.emplace_back(id, bcell_chain(), get_10x_chain(&row), row.get("raw_clonotype_id"));
			} else if (chain == "IGH") {
				cells.emplace_back(id, get_10x_chain(&row), bcell_chain(), row.get("raw_clonotype_id"));
			}
		} else {
			if (chain == "IGL" or chain == "IGK") {
				if (!cells[match].light.valid) {
					cells[match].light = get_10x_chain(&row);
				}
			} else if (chain == "IGH") {
				if (!cells[match].heavy.valid) {
					cells[match].heavy = get_10x_chain(&row);
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
	return true;
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
				codon c {char(toupper(nucseq[j])), char(toupper(nucseq[j+1])), char(toupper(nucseq[j+2]))};
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
			


void read_dist_matrix(string path, float* matrix, int dim) {
	ifstream ifile(path);
	string buff;
	getline(ifile, buff);
	for (int i = 0; i < dim; i ++) {
		ifile >> buff;
		for (int j = 0; j < i; j ++) {
			float val = -2;
			ifile >> val;
			matrix[i*dim + j] = val;
		}
	}
	
	for (int i = 0; i < dim; i ++) {
		matrix[i*dim + i] = 0;
	}
	
	for (int i = 0; i < dim; i ++) {
		for (int j = i+1; j < dim; j ++) {
			matrix[i*dim + j] = matrix[j*dim + i];
		}
	}
}



#endif
