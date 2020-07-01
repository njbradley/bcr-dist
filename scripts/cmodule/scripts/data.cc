#ifndef DATA_CC
#define DATA_CC

#include "data.h"
#include "table.h"

void load_persistant_data(string path) {
	itablestream vgenetable(path + "data/ighv_genes_final.csv");
	tablerow row(&vgenetable);
	while (!row.eof) {
		vgenes_to_full_cdr[row.get("vgene")] = row.get("fullseq");
		//cout << row.get("vgene") << ' ' << row.get("fullseq") << endl;
		row = tablerow(&vgenetable);
	}
	
	for (pair<string,string> kvpair : vgenes_to_full_cdr) {
		string vgene = kvpair.first;
		string seq = kvpair.second;
		string cdr1 = seq.substr(cdr_params::start_cdr1, cdr_params::len_cdr1);
		string cdr2 = seq.substr(cdr_params::start_cdr2, cdr_params::len_cdr2);
		//cout << vgene << ' ' << cdr1 << ' ' << cdr2 << endl;
		vgenes_to_cdrs[vgene] = pair<string,string>(cdr1,cdr2);
	}
	
	load_blosum(path);
}



void load_blosum(string path) {
	ifstream ifile(path + "data/blosum_scoring.txt");
	char lett;
	string buff;
	ifile >> lett;
	while (lett == '#') {
		getline(ifile, buff);
		ifile >> lett;
	}
	amino_acids.push_back(lett);
	for (int i = 0; i < 23; i ++) {
		ifile >> lett;
		amino_acids.push_back(lett);
	}
	
	char aa1;
	ifile >> aa1;
	while (!ifile.eof()) {
		for (char aa2 : amino_acids) {
			double score;
			double distance;
			ifile >> score;
			if (aa1 == aa2) {
				distance = 0;
			} else if (score < 0) {
				distance = 4;
			} else {
				distance = 4 - score;
			}
			blosum_distances[pair<char,char>(aa1, aa2)] = distance;
			//cout << aa1 << aa2 << ' ' << score << ' ' << distance << endl;
		}
		ifile >> aa1;
	}
}

int main() {
	load_persistant_data();
	return 0;
}

#endif
