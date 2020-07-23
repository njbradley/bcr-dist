#ifndef DATA_CC
#define DATA_CC

#include "data.h"
#include "table.h"

unordered_map<string,pair<string,string> > vgenes_to_cdrs;
unordered_map<string,string> vgenes_to_full_cdr;
unordered_map<pair<char,char>,double,pair_hash> blosum_distances;
vector<char> amino_acids;
unordered_map<codon,char,codon::hash> nuc_to_aa;

namespace dist_params {
	int gap_penalty = 4;
	int cdr3_weight = 3;
	int v_weight = 1;
}


namespace cdr_params {
	int start_cdr1 = 26;
	int len_cdr1 = 12;
	int start_cdr2 = 55;
	int len_cdr2 = 10;
}


void load_persistant_data(string path) {
	itablestream vgenetable(path + "data/ig_all_genes.csv");
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
	
	ifstream codonfile(path + "data/codons.txt");
	char x,y,z,aa;
	while(codonfile >> x >> y >> z >> aa) {
		nuc_to_aa[{x,y,z}] = aa;
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
	string linestr;
	getline(ifile, linestr);
	stringstream line(linestr);
	line >> lett;
	while (!line.eof()) {
		amino_acids.push_back(lett);
		line >> lett;
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


#endif
