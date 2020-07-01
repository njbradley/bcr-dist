#ifndef DATA_H
#define DATA_H

#include "classes.h"

struct pair_hash {
	template <class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2> &pair) const {
		return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
	}
};

unordered_map<string,pair<string,string> > vgenes_to_cdrs;
unordered_map<string,string> vgenes_to_full_cdr;
unordered_map<pair<char,char>,double,pair_hash> blosum_distances;
vector<char> amino_acids;

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

void load_persistant_data(string path = "");
void load_blosum(string path = "");

#endif
