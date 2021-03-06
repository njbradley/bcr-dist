#ifndef DATA_H
#define DATA_H

#include "classes.h"

struct pair_hash {
	template <class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2> &pair) const {
		return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
	}
};

struct codon {
	char x;
	char y;
	char z;
	friend bool operator == (const codon& a, const codon& b) {
		return a.x == b.x and a.y == b.y and a.z == b.z;
	}
	struct hash {
		size_t operator () (const codon& c) const {
			std::hash<char> char_hash;
			return char_hash(c.x) ^ char_hash(c.y) ^ char_hash(c.z);
		}
	};
};


extern unordered_map<string,pair<string,string> > vgenes_to_cdrs;
extern unordered_map<string,string> vgenes_to_full_cdr;
extern unordered_map<pair<char,char>,double,pair_hash> blosum_distances;
extern vector<char> amino_acids;
extern unordered_map<codon,char,codon::hash> nuc_to_aa;

namespace dist_params {
	extern int gap_penalty;
	extern int cdr3_weight;
	extern int v_weight;
}

namespace cdr_params {
	extern int start_cdr1;
	extern int len_cdr1;
	extern int start_cdr2;
	extern int len_cdr2;
}

void load_persistant_data(string path = "");
void load_blosum(string path = "");

#endif
