#ifndef CELL_CC
#define CELL_CC

#include "cell.h"
#include "data.h"

bcell::bcell(string newid, string v_gene, string newcdr3): id(newid), cdr3(newcdr3) {
	//cout << newid << ' ' << v_gene << ' ' << newcdr3 << endl;
	stringstream vgeness(v_gene);
	while (!vgeness.eof() and vgenes_to_cdrs.find(v_gene) == vgenes_to_cdrs.end()) {
		getline(vgeness, v_gene, ';');
	}
	//cout << v_gene << endl;
	pair<string,string> cdrs = vgenes_to_cdrs[v_gene];
	cdr1 = cdrs.first;
	cdr2 = cdrs.second;
	//cout << cdr1 << ' ' << cdr2 << endl;
}

bcell::bcell(string newid, string newcdr1, string newcdr2, string newcdr3): id(newid), cdr1(newcdr1), cdr2(newcdr2), cdr3(newcdr3) {
	
}

double bcell::aadist(string& seq1, string& seq2) {
	double distance = 0;
	for (int i = 0; i < seq1.length(); i ++) {
		distance += blosum_distances[pair<char,char>(seq1[i],seq2[i])];
	}
	//cout << "result " << seq1 << ' ' << seq2 << ' ' << distance << endl;
	return distance;
}

double bcell::unaligned_dist(string seq1, string seq2) {
	if (seq1.length() == seq2.length()) {
		return aadist(seq1, seq2);
	} else if (seq1.length() < seq2.length()) {
		align_aa(seq2, seq1);
		return aadist(seq1, seq2);
	} else {
		align_aa(seq1, seq2);
		return aadist(seq1, seq2);
	}
}

void bcell::align_aa(string& longer, string& shorter) {
	int num_gaps = longer.length() - shorter.length();
	
	vector<aa_match> matches;
	
	int start_off = 0;
	
	for (int i = 0; i < shorter.length(); i ++) {
		int off = start_off;
		while (shorter[i] != longer[i+off] and off < num_gaps) {
			off ++;
		}
		if (shorter[i] == longer[i+off]) {
			//cout << "doing it " << i << ' ' << off << ' ' << shorter[i] << ' ' << longer[i+off] << endl;
			aa_match match {i, i+off, 2};
			while (shorter.substr(match.seq1start, match.len) == longer.substr(match.seq2start, match.len) and match.seq2start + match.len <= longer.length()) {
				match.len ++;
			}
			match.len --;
			start_off = off;
			i += match.len-1;
			matches.push_back(match);
			cout << "match " << match.seq1start << ' ' << match.seq2start << ' ' << match.len << ' ' << shorter.substr(match.seq1start, match.len) << ' ' << off << endl;
		}
		// cout << start_off << endl;
		// for (int off = start_off + tmp_start_off; off <= num_gaps; off ++) {
		// 	if (shorter[i] == longer[i+off]) {
		// 		//cout << "doing it " << i << ' ' << off << ' ' << shorter[i] << ' ' << longer[i+off] << endl;
		// 		aa_match match {i, i+off, 2};
		// 		while (shorter.substr(match.seq1start, match.len) == longer.substr(match.seq2start, match.len) and match.seq2start + match.len <= longer.length()) {
		// 			match.len ++;
		// 		}
		// 		match.len --;
		// 		//i += match.len-1;
		// 		start_off = off;
		// 		tmp_start_off = match.len;
		// 		off += match.len-1;
		// 		bool exists = false;
		// 		matches.push_back(match);
		// 		cout << "match " << match.seq1start << ' ' << match.seq2start << ' ' << match.len << ' ' << shorter.substr(match.seq1start, match.len) << ' ' << off << endl;
		// 	}
		// }
	}
	
	int max_len = 0;
	
	for (int i = 0; i < matches.size(); i ++) {
		for (int j = i+1; j < matches.size(); j ++) {
			cout << "1 " << matches[i].seq1start << ' ' << matches[i].seq2start << ' ' << matches[i].len << ' ' << shorter.substr(matches[i].seq1start, matches[i].len) << endl;
			cout << "2 " << matches[j].seq1start << ' ' << matches[j].seq2start << ' ' << matches[j].len << ' ' << shorter.substr(matches[j].seq1start, matches[j].len) << endl;
			if (matches[i].seq2start - matches[i].seq1start > matches[j].seq2start - matches[j].seq1start) {
				cout << "problem!" << endl;
				if (matches[i].len < matches[j].len) {
					matches.erase(matches.begin()+i);
					cout << 1 << endl;
					i --;
					break;
				} else {
					cout << 2 << endl;
					matches.erase(matches.begin()+j);
					j --;
				}
			}
		}
	}
	
	int placed_gaps = 0;
	for (aa_match match : matches) {
		cout << "final match" << match.seq1start << ' ' << match.seq2start << ' ' << match.len << ' ' << shorter.substr(match.seq1start, match.len) << endl;
		int new_gaps = match.seq2start - match.seq1start - placed_gaps;
		for (int i = 0; i < new_gaps; i ++) {
			shorter.insert(shorter.begin() + match.seq1start + placed_gaps, '.');
		}
		placed_gaps += new_gaps;
	}
	
	for (int i = shorter.length(); i < longer.length(); i ++) {
		shorter.push_back('.');
	}
	
	cout << longer << endl << shorter << endl << endl;
}

double bcell::distance(bcell* other) {
	int distance = 0;
	distance += aadist(cdr1, other->cdr1) * dist_params::v_weight;
	distance += aadist(cdr2, other->cdr2) * dist_params::v_weight;
	distance += unaligned_dist(cdr3, other->cdr3) * dist_params::cdr3_weight;
	//cout << "total result " << distance << endl;;
	return distance;
}

#endif
