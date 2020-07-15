#ifndef CELL_CC
#define CELL_CC

#include "cell.h"
#include "data.h"

bcell_chain::bcell_chain(string v_gene, string newcdr3): cdr3(newcdr3) {
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

bcell_chain::bcell_chain(string newcdr1, string newcdr2, string newcdr3): cdr1(newcdr1), cdr2(newcdr2), cdr3(newcdr3) {
	
}

bcell_chain::bcell_chain(istream& ifile) {
	ifile >> cdr1 >> cdr2 >> cdr3;
}

void bcell_chain::to_file(ostream& ofile) {
	ofile << cdr1 << '\t' << cdr2 << '\t' << cdr3 << '\t';
}

double bcell_chain::aadist(string& seq1, string& seq2) {
	double distance = 0;
	for (int i = 0; i < seq1.length(); i ++) {
		distance += blosum_distances[pair<char,char>(seq1[i],seq2[i])];
	}
	//cout << "result " << seq1 << ' ' << seq2 << ' ' << distance << endl;
	return distance;
}

double bcell_chain::unaligned_dist(string seq1, string seq2) {
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


void bcell_chain::align_aa(string& longer, string& shorter) {
	int num_gaps = longer.length() - shorter.length();

	vector<aa_match> matches;

	int start_off = 0;

	for (int i = 0; i < shorter.length(); i ++) {
		aa_match final_match {0,0,0};
		int final_new_end = 0;
		int final_off = 0;
		int score = 0;
		
		for (int off = 0; off <= num_gaps; off ++) {
			if (shorter[i] == longer[i+off]) {
				aa_match match {i, i+off, 2};
				while (shorter.substr(match.seq1start, match.len) == longer.substr(match.seq2start, match.len) and match.seq2start + match.len <= longer.length()) {
					match.len ++;
				}
				match.len --;
				
				int deleted = 0;
				int size = matches.size();
				int new_end = matches.size()-1;
				int new_start_off = start_off;
				while (off < new_start_off and matches.size() > 0 and match.len > matches[new_end].len + deleted) {
					//cout << off << ' ' << new_start_off << ' ' << deleted << ' ' << new_end << endl;
					deleted += matches[new_end].len;
					new_end --;
					if (new_end == -1) {
						new_start_off = 0;
					} else {
						new_start_off = matches[new_end].seq2start - matches[new_end].seq1start;
					}
				}
				//cout << "end " << off << ' ' << new_start_off << ' ' << deleted << ' ' << new_end << endl;
				
				if (off >= new_start_off) {
					off += match.len-1;
					//cout << "match " << match.seq1start << ' ' << match.seq2start << ' ' << match.len << ' ' << shorter.substr(match.seq1start, match.len) << endl;
					int newscore = match.len - deleted;
					//cout << "score " << newscore << endl;
					if (newscore > score) {
						final_match = match;
						score = newscore;
						final_off = off;
						final_new_end = new_end;
					}
				}
			}
		}
		if (final_match.len > 0) {
			if (final_new_end+1 < matches.size()) {
				matches.erase(matches.begin() + final_new_end+1, matches.end());
			}
			start_off = final_off;
			i += final_match.len-1;
			matches.push_back(final_match);
			//cout << "good match " << final_match.seq1start << ' ' << final_match.seq2start << ' ' << final_match.len << ' ' << shorter.substr(final_match.seq1start, final_match.len) << endl;
		}
	}

	int placed_gaps = 0;
	for (aa_match match : matches) {
		int new_gaps = match.seq2start - match.seq1start - placed_gaps;
		for (int i = 0; i < new_gaps; i ++) {
			shorter.insert(shorter.begin() + match.seq1start + placed_gaps, '.');
		}
		placed_gaps += new_gaps;
	}

	for (int i = shorter.length(); i < longer.length(); i ++) {
		shorter.push_back('.');
	}
}

double bcell_chain::distance(bcell_chain* other) {
	int distance = 0;
	distance += aadist(cdr1, other->cdr1) * dist_params::v_weight;
	distance += aadist(cdr2, other->cdr2) * dist_params::v_weight;
	distance += unaligned_dist(cdr3, other->cdr3) * dist_params::cdr3_weight;
	return distance;
}




string get_id(istream& ifile) {
	string newid;
	ifile >> newid;
	return newid;
}


bcell_single::bcell_single(string newid, bcell_chain newchain): id(newid), chain(newchain) {
	
}

bcell_single::bcell_single(istream& ifile): id(get_id(ifile)), chain(ifile) {
	
}

void bcell_single::to_file(ostream& ofile) {
	ofile << id << '\t';
	chain.to_file(ofile);
}

double bcell_single::distance(bcell_single* other) {
	return chain.distance(&other->chain);
}




bcell_double::bcell_double(string newid, bcell_chain newheavy, bcell_chain newlight): id(newid), heavy(newheavy), light(newlight) {
	
}

bcell_double::bcell_double(istream& ifile): id(get_id(ifile)), heavy(ifile), light(ifile) {
	
}

void bcell_double::to_file(ostream& ofile) {
	ofile << id << '\t';
	heavy.to_file(ofile);
	light.to_file(ofile);
}

double bcell_double::distance(bcell_double* other) {
	return heavy.distance(&other->heavy) + light.distance(&other->light);
}

#endif
