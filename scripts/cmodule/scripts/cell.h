#ifndef CELL_H
#define CELL_H

#include "classes.h"

struct bcell {
	string id;
	string cdr1;
	string cdr2;
	string cdr3;
	bcell(string newid, string v_gene, string newcdr3);
	bcell(string newid, string newcdr1, string newcdr2, string newcdr3);
	double distance(bcell* other);
	static double aadist(string& seq1, string& seq2);
	static double unaligned_dist(string seq1, string seq2);
	static void align_aa(string& longer, string& shorter);
};

struct aa_match {
	int seq1start;
	int seq2start;
	int len;
};

#endif
