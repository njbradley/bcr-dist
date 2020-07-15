#ifndef CELL_H
#define CELL_H

#include "classes.h"

struct bcell_chain {
	string cdr1;
	string cdr2;
	string cdr3;
	bcell_chain(string v_gene, string newcdr3);
	bcell_chain(string newcdr1, string newcdr2, string newcdr3);
	bcell_chain(istream& ifile);
	void to_file(ostream& ofile);
	double distance(bcell_chain* other);
	static double aadist(string& seq1, string& seq2);
	static double unaligned_dist(string seq1, string seq2);
	static void align_aa(string& longer, string& shorter);
};

struct bcell_single {
	string id;
	bcell_chain chain;
	double distance(bcell_single* other);
	bcell_single(string newid, bcell_chain newchain);
	bcell_single(istream& ifile);
	void to_file(ostream& ofile);
};

struct bcell_double {
	string id;
	bcell_chain heavy;
	bcell_chain light;
	double distance(bcell_double* other);
	bcell_double(string newid, bcell_chain newheavy, bcell_chain newlight);
	bcell_double(istream& ifile);
	void to_file(ostream& ofile);
};

struct aa_match {
	int seq1start;
	int seq2start;
	int len;
};

#endif
