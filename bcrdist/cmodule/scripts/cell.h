#ifndef CELL_H
#define CELL_H

#include "classes.h"

struct bcell {
	string id;
	string clonotype;
	bcell(string newid, string newclono);
	bcell(istream& ifile);
	virtual ~bcell();
	virtual double distance(bcell* other) = 0;
	virtual void to_file(ostream& ofile);
};

struct bcell_chain {
	string cdr1;
	string cdr2;
	string cdr3;
	bool valid;
	bcell_chain();
	bcell_chain(string v_gene, string newcdr3);
	bcell_chain(string newcdr1, string newcdr2, string newcdr3);
	bcell_chain(istream& ifile);
	void to_file(ostream& ofile);
	double distance(bcell_chain* other);
	static double aadist(string& seq1, string& seq2);
	static double unaligned_dist(string seq1, string seq2);
	static void align_aa(string& longer, string& shorter);
};

struct ssbcell: bcell {
	bcell_chain chain;
	double distance(bcell* other);
	ssbcell(string newid, bcell_chain newchain, string newclono = "-");
	ssbcell(istream& ifile);
	void to_file(ostream& ofile);
};

struct dsbcell: bcell {
	bcell_chain heavy;
	bcell_chain light;
	double distance(bcell* other);
	dsbcell(string newid, bcell_chain newheavy, bcell_chain newlight, string newclono = "-");
	dsbcell(istream& ifile);
	void to_file(ostream& ofile);
};

struct aa_match {
	int seq1start;
	int seq2start;
	int len;
};

#endif
