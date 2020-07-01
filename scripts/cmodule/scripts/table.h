#ifndef TABLE_H
#define TABLE_H

#include "classes.h"

struct tablerow {
	bool eof;
	unordered_map<string,string> items;
	tablerow();
	tablerow(itablestream* itable);
	void add(string header, string value);
	void add(string header, int value);
	void add(string header, double value);
	string get(string header);
};

struct itablestream {
	ifstream ifile;
	vector<string> headers;
	char delim;
	itablestream(string path);
	void readline(tablerow* row);
};

struct otablestream {
	ofstream ofile;
	vector<string> headers;
	char delim;
	otablestream(string path, vector<string>* newheaders);
	void writeline(tablerow* row);
};
	

#endif
