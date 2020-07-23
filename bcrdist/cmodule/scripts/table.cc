#ifndef TABLE_CC
#define TABLE_CC

#include "table.h"


tablerow::tablerow(): eof(false) {
	
}

tablerow::tablerow(itablestream* itable): eof(false) {
	itable->readline(this);
}

void tablerow::add(string header, string value) {
	items[header] = value;
}

void tablerow::add(string header, int value) {
	items[header] = to_string(value);
}

void tablerow::add(string header, double value) {
	if (int(value) == value) {
		items[header] = to_string(int(value));
	} else {
		items[header] = to_string(int(value));
	}
}

string tablerow::get(string header) {
	return items[header];
}


char get_delim(string path) {
	stringstream ss(path);
	string ext;
	while (!ss.eof()) {
		getline(ss, ext, '.');
	}
	if (ext == "csv") {
		return ',';
	} else if (ext == "tsv") {
		return '\t';
	}
	cout << "unrecognised file extension, defaulted to csv" << endl;
	return ',';
}

itablestream::itablestream(string path): ifile(path), delim(get_delim(path)) {
	if (!ifile.good()) {
		cout << "ERR: file " << path << " does not exist" << endl;
		exit(1);
	}
	string line_str;
	getline(ifile, line_str);
	stringstream line(line_str);
	string word;
	do {
		getline(line, word, delim);
		if (word != "") {
			headers.push_back(word);
		}
	} while (!line.eof());
}

void itablestream::readline(tablerow* row) {
	string line_str;
	getline(ifile, line_str);
	if (line_str == "") {
		row->eof = true;
	} else {
		stringstream line(line_str);
		string word;
		for (int i = 0; i < headers.size(); i ++) {
			getline(line, word, delim);
			row->add(headers[i], word);
		}
	}
}


otablestream::otablestream(string path, vector<string>* newheaders): ofile(path), headers(*newheaders), delim(get_delim(path)) {
	for (string head : headers) {
		ofile << head << ',';
	}
	ofile << endl;
}

void otablestream::writeline(tablerow* row) {
	for (string head : headers) {
		ofile << row->get(head) << delim;
	}
	ofile << endl;
}


#endif
