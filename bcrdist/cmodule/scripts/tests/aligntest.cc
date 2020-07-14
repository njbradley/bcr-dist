#include "../cell.h"

int main(int numargs, const char** args) {
	
	vector<pair<string,string> > tests {
		{"XABC",
		 "ABC"},
		{"XAXBXC",
		 "ABCB"},
		{"XAXAXBCDXX",
		 "A_A__BCD"},
		{"VFVDVDVD",
		 "WEWEW"},
		{"HOWDIDWEDO",
		 "HWXDID"},
		{"XAXAXAXAXAXAXAXAXAXBCDXX",
		 "AGAGAGAGAGAGAGAGAGGBCD"}
	 };
	
	
	// string seq1 = args[1];
	// string seq2 = args[2];
	int start = clock();
	for (pair<string,string> seqs : tests) {
		string longer = seqs.first;
		string shorter = seqs.second;
		cout << longer << endl << shorter << endl;
		bcell_chain::align_aa(longer, shorter);
		cout << longer << endl << shorter << endl << endl;
	}
	cout << (clock()-start)/double(CLOCKS_PER_SEC) << endl;
	
}
