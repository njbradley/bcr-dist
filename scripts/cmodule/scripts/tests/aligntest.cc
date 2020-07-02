#include "../cell.h"

int main(int numargs, const char** args) {
	
	vector<pair<string,string> > {
		{"XABC",
		 "ABC"},
		{"XAXBXC",
		 "ABCB"},
		{"XAXAXAXBCDXX",
		 "A_A_A__BCD"},
		{"VFVDVDVD",
		 "WEWEW"},
		{
	
	
	string seq1 = "sjkfksdlhfajkshdlfka";
	string seq2 = "HJGJK";
	
	// string seq1 = args[1];
	// string seq2 = args[2];
	
	cout << seq1 << endl << seq2 << endl;
	
	int start = clock();
	for (int i = 0; i < 1000; i ++) {
		string longer = seq1;
		string shorter = seq2;
		bcell::align_aa(longer, shorter);
	}
	cout << (clock()-start)/double(CLOCKS_PER_SEC) << endl;
	
	cout << seq1 << endl << seq2 << endl;
}
