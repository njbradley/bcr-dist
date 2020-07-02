#include "../cell.h"

int main(int numargs, const char** args) {
	
	string seq1 = "GAGAGABCGGABC";
	string seq2 = "AXAXXABC";
	
	// string seq1 = args[1];
	// string seq2 = args[2];
	
	cout << seq1 << endl << seq2 << endl;
	
	bcell::align_aa(seq1, seq2);
	
	cout << seq1 << endl << seq2 << endl;
}
