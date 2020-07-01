#include "../cell.h"

int main() {
	
	string seq1 = "ABAsdlfjsld";
	string seq2 = "sdlXXA";
	
	cout << seq1 << endl << seq2 << endl;
	
	bcell::align_aa(seq1, seq2);
	
	cout << seq1 << endl << seq2 << endl;
}
