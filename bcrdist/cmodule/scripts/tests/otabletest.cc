#include "../table.h"



int main() {
	vector<string> headers {"id", "rand", "num"};
	otablestream otable("tmp.csv", &headers);
	for (int i = 0; i < 1000; i ++) {
		tablerow row;
		row.add("id", i);
		row.add("rand", rand());
		row.add("num", i*i);
		otable.writeline(&row);
	}
}
