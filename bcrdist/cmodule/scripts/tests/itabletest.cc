#include "../table.h"

int main() {
	itablestream itable("tmp.csv");
	tablerow row(&itable);
	while (!row.eof) {
		for (string head : itable.headers) {
			cout << head << ':' << row.get(head) << endl;
		}
		row = tablerow(&itable);
	}
}
