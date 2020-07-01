#ifndef FILEIO_H
#define FILEIO_H

#include "classes.h"

void load_bd_data(string path, vector<bcell>& cells);

void save_dist_matrix(vector<bcell> cells);

#endif
