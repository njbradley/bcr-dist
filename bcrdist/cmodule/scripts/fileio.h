#ifndef FILEIO_H
#define FILEIO_H

#include "classes.h"

void load_bd_data(string path_heavy, string path_light, vector<bcell_double>& cells);
void load_bd_data(string path, vector<bcell_single>& cells);

void load_10x_data(string path, vector<bcell_double>& cells);

void load_dekosky_data(string path, vector<bcell_double>& cells);

void save_dist_matrix(string path, vector<bcell_double>& cells);
void save_dist_matrix(string path, vector<bcell_single>& cells);


#endif
