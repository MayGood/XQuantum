/*
 * util.h
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */

#ifndef UTIL_H_
#define UTIL_H_

#define BOHR2ANGSTROM 0.529177249000
#define ANGSTROM2BOHR 1.889725988579

#define HARTREE2EV    27.2113957

#define DEBYE2AU 0.393430307
#define AU2DEBYE 2.54174623

int time_str(int max, char *str);
int findf(FILE *fd, int n, char *str, ...);
#define SYMB_SHORTNAME 0
#define SYMB_LONGNAME  1
int sym2Z(const char *sym, int type);
void Z2Sym(int Z, char *sym);
void Z2SymShort(int Z, char *sym);

#endif /* UTIL_H_ */
