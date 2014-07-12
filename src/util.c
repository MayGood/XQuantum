/*
 * util.c
 *
 *  Created on: Mar 31, 2014
 *      Author: xxhx
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "util.h"

#define TIME_FORMAT "%b %d, %Y - %H:%M:%S"
int time_str(int max, char *str){
	time_t t;
	struct tm *tmp;

	t = time(NULL);
	tmp = localtime(&t);
	if (tmp == NULL){
		perror("localtime");
		exit(EXIT_FAILURE);
	}

	return strftime(str, max, TIME_FORMAT, tmp);
}

int findf(FILE *fd, int n, char *str, ...){

	int  s;
	int  ret;
	char buf[1024];

	for (s = 0; (ret = fscanf(fd, "%s", buf)) == 1;){

		if (strncmp(buf, *(&str + s), 1024) == 0 ||
			strcmp("*", *(&str + s)) == 0) s++; else s = 0;

		if (s == n) break;
	}

	return ret;
}

struct AtomName_t{
	int   Z;
	char  *shortName;
	char  *longName;
};

#define MAX_PERIODIC_ATOM 103
static struct AtomName_t PeriodicName[MAX_PERIODIC_ATOM] =
{
	{1, "H", "HYDROGEN"},  \
	{2, "HE", "HELIUM"}, \
	{3, "LI","LITHIUM"},   \
	{4, "BE","BERYLLIUM"}, \
	{5, "B", "BORON"}, \
	{6, "C", "CARBON"}, \
	{7, "N", "NITROGEN"}, \
	{8, "O", "OXYGEN"}, \
	{9, "F", "FLUORINE"}, \
	{10, "NE", "NEON"}, \
	{11, "NA", "SODIUM"}, \
	{12, "MG", "MAGNESIUM"}, \
	{13, "AL", "ALUMINUM"}, \
	{14, "SI", "SILICON"}, \
	{15, "P", "PHOSPHOROUS"}, \
	{16, "S", "SULFUR"}, \
	{17, "CL", "CHLORINE"}, \
	{18, "AR", "ARGON"}, \
	{19, "K", "POTASSIUM"}, \
	{20, "CA", "CALCIUM"}, \
	{21, "SC", "SCANDIUM"}, \
	{22, "TI", "TITANIUM"}, \
	{23, "V", "VANADIUM"}, \
	{24, "CR", "CHROMIUM"}, \
	{25, "MN", "MANGANESE"}, \
	{26, "FE", "IRON"}, \
	{27, "CO", "COBALT"}, \
	{28, "NI", "NICKEL"}, \
	{29, "CU", "COPPER"}, \
	{30, "ZN", "ZINC"}, \
	{31, "GA", "GALLIUM"}, \
	{32, "GE", "GERMANIUM"}, \
	{33, "AS", "ARSENIC"}, \
	{34, "SE", "SELENIUM"}, \
	{35, "BR", "BROMINE"}, \
	{36, "KR", "KRYPTON"}, \
	{37, "RB", "RUBIDIUM"}, \
	{38, "SR", "STRONTIUM"}, \
	{39, "Y", "YTTRIUM"}, \
	{40, "ZR", "ZIRCONIUM"}, \
	{41, "NB", "NIOBIUM"}, \
	{42, "MO", "MOLYBDENUM"}, \
	{43, "TC", "TECHNETIUM"}, \
	{44, "RU", "RUTHENIUM"}, \
	{45, "RH", "RHODIUM"}, \
	{46, "PD", "PALLADIUM"}, \
	{47, "AG", "SILVER"}, \
	{48, "CD", "CADMIUM"}, \
	{49, "IN", "INDIUM"}, \
	{50, "SN", "TIN"}, \
	{51, "SB", "ANTIMONY"}, \
	{52, "TE", "TELLURIUM"}, \
	{53, "I", "IODINE"}, \
	{54, "XE", "XENON"}, \
	{55, "CS", "CESIUM"}, \
	{56, "BA", "BARIUM"}, \
	{57, "LA", "LANTHANUM"}, \
	{58, "CE", "CERIUM"}, \
	{59, "PR", "PRASEODYMIUM"}, \
	{60, "ND", "NEODYMIUM"}, \
	{61, "PM", "PROMETHIUM"}, \
	{62, "SM", "SAMARIUM"}, \
	{63, "EU", "EUROPIUM"}, \
	{64, "GD", "GADOLIUM"}, \
	{65, "TB", "TERBIUM"}, \
	{66, "DY", "DYSPROSIUM"}, \
	{67, "HO", "HOLMIUM"}, \
	{68, "ER", "ERBIUM"}, \
	{69, "TM", "THULIUM"}, \
	{70, "YB", "YTTERBIUM"}, \
	{71, "LU", "LUTETIUM"}, \
	{72, "HF", "HAFNIUM"}, \
	{73, "TA", "TANTALUM"}, \
	{74, "W", "TUNGSTEN"}, \
	{75, "RE", "RHENIUM"}, \
	{76, "OS", "OSMIUM"}, \
	{77, "IR", "IRIDIUM"}, \
	{78, "PT", "PLATINUM"}, \
	{79, "AU", "GOLD"}, \
	{80, "HG", "MERCURY"}, \
	{81, "TL", "THALLIUM"}, \
	{82, "PB", "LEAD"}, \
	{83, "BI", "BISMUTH"}, \
	{84, "PO", "POLONIUM"}, \
	{85, "AT", "ASTATINE"}, \
	{86, "RN", "RADON"}, \
	{87, "FR", "FRANCIUM"}, \
	{88, "RA", "RADIUM"}, \
	{89, "AC", "ACTINIUM"}, \
	{90, "TH", "THORIUM"}, \
	{91, "PA", "PROTACTINIUM"}, \
	{92, "U", "URANIUM"}, \
	{93, "NP", "NEPTUNIUM"}, \
	{94, "PU", "PLUTONIUM"}, \
	{95, "AM", "AMERICIUM"}, \
	{96, "CM", "CURIUM"}, \
	{97, "BK", "BERKELIUM"}, \
	{98, "CF", "CALIFORNIUM"}, \
	{99, "ES", "EINSTEINIUM"}, \
	{100, "FM", "FERMIUM"}, \
	{101, "MD", "MENDELEVIUM"}, \
	{102, "NO", "NOBELIUM"}, \
	{103, "LR", "LAWRENCIUM"}
};

#ifndef SYMB_SHORTNAME
#define SYMB_SHORTNAME 0
#endif
#ifndef SYMB_LONGNAME
#define SYMB_LONGNAME  1
#endif
int sym2Z(const char *sym, int type){
	char str[1024];
	int i;

	// make sym upper case
	for (i = 0; i < strlen(sym); i++)
		str[i] = toupper(sym[i]);
	// null termination
	str[i] = '\0';

	for (i = 0; i < MAX_PERIODIC_ATOM; i++){
		switch (type){
		case SYMB_SHORTNAME:
			if (strcmp(str, PeriodicName[i].shortName) == 0)
				return PeriodicName[i].Z;
			break;
		case SYMB_LONGNAME:
			if (strcmp(str, PeriodicName[i].longName) == 0)
				return PeriodicName[i].Z;
			break;
		}
	}
	return 0;
}

void Z2Sym(int Z, char *sym){
	if (0 < Z && Z <= MAX_PERIODIC_ATOM)
		strcpy(sym, PeriodicName[Z - 1].longName);
	else{
		printf("Z2sym - Error : Z = %d\n", Z);
		exit(EXIT_FAILURE);
	}
}

void Z2SymShort(int Z, char *sym){
	if (0 < Z && Z <= MAX_PERIODIC_ATOM)
		strcpy(sym, PeriodicName[Z - 1].shortName);
	else{
		printf("Z2sym - Error : Z = %d\n", Z);
		exit(EXIT_FAILURE);
	}
}
