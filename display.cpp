#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>

void saveGif(char* endf, char* endr,double dx, double dt,double tick) {
	char gnuplotparam[200];
	sprintf(gnuplotparam, "start gnuplot  -e dt=%f -e dx=%f -e tick=%f -e endf=\'%s\' -e endr=\'%s\' -p  \"scripts\\salvarGif.txt\" ", dt, dx,tick, endf,endr);
	system(gnuplotparam);
}

void exibirGif(char* endf, double dx, double dt) {
	char gnuplotparam[200];
	sprintf(gnuplotparam, "gnuplot  -e dt=%f -e dx=%f -e end=\'%s\' -p  \"scripts\\exibirGif.txt\" ", dt, dx, endf);
	system(gnuplotparam);
}

void saveFoto(char* endf, char* endr, double dx, double dt,int tick) {
	char gnuplotparam[200];
	sprintf(gnuplotparam, "gnuplot  -e dt=%f -e dx=%f -e endf=\'%s\' -e endr=\'%s\' -e tick=%d -p  \"scripts\\salvarFoto.txt\" ", dt, dx, endf, endr,tick);
	system(gnuplotparam);

}
void exibirFoto(char* endf, double dx, double dt,double dtTarget) {
	char gnuplotparam[200];
	sprintf(gnuplotparam, "gnuplot  -e dt=%f -e dx=%f -e end=\'%s\' -e dtT=\'%d\'  -p  \"scripts\\exibirFoto.txt\" ", dt, dx, endf,dtTarget);
	system(gnuplotparam);
}
void saveDl(char* end, char* endf,std::string nome1, std::string nome2) {
	char gnuplotparam[200];
	sprintf(gnuplotparam, "gnuplot  -e end=\'%s\'  -e endf=\'%s\' -e n1=\'%s\' -e n2=\'%s\'  -p  \"scripts\\saveDoubleLine.txt\" ", end, endf,nome1,nome2);
	system(gnuplotparam);
}
