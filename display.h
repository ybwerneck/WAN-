#pragma once
#ifndef DISPLAY_H
#define DISPLAY_H
void saveGif(char* endf, char* endr, double dx, double dt,double tick);
void saveDl(char* endf, char* end,std::string,std::string);
void saveTl(char* endf, char* end, std::string, std::string,std::string);
void exibirGif(char* endr, double dx, double dt);
void saveFoto(char* endf, char* endr, double dx, double dt,int tick);
void exibirFoto(char* endf, double dx, double dt,double dtTarget);
#endif
