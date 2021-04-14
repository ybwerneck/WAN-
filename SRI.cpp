// SRI.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <thread>         // std::thread
#include "utils.h"// std::thread
#include "matrix.h"
#include <condition_variable>
#include "display.h"
#include <omp.h>
#include <math.h>       /* sin */
#include <chrono>
using namespace std::chrono_literals;
#include <chrono>
#include <cstdlib>

Matrix* S, * I;

//--------Parametros simulacao
int L = 60;
double dt = (1.0 / 40.0), dx = 1.0 / 3.0;
//possivel problema com long double
int tam = L/dx;
//--------- Parametros modelo
Matrix *D11 ,*D22;
long double R0 = 1.17 , Rd= 1.3 , v =0.15;
long double h = 1.0 / 3.0;
long double v0 = 0.01;
long double u0 = 0.005;
//---------
FILE* Sarq, * Iarq, * Integralsarq, * D1arq,*D2arq,*Rarq;
std::string subpasta = "resultado";

std::string filename0 = "data/Int.txt";
std::string filename = "data/S.txt";
std::string filename2 = "data/I.txt";
std::string filename3 = "data/D1.txt";
std::string filename5 = "data/D2.txt";
std::string filename4 = "data/R.txt";

int n = 200000;
int random(int min, int max) {
    return rand() % (max + 1 - min) + min;

}
void initValues() {
  
    S = new Matrix(tam, tam,0);
    I = new Matrix(tam, tam, 0);
    
    for (int i = 0; i < tam; i++){
        for (int j = 0; j < tam; j++)
        {
            long double r1 = random(80,100)*0.01;
            long double r2 = random(80, 100) * 0.01;
            S ->operator()(i, j,r1*v0);
            I ->operator()(i,j,r2*u0);

        }
    }

    double centrox = tam / 2, centroy = tam / 2;
    double maxdistancia= pow(pow(centrox - 0, 2) + pow(centroy - 0, 2), 0.5);

    D11 = new Matrix(tam, tam, 0.01);
    D22 = new Matrix(tam, tam, 0.25);
    for (int i = 0; i < tam; i++) {

        for (int j = 0; j < tam; j++)
        {
            double distancia = pow(pow(centrox-i,2)  + pow(centroy - j,2) ,0.5);
            double frac = (distancia / maxdistancia);
            //D11 ->operator()(i, j, 0.05 - pow(frac, 0.5) * 0.04);
            // D22 ->operator()(i, j, 1.0 - pow(frac, 0.5) * 0.9);

        }
    }


       


}
long double f(long double S, long double I) {


    long double t = v * Rd * (S + I) * (1.0 - (S + I)) - R0 * ((S*I) / (S + I)) - v * S;
    return t;
}
long double g(long double S, long double I) {

   long double t=R0*((S*I)/(S+I)) - I ;    //  1.212 * ((S*X)/(S+X)) - X ,S =0.003630726111,X = 0.0007697139352
   //Problema numerico
   return t; //t deveria ser igual a 0 de acordo com resultado da calculadora
}
long double difussion(Matrix* ua, int x, int y, long double coef) {

    Matrix* S = ua;
    long double difx, dify;
    long double Y = coef/(dx*dx);

    difx = (-2 * S->operator()(x, y) + S->operator()(x - 1, y) + S->operator()(x + 1, y)) * Y;

    dify = (-2 * S->operator()(x, y) + S->operator()(x, y - 1) + S->operator()(x, y + 1)) * Y;

    return difx + dify;

}
void attD(int k) {

}
void attDr(int k) {
}
void step() {


    int i, j;
    #pragma omp parallel for
        for (i = 1; i < tam - 1; i++) {
        for (j = 1; j < tam - 1; j++) {

            double d11 = D11->operator()(i, j);
            double d22 = D22->operator()(i, j);

            long double r1 = f(S->operator()(i, j), I->operator()(i, j));
            long double r2 = g(S->operator()(i, j), I->operator()(i, j));
            long double deltaS = (1/40.0) *(r1 + difussion(S, i, j, d11));
            long double deltaI = (1/40.0) * (r2 + difussion(I, i, j, d22) );
            S->operator()(i, j, deltaS + S->operator()(i, j));
            I->operator()(i, j, deltaI + I->operator()(i, j));
        }
    }

    //condição de borda de neuman
#pragma omp parallel for 
    for (int i = 0; i < tam; i++)
    {
        S->operator()(i, 0, S->operator()(i, 1));
        I->operator()(i, 0, I->operator()(i, 1));

        S->operator()(i, tam - 1, S->operator()(i, tam - 2));
        I->operator()(i, tam - 1, I->operator()(i, tam - 2));
    }
#pragma omp parallel for 
    for (int j = 0; j < tam; j++)
    {
        S->operator()(0, j, S->operator()(1, j));
        I->operator()(0, j, I->operator()(1, j));

        S->operator()(tam - 1, j, S->operator()(tam - 2, j));
        I->operator()(tam - 1, j, I->operator()(tam - 2, j));
    }


}
void printIntegrals(Matrix* S, Matrix* I, FILE* Iarq,double t){
           
       fprintf(Iarq, "%f %f %f",S->sum(),I->sum(),t);      
       fprintf(Iarq, "\n");

    };
void printCoeficients(FILE* Darq,FILE* Rarq, double t) {

    

    fprintf(Rarq, "%f %f %f", R0, Rd, t);
    fprintf(Rarq, "\n");
}; 
void printMatrixtoFile(Matrix* S, Matrix* I,  FILE* Uarq, FILE* Varq) {

    for (int v = 0; v < tam; v++)
    {
        int dx = 1;
        for (int j = 0; j < tam; j++)
        {

            fprintf(Uarq, "%f ", S->operator()(v, j)); //att
            fprintf(Varq, "%f ", I->operator()(v, j)); //att


        }
        fprintf(Uarq, "\n");
        fprintf(Varq, "\n");

     };




}
void makeResultFiles(std::string sub) {

    subpasta = sub;
    filename0 = subpasta +"/"+ filename0;
    filename = subpasta + "/" + filename;
    filename2 = subpasta + "/" + filename2;
    filename3 = subpasta + "/" + filename3;
    filename4 = subpasta + "/" + filename4;
    filename5 = subpasta + "/" + filename5;

    char param[200];
   
    sprintf(param, "mkdir %s",subpasta);
    system(param);
    sprintf(param, "mkdir %s\\data", subpasta);
    system(param);

    sprintf(param, "mkdir %s\\result", subpasta);
    system(param);

    D1arq = fopen(filename3.c_str(), "w");
    D2arq = fopen(filename5.c_str(), "w");

    Sarq = fopen(filename.c_str(), "w");
    Iarq = fopen(filename2.c_str(), "w");
    Integralsarq = fopen(filename0.c_str(), "w");
    Rarq = fopen(filename4.c_str(), "w");


}
std::atomic<bool> report = false;
void reporta(int * k, int parada) {
    while (report) {
        double o = parada == -1 ? n :parada;
        system("cls");
        printf("\nDT: %.10f DX:%.10f TAM= %d \n Quadro %d \ %d   completed %.2f%% \n", dt, dx, tam, n, *k, (*k/o) * 100.0);

        std::this_thread::sleep_for(300ms);
    }
}
int kt=0;
int main()
{
    initValues();
    
    makeResultFiles("heterogenea");
    report = true;
    int frames = 100; // quantos quadros terão na animação resultante, afeta muito o desempenho do gnuplot 

    double tick = (1.0 /frames) * n;
    int parada = -1;
    std::thread printer = std::thread(reporta, &kt, parada);

    for (int i = 0; i < n; i++) {
        
        //passo no tempo
        step();
        attD(i);
        attDr(i);
        kt = i;

        //rotinas de impressão no arquivo 
        if(i%40==0)
        printIntegrals(S, I, Integralsarq,  i);
        printCoeficients(D1arq,Rarq,i);

        if (parada==-1&&(i % (n/frames) == 0))
        {
            printMatrixtoFile(D11, D22, D1arq, D2arq);

            printMatrixtoFile(S, I,Sarq, Iarq);

        if (i != n - 1)
        {
            fprintf(Sarq, "\n\n");
            fprintf(Iarq, "\n\n");
            fprintf(D1arq, "\n\n");
            fprintf(D2arq, "\n\n");
         }
       }
        if (parada != -1 && i == parada)
        {

            printMatrixtoFile(S, I, Sarq, Iarq);


        }

    }

    report = false;
    fclose(D1arq);
    fclose(D2arq);
    fclose(Sarq);
    fclose(Iarq);
    fclose(Integralsarq);  


    printer.join(); 
    if (parada == -1) {
        saveGif((char*) filename.c_str(), (char*)(subpasta + "/result/S.gif").c_str(), dx, dt, tick);
        saveGif((char*) filename2.c_str(), (char*)(subpasta + "/result/I.gif").c_str(), dx, dt, tick);
        saveGif((char*)filename3.c_str(), (char*)(subpasta + "/result/D1.gif").c_str(), dx, dt, tick);
        saveGif((char*)filename5.c_str(), (char*)(subpasta + "/result/D2.gif").c_str(), dx, dt, tick);

    }
    else
    {
        saveFoto((char*) filename.c_str(), (char*)(subpasta+"/result/S.png").c_str(), dx,dt, parada);
        saveFoto((char*) filename2.c_str(), (char*)(subpasta + "/result/I.png").c_str(), dx, dt, parada);

    }
    // saveDl((char*)filename3.c_str(), (char*)(subpasta + "/result/Dif.png").c_str(), "d1", "d2");
    
    saveDl((char*)filename0.c_str(), (char*)(subpasta + "/result/Int.png").c_str(), "S", "I");
    saveDl((char*)filename4.c_str(), (char*)(subpasta + "/result/R.png").c_str(), "R0", "Rd");


}

