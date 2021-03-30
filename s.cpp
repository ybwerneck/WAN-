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
#include "../../../Documents/bolsa/adi/ADIFHN/display.h"

Matrix* S, * R, * I;

//--------Parametros simulacao
int L = 10;
int t = 2;
float dt = 0.01, dx = 0.1;

int tam = L / dx;
//--------- Parametros modelo
double d11 = 4, d12 = 10 , d22 = 0.5, d33 = 1;
double B = 0.01;
double mi = 0.01;
//---------

int n = t / dt;
float random(int min, int max) {
    return rand() % (max + 1 - min) + min;

}
void initValues() {

    S = new Matrix(tam, tam, 0);
    I = new Matrix(tam, tam, 0);
    R = new Matrix(tam, tam, 0.0);

    for (int i = 0; i < tam; i++)
        for (int j = 0; j < tam/2; j++)
        {
            S ->operator()(i, j, random(0.3, 1) * 1000);
            I ->operator()(i, j, random(0, 1) * 10);

        }
}

double f(double S, double I) {

    return  -B * S * I * I;
}
double g(double S, double I) {

    return B * S * I * I - (mi)*I;

}

double difussion(Matrix* ua, int x, int y, float coef) {

    //Difusão de Euler usada na aproximação de Ut+1 usada no calculo de Ru
    Matrix* u = ua;
    double difx, dify;
    double Y =  coef ;

    difx = (-2 * u->operator()(x, y) + u->operator()(x - 1, y) + u->operator()(x + 1, y)) * Y;

    dify = (-2 * u->operator()(x, y) + u->operator()(x, y - 1) + u->operator()(x, y + 1)) * Y;

    return difx + dify;

}
void step() {


    int i, j;
    for (i = 1; i < tam - 1; i++) {
        for (j = 1; j < tam - 1; j++) {

            double deltaS = (f(S->operator()(i, j), I->operator()(i, j)) + difussion(S, i, j, d11) + difussion(I, i, j, d12));
            double deltaI = (g(S->operator()(i, j), I->operator()(i, j)) + difussion(I, i, j, d22));
            double deltaR = (mi * I->operator()(i, j) + difussion(R, i, j, d33));
            S->operator()(i, j,dt* deltaS + S->operator()(i, j));
            R->operator()(i, j,dt* deltaR + R->operator()(i, j));
            I->operator()(i, j,dt* deltaI + I->operator()(i, j));
        }
    }


    for (int i = 0; i < tam; i++)
    {
        S->operator()(i, 0, S->operator()(i, 1));
        R->operator()(i, 0, R->operator()(i, 1));
        I->operator()(i, 0, I->operator()(i, 1));

        S->operator()(i, tam - 1, S->operator()(i, tam - 2));
        R->operator()(i, tam - 1, R->operator()(i, tam - 2));
        I->operator()(i, tam - 1, I->operator()(i, tam - 2));
    }

    for (int j = 0; j < tam; j++)
    {
        S->operator()(0, j, S->operator()(1, j));
        R->operator()(0, j, R->operator()(1, j));
        I->operator()(0, j, I->operator()(1, j));

        S->operator()(tam - 1, j, S->operator()(tam - 2, j));
        R->operator()(tam - 1, j, R->operator()(tam - 2, j));
        I->operator()(tam - 1, j, I->operator()(tam - 2, j));
    }


}
void printMatrixtoFile(Matrix* I, Matrix* S, Matrix* R, FILE* Iarq, FILE* Rarq, FILE* Sarq) {

    for (int v = 0; v < tam; v++)
    {
        int dx = 1;
        for (int j = 0; j < tam; j++)
        {

            fprintf(Iarq, "%f ", I->operator()(v, j)); //att
            fprintf(Rarq, "%f ", R->operator()(v, j)); //att
            fprintf(Sarq, "%f ", S->operator()(v, j)); //att


        }
        fprintf(Iarq, "\n");
        fprintf(Sarq, "\n");
        fprintf(Rarq, "\n");

    };




}

int main()
{
    initValues();
    FILE* Iarq, * Sarq, * Rarq;
    char* filename = (char*)"I.txt";
    Iarq = fopen(filename, "w");
    char* filename2 = (char*)"S.txt";
    Sarq = fopen(filename2, "w");
    char* filename3 = (char*)"R.txt";
    Rarq = fopen(filename3, "w");

    for (int i = 0; i < n; i++) {

        step();
        if (i % 10 == 0)
        {

            printMatrixtoFile(I, S, R, Iarq, Rarq, Sarq);

            if (i != n - 1)
            {
                fprintf(Iarq, "\n\n");
                fprintf(Sarq, "\n\n");
                fprintf(Rarq, "\n\n");
            }
        }

    }
    fclose(Iarq);
    fclose(Sarq);
    fclose(Rarq);
    saveGif((char*)"I.txt", (char*)"I.gif", 0.1, 0.01);
    saveGif((char*)"S.txt", (char*)"S.gif", 0.1, 0.01);

}

