#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

using namespace std;


unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;
#define NormRANu (2.3283063671E-10F)
#define Pi 3.14159265

#define densi 0.6
#define Tfin 1000
#define rad 1

const double a = Pi*rad*rad, d=0.2*rad;
const int NPart = 1800, NBoxes= sqrt(NPart/2);
const double l= sqrt(2*a/densi), L = NBoxes*l;

//Generador de números de Parisi-Rapuano.
void ini_ran(int SEMILLA)
{
    int INI,FACTOR,SUM,i;
    srand(SEMILLA);
    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;
    for(i=0;i<256;i++)
    {
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;
    }
    ind_ran=ig1=ig2=ig3=0;
}

float Random(void)
{
    float r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;
    return r;
}

double Distancia(double x, double y, double x1, double y1){
if(fabs(x-x1)<=L/2 && fabs(y-y1)<= L/2)
    return sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1));
if(fabs(x-x1)<=L/2){
    if(y>=L/2){
        return sqrt((x-x1)*(x-x1) + (y-y1-L)*(y-y1-L));
    }
    else{
        return sqrt((x-x1)*(x-x1) + (y-y1+L)*(y-y1+L));
    }
}
if(fabs(y-y1)<=L/2){
    if(x>=L/2){
        return sqrt((y-y1)*(y-y1) + (x-x1-L)*(x-x1-L));
    }
    else{
        return sqrt((y-y1)*(y-y1) + (x-x1+L)*(x-x1+L));
    }
}
else{
    if(y>=L/2 && x>=L/2){
        return sqrt((x-x1-L)*(x-x1-L) + (y-y1-L)*(y-y1-L));
    }
     if(y>=L/2 && x<L/2){
        return sqrt((x-x1+L)*(x-x1+L) + (y-y1-L)*(y-y1-L));
    }
     if(y<L/2 && x>=L/2){
        return sqrt((x-x1-L)*(x-x1-L) + (y-y1+L)*(y-y1+L));
    }
     if(y<L/2 && x<L/2){
        return sqrt((x-x1+L)*(x-x1+L) + (y-y1+L)*(y-y1+L));
    }
}
}

void Contorno(double &xnew,double &ynew, int (&Contorno)[NPart][2], int n){
if(xnew >= L){
    xnew = xnew - L;
    Contorno[n][0]++;
}
if(xnew < 0){
    xnew = xnew + L;
    Contorno[n][0]=Contorno[n][0]-1;
}
if(ynew>=L){
    ynew = ynew -L;
    Contorno[n][1]++;
}
if(ynew<0){
    ynew = ynew + L;
    Contorno[n][1]=Contorno[n][1]-1;
}
}

int Aceptar(double (&r)[NPart][2], double xnew, double ynew, int n){
for(int i=0; i<n; i++){
    if(Distancia(r[i][0], r[i][1], xnew, ynew)<= 2*rad)
        return 1;
}
for(int i=n+1; i<NPart; i++){
    if(Distancia(r[i][0], r[i][1], xnew, ynew)<= 2*rad)
        return 1;
}
return 0;
}


void PasoMonte(double (&r)[NPart][2], int (&NVueltas)[NPart][2]){
for(int i=0; i<NPart; i++){
    int n= NPart*Random(), flag;
    double xnew = r[n][0] + d*(Random()-0.5), ynew = r[n][1] + d*(Random()-0.5);
    Contorno(xnew, ynew, NVueltas, n);
    flag=Aceptar(r, xnew, ynew, n);
    if(flag ==0){
        r[n][0] = xnew;
        r[n][1] = ynew;
    }
}
}


void Guardardatos(double (&r)[NPart][2], int (&NVueltas)[NPart][2],int i){
char filename[50];
    sprintf(filename,"Evolucion/t%d.txt", i);
    FILE *f=fopen(filename, "wt");
    if(f==NULL)
        printf("No se pudo abrir el fichero \n");
    else{
        for(int j=0; j<NPart; j++){
            fprintf(f, "%lf \t %lf \t\n", r[j][0], r[j][1]);
            }
        fclose(f);
    }
}


void Inicializer(double (&r)[NPart][2], int (&NVueltas)[NPart][2] ){
for(int i=0; i<NBoxes; i++){
    for(int j=0; j<NBoxes; j++){
        r[i*NBoxes+j][0]=j*l;
        r[i*NBoxes+j][1]=i*l;
        r[NBoxes*NBoxes+i*NBoxes+j][0]=j*l+l/2;
        r[NBoxes*NBoxes+i*NBoxes+j][1]=i*l+l/2;
    }
}
for(int i=0; i< NPart; i++){
    NVueltas[i][0]=0;
    NVueltas[i][1]=0;
}
}

void Metropolis(void){
double particulas[NPart][2];
int NVueltas[NPart][2];
Inicializer(particulas, NVueltas);
Guardardatos(particulas, NVueltas, 0);
//Termalizo el sistema
for(int i=0; i<0; i++)
    PasoMonte(particulas, NVueltas);
for(int i=0; i<Tfin; i++){
    PasoMonte(particulas, NVueltas);
    Guardardatos(particulas, NVueltas, i+1);
}
}

void MSD(void){
FILE *fprom=fopen("MSD/delta = 0.3.txt", "wt");
    if(fprom != NULL){
        for(int i=0; i<=Tfin; i++){
        FILE *f0 = fopen("Evolucion/t0.txt", "rt");
        if(f0!=NULL){
            char filename[20];
            sprintf(filename,"Evolucion/t%d.txt", i);
            FILE *f=fopen(filename, "rt");
            if(f==NULL){
                printf("No se pudo abrir el fichero: %s.\n", filename);
            }
            else{
                double x0, y0, x, y, prom=0, desv=0;
                for(int j=0; j< NPart; j++){
                    fscanf(f0, "%lf\t%lf\n", &x0, &y0);
                    fscanf(f, "%lf\t%lf\n", &x, &y);
                    prom=prom+(x-x0)*(x-x0)+(y-y0)*(y-y0);
                    desv = desv + (x-x0)*(x-x0)+(y-y0)*(y-y0)*(x-x0)*(x-x0)+(y-y0)*(y-y0);
                }
                fprintf(fprom,"%d\t%lf\n", i, prom/NPart);
                fclose(f);
            }
            fclose(f0);
        }
        else{
            printf("No se pudo abrir el fichero: Evolucion/t0.txt");
            }
        }
    fclose(fprom);
}
else{
    printf("No se pudo abrir el fichero.");
}
}

int main(){
ini_ran(time(NULL));
Metropolis();
MSD();
return 0;
}
