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

#define densi 0.9
#define Tfin 1000
#define rad 1
#define beta 1
#define g 0.01

const double a = Pi*rad*rad, d=0.5*rad;
const int NPart = 1000, NBoxesX= 50, NBoxesY=10;
const double lx= sqrt(2*a/(sqrt(3)*densi)), ly = sqrt(3)*lx, LX = NBoxesX*lx, LY=10*LX;

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
if(fabs(x-x1)<=LX/2)
    return sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1));
else{
    if(x>= LX/2)
        return sqrt((x-x1-LX)*(x-x1-LX) + (y-y1)*(y-y1));
    else
        return sqrt((x-x1+LX)*(x-x1+LX) + (y-y1)*(y-y1));
}
}

void Contorno(double &xnew,double &ynew, int (&Contorno)[NPart][2], int n){
if(xnew >= LX){
    xnew = xnew - LX;
    Contorno[n][0]++;
}
if(xnew < 0){
    xnew = xnew + LX;
    Contorno[n][0]=Contorno[n][0]-1;
}
}

int Aceptar(double (&r)[NPart][2], double xnew, double ynew, int n){
if(ynew <= rad || ynew>= LY-rad)
    return 1;
if(Random()>exp(beta*g*(r[n][1]-ynew)))
    return  1;
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
    double xnew = r[n][0] + d*(Random()-0.5), ynew = r[n][1] + 5*d*(Random()-0.5);
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
FILE *fini=fopen("Evolucion/t1000.txt", "rt");
if(fini!= NULL){
    for(int i=0; i<NPart; i++){
        fscanf(fini, "%lf \t%lf\n", &(r[i][0]), &(r[i][1]));
    }
  fclose(fini);
}
else{
    printf("No se ha podido encontrar el fichero con las condiciones iniciales.\nSe crearan unas condiciones iniciales cristalinas.\n");
    for(int i=0; i<NBoxesY; i++){
    for(int j=0; j<NBoxesX; j++){
        r[i*NBoxesX+j][0]=j*lx;
        r[i*NBoxesX+j][1]=i*ly+rad;
        r[NBoxesX*NBoxesY+i*NBoxesX+j][0]=j*lx+lx/2;
        r[NBoxesX*NBoxesY+i*NBoxesX+j][1]=i*ly+ly/2+rad;
    }
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
for(int i=0; i<5000; i++)
    PasoMonte(particulas, NVueltas);
printf("He acabado la termalizacion del sistema.\n");
for(int i=0; i<Tfin; i++){
    for(int n = 0; n<200; n++)
        PasoMonte(particulas, NVueltas);
    Guardardatos(particulas, NVueltas, i+1);
}
}

void AverageY(void){
char filename[20];
FILE *fop=fopen("AverageY.txt", "wt");
if(fop!=NULL){
    double prom =0, x, y;
    for(int i=1; i<=Tfin; i++){
        prom=0;
        sprintf(filename, "Evolucion/t%d.txt", i);
        FILE *f = fopen(filename, "rt");
        if(f!=NULL){
            for(int j=0; j<NPart; j++){
                fscanf(f,"%lf\t%lf\n", &x, &y);
                prom = prom+y;
            }
            fprintf(fop, "%d\t%lf\n", i, prom/NPart);
            fclose(f);
        }
        else{
            printf("No se pudo abrir el fichero\n");
        }
    }
    fclose(fop);
}

}

void Perfil(void){
FILE *fperf = fopen("Perfil/T=1,g=1.txt", "wt");
double Perfil[100], PerfilProm[100], s;
if(fperf != NULL) {
    char filename[20];
    double x, y;
    for(int i=1; i<=int(Tfin); i++){
        sprintf(filename, "Evolucion/t%d.txt", i);
        FILE *fact=fopen(filename, "rt");
        if(fact != NULL){
            for(int j=0; j<100; j++)
                Perfil[j]=0;
            for(int n=0; n<NPart; n++){
                fscanf(fact, "%lf\t%lf\n", &x, &y);
                for(int j=0; y>j*10 || j<100; j++){
                    if((j+1)*10>=y+rad && (j)*10<=y-rad){
                        Perfil[j]=Perfil[j]+a/(LX*10);
                        break;
                    }
                    if(y-rad<j*10){
                        s= a*2*acos(y-j*10.0)/(2*Pi)-(y-j*10.0)*sqrt(((y-j*10.0))*((y-j*10.0)) +1.0);
                        if(fabs(y-j*10) > 1){
                            printf("%e\t %lf\t %lf\n", y,j*10.0, s);
                        }
                        Perfil[j] = Perfil[j]+ (a-s)*1.0/(LX*10.0);
                        Perfil[j-1]=Perfil[j-1] + s*1.0/(LX*10.0);
                        break;
                    }
                    if(y+rad>(j+1)*10 && y<(j+1)*10 ){
                        s= a*2*acos((j+1)*10.0-y)/(2*Pi)-((j+1)*10.0-y)*sqrt(((j+1)*10.0-y)*((j+1)*10.0-y) +1.0);
                        Perfil[j] = Perfil[j]+ (a-s)/(LX*10);
                        Perfil[j+1]=Perfil[j+1] + s/(LX*10);
                        break;
                    }
                }
            }
            for(int j=0; j<100; j++){
                PerfilProm[j] = PerfilProm[j]+Perfil[j];
            }
            fclose(fact);
        }
        else{
            printf("No se pudo abrir el fichero: %s\n", filename);
        }
    }
    for(int j=0; j<100; j++){
        fprintf(fperf, "%lf\t %lf\n", j*10+5.0, PerfilProm[j]/int(Tfin));
    }
fclose(fperf);
}
else{
    printf("No se pudo abrir el fichero.\n");
}
}

int main(){
ini_ran(time(NULL));
printf("%lf \n", NPart*a/(LX*LY));
Metropolis();
printf("He acabado la simulacion.\nEstoy calculando el perfil.\n");
AverageY();
Perfil();
return 0;
}
