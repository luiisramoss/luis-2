#include <stdio.h>
#include <cmath>
#include <iostream>

void aceleracion(double m[9], double x[9], double y[9], double ax[9], double ay[9]);
void QUITARreescalar(double x[9], double y[9]);
void shorten(double vy[9]);
void wakeup(double vy[9]);
void quitar_reescalarR(double* V);
void posicion(double x[9], double y[9], double vx[9], double vy[9], double ax[9], double ay[9], double h, int i);
void velocidadauxiliar(double vx[9],double wx[9], double wy[9], double vy[9], double ax[9], double ay[9], double h, int i);
void velocidad(double vx[9], double vy[9], double wx[9], double wy[9], double ax[9], double ay[9], double h, int i);


//voy a crear un fichero (input) con las condiciones iniciales de los planetas y del sol
//a partir de los datos de ese fichero y del algoritmo de Verlet voy a calcular las siguientes variables
//las voy a ir guardando paso a paso en el nuevo fichero (output), que va a ser el que utilice en el .py
//mi fichero imput esta estructurado así:
//
//    xsol ysol vsolx vsoly masasol
//    xmerc ymerc vmercx vmercy masamerc
//    xvenus yvenus vvenusx vvenusy masaven
//    ...
//
//
//    xnept ynept vneptx vnepty masanept
//
//  todos los datos de mi fichero imput han sido modificados para evitar un error tan grande


int main() {
    // Creo y abro los ficheros que voy a usar
    FILE *input;
    FILE *output;
    FILE *output2;
    FILE *output3;

    input = fopen("datosREESCALADOS.txt", "r");
    output = fopen("valoresfinal.txt", "w");
    output2= fopen ("Energias.txt", "w");
    output3= fopen ("Periodos.txt", "w");

    double x[9], y[9], vx[9], vy[9], ax[9], ay[9], wx[9], wy[9], m[9];
    double E[9], Ec[9], Ep[9];
    double xauxant[9], xaux[9], periodo[9], cont[9];


    double t = 0, h = 0.01; // Inicializo t y defino el paso h
    int i, j;

    // Inicializo mis vectores con los valores iniciales correspondientes al sol y a cada planeta
    for(i=0; i<9; i++){
        fscanf(input, "%lf\t%lf\t%lf\t%lf\t%lf",&x[i],&y[i],&vx[i],&vy[i],&m[i]);

        //hago que xauxant sea el primer x[i]
            xauxant[i]=x[i];

        //inicializo mi periodo y mi contador
        periodo[i]=0.;
        cont[i]=0.;

    //Incluyo la situacion inicial en los resultados para que empiecen desde el mismo lugar
        fprintf(output, "%lf,%lf\n", x[i], y[i]);
    }
     fprintf(output, "\n");

    //Calculo la aceleracion inicial
    aceleracion(m,x,y,ax,ay);

    while (t <100) {
        
       
        // Inicio el bucle para guardar los resultados en los vectores
        for (i = 0;i<9; i++) {
            
           
            // Posición en x e y cuando t+h
            posicion(x, y, vx, vy, ax, ay, h, i);

            //vector auxiliar
            velocidadauxiliar(vx, wx, wy, vy, ax, ay, h, i);

            //actualizo la aceleracion
            aceleracion(m, x, y, ax, ay);

            // Velocidades en t+h
            velocidad(vx, vy, wx, wy, ax, ay, h, i);
            if((xaux[i]>xauxant[i])&&(xaux[i]>x[i])){
                periodo[i]=t; //me quedo con el momento en el que se da mi ultimo periodo
                cont[i]++; //contador para despues hacer la media con la sumatoria
            }

            //actualizo mis xant y xahora para que el bucle del periodo siga corriendo
            xauxant[i]=x[i]-wx[i]*h;
            xaux[i]=x[i];

            x[0]=0.;
            y[0]=0.; //el sol no se mueve

            
            fprintf(output, "%lf, %lf\n", x[i], y[i]); //imprimo las orbitas
        }


        for (i = 0;i<9; i++) {
            // Energía de cada planeta
    
            Ec[i]=0.5*m[i]*(pow(vx[i],2.)+pow(vy[i],2.));
            Ep[i]=0.;
            for(j=0;j<9;j++){
                if(i!=j){
                     Ep[i]= Ep[i]-(1.*m[j]*m[i]/(1.*pow((pow((x[i]-x[j]),2.)+pow((y[i]-y[j]),2.)),0.5)));
                }
            }
            E[i]=(Ec[i]+Ep[i])*10e5;
            fprintf(output2, "%e, %e\n", t, E[i]);
        }
        fprintf(output, "\n");
        fprintf(output2, "\n");

       
        t = t + h;
    }

    //Imprimimos los periodos
        for (i = 0;i<9; i++){
            fprintf(output3, "%lf\n", periodo[i]*58.1/cont[i]);

        }

    fclose(input);
    fclose(output);

    return 0;
}

void aceleracion(double m[9], double x[9], double y[9], double ax[9], double ay[9]){
    int i, j;
    
    //Calculo la aceleracion de cada elemento i causada por el resto de elementos j
    for(i=0; i<9; i++){
        ax[i]=0.;
        ay[i]=0.;
        for(j=0; j<9; j++){
            if(j!=i){
            ax[i]=ax[i]+((-m[j]*(x[i]-x[j]))/(1.*pow((pow((x[i]-x[j]),2.)+pow((y[i]-y[j]),2.)),1.5)));
            ay[i]=ay[i]+((-m[j]*(y[i]-y[j]))/(1.*pow((pow((x[i]-x[j]),2.)+pow((y[i]-y[j]),2.)),1.5)));
            }
            
        }
    }
    return;
}



void quitar_reescalarR(double* V)
{
    double c=1496e+8; // distancia entre la tierra y el sol en metros
    for (int i = 0; i < 9; i++)
    {
        V[i] = V[i] * c;
    }
    return;
}

void shorten(double vy[9])
{
    double c=1496e+8; 
    for (int i = 0; i < 9; i++)
    {
        vy[i] = vy[i] / c;
    }
}
void wakeup(double vy[9])
{
    double c=1496e+8; // distancia entre la tierra y el sol en metros
    double G =6.674e-11;
    for (int i = 0; i < 9; i++)
    {
        vy[i] = vy[i] / sqrt((G * (1988475e+24)) / pow(c, 3));
    }
}
void reescalarVelocidady(double vy[9])
{
    double c=1496e+8; // distancia entre la tierra y el sol en metros
    double G =6.674e-11;
    for (int i = 0; i < 9; i++)
    {
        vy[i] = (vy[i]/c)/sqrt((G * (1988475e+24)) / pow(c, 3)); // Aplicar shorten
    }
}
void QUITARreescalar(double x[9], double y[9]){
    double c=1496e+8; // distancia entre la tierra y el sol en metros
    for (int i=0; i<9;i++){
        x[i]=x[i]*1496e+8;
        y[i]=y[i]*1496e+8;
    }
}
void posicion(double x[9], double y[9], double vx[9], double vy[9], double ax[9], double ay[9], double h, int i){
            // Posición en x e y cuando t+h
            for (i = 0;i<9; i++) {
            x[i] = x[i] + h * vx[i] + h * h / 2. * ax[i]; //defino mi x posterior
            y[i] = y[i] + h * vy[i] + h * h / 2. * ay[i];
           
        }
        return;
}
void velocidadauxiliar(double vx[9],double wx[9], double wy[9], double vy[9], double ax[9], double ay[9], double h, int i){
            //vector auxiliar
            for (i = 0;i<9; i++) {
                wx[i] = vx[i] + h / 2. * ax[i];
                wy[i] = vy[i] + h / 2. * ay[i];

        }
        return;
}
void velocidad(double vx[9], double vy[9], double wx[9], double wy[9], double ax[9], double ay[9], double h, int i){
            // Velocidades en t+h
            for (i = 0;i<9; i++) {
                vx[i] = wx[i] + h/2. * ax[i];
                vy[i] = wy[i] + h/2. * ay[i];
        }
        return;
}
