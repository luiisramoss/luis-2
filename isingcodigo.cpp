#include <iostream>
#include <vector>
#include <cstdlib> // Para rand() y srand()
#include <ctime>   // Para time()
#include <cmath>   // Para exp()
#include <fstream> // Para guardar la matriz en un archivo
using namespace std;
#define N 100 // Definir el tamaño de la matriz

double calcular_energia(int ma[N][N]) {
  double energia = 0;

  // Iterar sobre cada elemento de la matriz
  for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
          // Condiciones de contorno periódicas
          int spin = ma[i][j];
          int derecha = ma[i][(j + 1) % N];       // Vecino a la derecha
          int izquierda = ma[i][(j - 1 + N) % N]; // Vecino a la izquierda
          int abajo = ma[(i + 1) % N][j];         // Vecino abajo
          int arriba = ma[(i - 1 + N) % N][j];    // Vecino arriba

          // Sumar las interacciones con los vecinos
          energia += spin * (derecha + izquierda + abajo + arriba);
      }
  }

  // Multiplicar por el factor -1/2
  energia *= -0.5;

  return energia;
}

void guardarEnergiaEnFichero(double energia, const string &nombreFichero) {
  ofstream fichero(nombreFichero, ios::app); // Abrir el fichero en modo de añadir (append)
  if (fichero.is_open()) {
      fichero << energia << endl; // Escribir la energía en una nueva línea
      fichero.close();
  } else {
      cerr << "Error al abrir el fichero " << nombreFichero << endl;
  }
}


void guardarMatrizEnFichero(const int ma[N][N], const string &nombreFichero) {
  ofstream fichero(nombreFichero, ios::app); // Abrir el fichero en modo de añadir (append)
  if (fichero.is_open()) {
      for (int i = 0; i < N; i++) {
          for (int j = 0; j < N; j++) {
              fichero << ma[i][j]; // Escribir el valor de la matriz
              if (j < N - 1) {
                  fichero << ", "; // Agregar una coma entre los valores, excepto al final de la fila
              }
          }
          fichero << endl; // Nueva línea después de cada fila
      }
      fichero << endl; // Línea en blanco para separar matrices
      fichero.close();
  } else {
      cerr << "Error al abrir el fichero " << nombreFichero << endl;
  }
}
double evaluarpuntos(int &n, int &m, double &p, int ma[N][N], double T) {
  double ex, x, de; // Definir la matriz de 5x5

  if (n==N-1 && m<N-1 && m>0) { // Si es un borde
    de=2*ma[n][m]*(ma[0][m]+ma[n-1][m]+ma[n][m+1]+ma[n][m-1]); // Evaluar el punto en la matriz
   
  } else if  (m==N-1 && n<N-1 && n>0) { // Si no es un borde
    // Calcular la energía del punto en la matriz
    de=2*ma[n][m]*(ma[n+1][m]+ma[n-1][m]+ma[n][0]+ma[n][m-1]); // Evaluar el punto en la matriz
 
 
  } else if (n==0 && n<N-1  && m<N-1 && m>0) { // Si no es un borde
  de=2*ma[n][m]*(ma[n+1][m]+ma[N-1][m]+ma[n][m+1]+ma[n][m-1]); // Evaluar el punto en la matriz
  
  } else if (m==0 && n<N-1 && n>0) { // Si no es un borde
    de=2*ma[n][m]*(ma[n+1][m]+ma[n-1][m]+ma[n][m+1]+ma[n][N-1]); // Evaluar el punto en la matriz
  
  } else if (n==0 && m==0) { // Si no es un borde
    de=2*ma[n][m]*(ma[n+1][m]+ma[N-1][m]+ma[n][m+1]+ma[n][N-1]); // Evaluar el punto en la matriz
  } else if (n==N-1 && m==N-1) { // Si no es un borde
    de=2*ma[n][m]*(ma[0][m]+ma[N-1][m]+ma[n][0]+ma[n][N-1]); // Evaluar el punto en la matriz
  } else if (n==0 && m==N-1) { // Si no es un borde
    // Calcular la energía del punto en la matriz
    de=2*ma[n][m]*(ma[0][m]+ma[N-1][m]+ma[n][0]+ma[n][N-1]); // Evaluar el punto en la matriz
  } else if (n==N-1 && m==0) { // Si no es un borde
    // Calcular la energía del punto en la matriz
    de=2*ma[n][m]*(ma[0][m]+ma[N-1][m]+ma[n][0]+ma[n][N-1]); // Evaluar el punto en la matriz
  } else de=2*ma[n][m]*(ma[n+1][m]+ma[n-1][m]+ma[n][m+1]+ma[n][m-1]); // Evaluar el punto en la matriz
  

  x = (-de)/T; // Calcular la energía del punto en la matriz
  ex= exp(x); // Calcular la exponencial
  if (1<ex) { // Si la exponencial es mayor que 1
   p=1.0; // Asignar 1 a p
  } else { // Si no
    p=ex; // Asignar la exponencial a p
    }
    return p;
  }
// Generar un número aleatorio decimal entre 0 y 1
double numeroDecimalAleatorio() {
  return static_cast<double>(rand()) / RAND_MAX;
}

int matrizespines( int ma[N][N]) {
 


    // Llenar la matriz con valores aleatorios
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            ma[i][j] = (rand() % 2 == 0) ? 1 : -1; // Generar -1 o 1 aleatoriamente
        }
    }

    // Imprimir la matriz
    cout << "Matriz " << N << "x" << N << " con valores aleatorios entre -1 y 1:" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << ma[i][j] << " ";
        }
        cout << endl;
    }
    return ma[N][N];
}
int cambiosigno(int &n, int &m, double p, double epsilon, int ma[N][N]) {
    // Cambiar el signo de los valores en la matriz
    if (p > epsilon) { // Si p es mayor que 0.5
       ma[n][m] = -ma[n][m]; // Cambiar el signo del valor en la matriz
    }
  return ma[N][N];
}
int puntoAleatorioEntero(int max, int &n) {
    // Generar valores aleatorios enteros entre 0 y max
   return n = rand() % (max);
   
}
/*void energía(int ma[N][N], int i, int j, int E){
  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      E=ma[i][j]*(ma[][])
    }

  }
}*/


int main() { 
    int t=0;
    int pasos_MC=N*N;
    int n, m,k; // Variables para almacenar los puntos aleatorios
    double p, T=2.0; // Inicializar la temperatura y la probabilidad
    int ma[N][N]; // Inicializar la matriz de 5x5
    srand(time(0));
    matrizespines(ma); // Llamar a la función para crear y mostrar la matriz
    guardarMatrizEnFichero(ma, "matriz.txt"); // Guardar la matriz en un archivo

  /*for (int t = 0; t < pasos_MC; t++) { // pasos_MC = número de pasos de Monte Carlo que quieras
    for (int k = 0; k < N*N; k++) {
        // Escoge un espín aleatorio
        puntoAleatorioEntero(N, n);
        puntoAleatorioEntero(N, m);

        // Calcula ΔE y la probabilidad de aceptación
        evaluarpuntos(n, m, p, ma, T);
        double epsilon = numeroDecimalAleatorio();

        // Aplica la regla de Metropolis
        cambiosigno(n, m, p, epsilon, ma);
    }

    // Guarda energía y configuración después de cada paso MC
    double energia = calcular_energia(ma);
    guardarEnergiaEnFichero(energia, "energia1.1.txt");
    guardarMatrizEnFichero(ma, "matriz1.1.txt");
}
 */

    for (int i = 0; i < N*N; i++) {
    for (int t = 0; t <= N*N; t++) {
    puntoAleatorioEntero(N, n);
    puntoAleatorioEntero(N,m); // Llamar a la función para obtener un punto aleatorio entero
    evaluarpuntos(n, m, p, ma, T); // Evaluar el punto en la matriz
    //cout << "Punto aleatorio  seleccionado: (" << n << ", " << m << ")" << endl;
    double epsilon = numeroDecimalAleatorio();
    //cout << "Número decimal aleatorio entre 0 y 1: " << epsilon << endl;
    cambiosigno(n, m, p, epsilon, ma); // Llamar a la función para cambiar el signo 
    // Cambiar el signo, si corresponde, de la posición de la matriz incial n y m
    //cout << "Matriz después de cambiar el signo:" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << ma[i][j] << " "; // Imprimir la matriz después de cambiar el signo
        }
        cout << endl;

    }
    cout << endl; // Imprimir una línea en blanco
   
    if (t==N*N){
    calcular_energia(ma); // Calcular la energía de la matriz
    guardarEnergiaEnFichero(calcular_energia(ma), "energia1.txt"); // Guardar la energía en un archivo
    guardarMatrizEnFichero(ma, "matriz1.2.txt"); // Guardar la matriz en un archivo
    }
    //cout << "Probabilidad: " << p << endl; // Imprimir la probabilidad
    //cout << "Epsilon: " << epsilon << endl; // Imprimir la temperatura
  }
 
}
return 0;
}

