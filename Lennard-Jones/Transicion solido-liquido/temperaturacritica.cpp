#include <iostream>
#include <random>
#include <cmath>
#include <fstream>

using namespace std;

const int N = 16;
const double L = 4.0;
const double dt = 0.004;
const int Pasos = 200;
const double sigma = 1.0; // Parámetro de Lennard-Jones
const double epsilon = 1.0; // Profundidad del pozo de potencial
const double m = 1.0; // Masa de las partículas
const double kB = 1.0; // Constante de Boltzmann
const double vi = 1.0; // Velocidad inicial de las partículas
double rand_double(double min, double max, mt19937 &gen) {
    uniform_real_distribution<double> dist(min, max);
    return dist(gen);
}
void reescalar_velocidades(double vx[N], double vy[N], double T_deseada) {
    double suma_vsqr = 0.0;
    for (int i = 0; i < N; ++i)
        suma_vsqr += vx[i]*vx[i] + vy[i]*vy[i];

    double T_actual = suma_vsqr / (2.0 * N);
    double factor = sqrt(T_deseada / T_actual);

    for (int i = 0; i < N; ++i) {
        vx[i] *= factor;
        vy[i] *= factor;
    }
}
void aplicarCondicionContornoDistancia(double& dx, double L) {
    if (dx > L / 2.0) dx -= L;
    else if (dx < -L / 2.0) dx += L;
}

void cotornoperiodicas(double &x, double L) {
    if (x < 0) x += L;
    if (x >= L) x -= L;
}

// Inicializa posiciones en red cuadrada y velocidades en reposo
void inicializar_cuadrado(double x[N], double y[N], double vx[N], double vy[N], double ax[N], double ay[N], double L) {
    int n = sqrt(N); // Para N=16, n_side=4
    double dx = L / n;
    int idx = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            x[idx] = (i + 0.5) * dx;
            y[idx] = (j + 0.5) * dx;
            vx[idx] = 0.0;
            vy[idx] = 0.0;
            ax[idx] = 0.0;
            ay[idx] = 0.0;
            idx++;
        }
    }
}

// Inicializa posiciones en red hexagonal y velocidades aleatorias 
void inicializar_hexagonal(double x[N], double y[N], double vx[N], double vy[N], double ax[N], double ay[N]) {
    int n = 4; // número de átomos por fila
    double dx = L / n;
    int idx = 0;

    for (int j = 0; j < n; ++j) { // filas (vertical)
        for (int i = 0; i < n; ++i) { // columnas (horizontal)
            x[idx] = i * dx + 0.5 * dx * (j % 2);
            y[idx] = j * dx * sqrt(3) / 2.0;
            vx[idx] = 0.0;
            vy[idx] = 0.0;
            ax[idx] = 0.0;
            ay[idx] = 0.0;
            idx++;
        }
    }
}

// Calcula fuerzas de Lennard-Jones y actualiza aceleraciones
void fuerza(double x[N], double y[N],
            double ax[N], double ay[N], double L) {
    for (int i = 0; i < N; ++i) {
        ax[i] = 0.0;
        ay[i] = 0.0;
    }
    for (int i = 0; i < N; ++i) {
        for (int j = i+1; j < N; ++j) {
                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                aplicarCondicionContornoDistancia(dx, L);
                aplicarCondicionContornoDistancia(dy, L);
                double r2 = dx*dx + dy*dy;
                if (r2 > 0.0*sigma && r2 < 9.0*sigma*sigma) {
                    double inv_r2 = 1.0 / r2;
                    double inv_r6 = inv_r2 * inv_r2 * inv_r2;
                    double inv_r12 = inv_r6 * inv_r6;
                    double f = 24.0 * epsilon * inv_r2 * (2.0 * inv_r12 - inv_r6);

                    ax[i] += f * dx;
                    ay[i] += f * dy;
                    ax[j] -= f * dx;
                    ay[j] -= f * dy;
                
            }
        }
    }
}
//Energia cinetica del sistema
double energia_cinetica(double vx[N], double vy[N]) {
    double ec = 0.0;
    for (int i = 0; i < N; ++i) {
        ec += 0.5 * m * (vx[i] * vx[i] + vy[i] * vy[i]);
    }
    return ec;
}
// Energia potencial del sistema
double energia_potencial(double x[N], double y[N], double L) {
    double ep = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            aplicarCondicionContornoDistancia(dx, L);
            aplicarCondicionContornoDistancia(dy, L);
            double r2 = dx * dx + dy * dy;
            if (r2 > 0.0 && r2 < 9.0 * sigma * sigma) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                ep += 4 * epsilon * (pow(sigma,12)/r12 - pow(sigma,6)/r6);
            }
        }
    }
    return ep;

}


int main() {
    double x[N], y[N], vx[N], vy[N], ax[N], ay[N];
    double ec[N], ep[N];
    double sumaT = 0.0;
    int conteoT = 0;
    inicializar_cuadrado(x,y,vx,vy,ax,ay,L); // Inicializa posiciones en red hexagonal
   // Guardar posición inicial de una partícula para medir fluctuación
    double x0_i = x[0], y0_i = y[0]; // Posición inicial partícula 0
    double x0_j = x[1], y0_j = y[1]; // Posición inicial partícula 1

    ofstream sep_ij("Separacion_ij.txt");

    
   //Inicializar fuerzas
    fuerza(x, y, ax, ay, L);
    ofstream vxfile("vx_sl.txt");
    ofstream vyfile("vy_sl.txt");
    ofstream vfile("v_sl.txt"); // módulo de la velocidad del sistema
    ofstream fout("posiciones_sl.txt");
    ofstream fout2("energias_sl.txt");
    ofstream fluct ( "Fluctuaciones.txt");

    //Algoritmo de integración de Verlet 
    for (double k = 0.0; k < Pasos ; k += dt) {
        double wx[N], wy[N];

         // REESCALAR velocidades en t = 20, 30, 35, 45
  if ( (k >= 60.000 && k < 60.000 + dt) ||
     (k >= 120.000 && k < 120.000 + dt) ||
     (k >= 180.000 && k < 180.000 + dt) ) {
    for (int i = 0; i < N; ++i) {
        vx[i] *= 1.1;
        vy[i] *= 1.1;
    }
    cout << "Reescalado por 1.1 en t = " << k << endl;
}

        
        // Comienzo del algoritmo de Verlet
        for (int i = 0; i < N; ++i) {
            x[i] += vx[i] * dt + 0.5 * ax[i] * dt * dt;
            y[i] += vy[i] * dt + 0.5 * ay[i] * dt * dt;
            cotornoperiodicas(x[i], L);
            cotornoperiodicas(y[i], L);

            wx[i] = vx[i] + 0.5 * ax[i] * dt;
            wy[i] = vy[i] + 0.5 * ay[i] * dt;
        }

        for (int i = 0; i < N; ++i)
            fout << x[i] << ", " << y[i] << endl;
        fout << endl;

        fuerza(x, y, ax, ay, L);

        for (int i = 0; i < N; ++i) {
            vx[i] = wx[i] + 0.5 * ax[i] * dt;
            vy[i] = wy[i] + 0.5 * ay[i] * dt;
        }


        // Guardar velocidades de un rango de tiempo
        // Guardar velocidades al inicio y en equilibrio
        if ( (k >= 70.0 && k <= 100.0)) {
        for (int i = 0; i < N; ++i) {
        double vmod = sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
        vxfile << vx[i] << endl;
        vyfile << vy[i] << endl;
        vfile << vmod << endl;
        }
        }

     
   
        
       
        // Guardar energías
        double ec_total = energia_cinetica(vx, vy);
        double ep_total = energia_potencial(x, y, L);
        double etotal = ec_total + ep_total;
        fout2 << ec_total << ", " << ep_total << ", " << etotal << ", " << k << endl;
        // Temperatura instantánea usando el teorema de equipartición
         double suma_vsqr = 0.0;
        for (int i = 0; i < N; ++i) {
            suma_vsqr += vx[i] * vx[i] + vy[i] * vy[i];
        }
        double temp_inst = suma_vsqr / (2.0 * N);
        if (k >= 70.0 && k <= 100.0) {
            sumaT += temp_inst;
            conteoT++;
        }

    double dx_ij = x[0] - x[1];
    double dy_ij = y[0] - y[1];
    aplicarCondicionContornoDistancia(dx_ij, L);
    aplicarCondicionContornoDistancia(dy_ij, L);
    double dr2_ij = dx_ij*dx_ij + dy_ij*dy_ij;

sep_ij << k << ", " << dr2_ij << ", " << temp_inst << endl;
}

    // Temperatura promedio temporal
    double temp_prom = sumaT / conteoT;
    std::cout << "Temperatura promedio temporal: T = " << temp_prom << endl;
    // Calcular temperatura final usando el teorema de equipartición    
    double suma_vsqr = 0.0;
    for (int i = 0; i < N; ++i) {
    suma_vsqr += vx[i] * vx[i] + vy[i] * vy[i];
    }
    double temperatura = suma_vsqr / (2.0 * N);
    std::cout << "Temperatura final (teorema de equipartición): T = " << temperatura << endl;
        fout2.close();
    fluct.close();
    vxfile.close();
    vyfile.close();
    vfile.close();
    sep_ij.close();
    fout.close();
    std:: cout<< "Simulación completada. Posiciones guardadas en posiciones.txt\n";
    return 0;
}
