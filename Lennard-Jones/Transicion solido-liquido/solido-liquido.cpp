#include <iostream>
#include <random>
#include <cmath>
#include <fstream>

using namespace std;

const int N = 16;
const double L = 4.0;
const double dt = 0.002;
const int Pasos = 100;
const double sigma = 1.0; // Parámetro de Lennard-Jones
const double epsilon = 1.0; // Profundidad del pozo de potencial
const double m = 1.0; // Masa de las partículas
const double kB = 1.0; // Constante de Boltzmann
const double vi = 1.0; // Velocidad inicial de las partículas
double rand_double(double min, double max, mt19937 &gen) {
    uniform_real_distribution<double> dist(min, max);
    return dist(gen);
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
void inicializar_solid(double x[N], double y[N], double vx[N], double vy[N], double ax[N], double ay[N], double L) {
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
    inicializar_solid(x, y, vx, vy, ax, ay, L);

    fuerza(x, y, ax, ay, L);
    ofstream vxfile("vx_sl.txt");
    ofstream vyfile("vy_sl.txt");
    ofstream vfile("v_sl.txt"); // módulo de la velocidad del sistema
    ofstream fout("posiciones_sl.txt");
    ofstream fout2("energias_sl.txt");
    ofstream fluct ( "Fluctuaciones.txt");
  // Guarda la posición inicial de la partícula 0
    double x0_inicial = x[0];
    double y0_inicial = y[0];
    for (int i = 0; i < N; ++i)
        fout << x[i] << ", " << y[i] << endl;
    fout << endl;
    //Algoritmo de integración de Verlet 
    for (double k = 0.0; k < Pasos ; k += dt) {
        double wx[N], wy[N];
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

        // Reescalado de velocidades en tiempos específicos
    if (
    abs(k - 20.0) < 0.5*dt ||
    abs(k - 30.0) < 0.5*dt ||
    abs(k - 35.0) < 0.5*dt ||
    abs(k - 45.0) < 0.5*dt
    ) {
    for (int i = 0; i < N; ++i) {
        vx[i] *= 1.5;
        vy[i] *= 1.5;
    }
    std::cout << "Velocidades reescaladas en t = " << k << std::endl;
}

// Calcular fluctuación de posición de la partícula 0
double dx_fluct = x[0] - x0_inicial;
double dy_fluct = y[0] - y0_inicial;

// Aplicar condiciones de contorno mínimas a la diferencia:
aplicarCondicionContornoDistancia(dx_fluct, L);
aplicarCondicionContornoDistancia(dy_fluct, L);

double r2_fluct = dx_fluct*dx_fluct + dy_fluct*dy_fluct;

// Guardar el tiempo y la fluctuación
fluct << k << " " << r2_fluct << endl;

        // Guardar velocidades de un rango de tiempo
        // Guardar velocidades al inicio y en equilibrio
        if ( (k >= 70.0 && k <= 1000.0)) {
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
        // Calcular temperatura instantánea
         // Calcular temperatura instantánea
         // Solo se guarda la temperatura en el rango de 70 a 100
       
         double suma_vsqr = 0.0;
        for (int i = 0; i < N; ++i) {
            suma_vsqr += vx[i] * vx[i] + vy[i] * vy[i];
        }
        double temp_inst = suma_vsqr / (2.0 * N);
        if (k >= 70.0 && k <= 100.0) {
            sumaT += temp_inst;
            conteoT++;
        }
    }
    fout2.close();
 

    vxfile.close();
    vyfile.close();
    vfile.close();

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
     

    fout.close();
    std:: cout<< "Simulación completada. Posiciones guardadas en posiciones.txt\n";
    return 0;
}
