#include <iostream>
#include <random>
#include <cmath>
#include <fstream>

using namespace std;

const int N = 20;
const double L = 10.0;
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

// Inicializa posiciones aleatorias sin superposición simple y velocidades aleatorias de módulo 1
void inicializarparticulas(double x[N], double y[N],
                           double vx[N], double vy[N],
                           double ax[N], double ay[N], double L) {
    mt19937 gen(random_device{}());
    double min_dist = 1.0;

    for (int i = 0; i < N; ++i) {
        bool valid;
        double xi, yi;
        do {
            valid = true;
            xi = rand_double(0, L, gen);
            yi = rand_double(0, L, gen);
            for (int j = 0; j < i; ++j) {
                double dx = xi - x[j];
                double dy = yi - y[j];
                dx -= L * round(dx / L);
                dy -= L * round(dy / L);
                if (sqrt(dx*dx + dy*dy) < min_dist) {
                    valid = false;
                    break;
                }
            }
        } while (!valid);

        double theta = rand_double(0, 2 * M_PI, gen);
        x[i] = xi;
        y[i] = yi;
        vx[i] = vi*cos(theta);
        vy[i] = vi*sin(theta);
        ax[i] = 0.0;
        ay[i] = 0.0;
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
   
    inicializarparticulas(x, y, vx, vy, ax, ay, L);

    fuerza(x, y, ax, ay, L);
    ofstream vxfile("vx.txt");
    ofstream vyfile("vy.txt");
    ofstream vfile("v.txt"); // módulo de la velocidad del sistema
    ofstream fout("posiciones3.txt");
    ofstream fout2("energias.txt");
  
    for (int i = 0; i < N; ++i)
        fout << x[i] << ", " << y[i] << endl;
    fout << endl;
    double sumaT = 0.0;
    int conteoT = 0;

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
        // Guardar velocidades de un rango de tiempo
        // Guardar velocidades al inicio y en equilibrio
        if ((k >= 70.0 && k <= 100.0)) {
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
        // Calcular y acumular temperatura instantánea SOLO entre k=20 y k=50
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
