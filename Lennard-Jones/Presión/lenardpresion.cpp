#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

const int N = 20;
const double L = 10.0;
const double h = 0.002;
const int Pasos = 100;
const double sigma = 1.0;
const double epsilon = 1.0;
const double m = 1.0;
const double kB = 1.0;
double x_ant[N], y_ant[N];
double cte = 0.1;

double rand_double(double min, double max, mt19937 &gen) {
    uniform_real_distribution<double> dist(min, max);
    return dist(gen);
}

double temperatura(double vx[N], double vy[N]) {
    double suma_vsqr = 0.0;
    for (int i = 0; i < N; ++i) {
        suma_vsqr += vx[i] * vx[i] + vy[i] * vy[i];
    }
    return suma_vsqr / (2.0 * N);
}

double presion(double x[N], double y[N], double vx[N], double vy[N], double L) {
    double delta_p = 0.0;
    for (int i = 0; i < N; ++i) {
        if (x_ant[i] > L - cte && x[i] < cte) delta_p +=  m * fabs(vx[i]);
        else if (x_ant[i] < cte && x[i] > L - cte) delta_p +=  m * fabs(vx[i]);
        if (y_ant[i] > L - cte && y[i] < cte) delta_p +=  m * fabs(vy[i]);
        else if (y_ant[i] < cte && y[i] > L - cte) delta_p += m * fabs(vy[i]);
    }
    return delta_p / (L * L * h);
}

void aplicarCondicionContornoDistancia(double& dx, double L) {
    if (dx > L / 2.0) dx -= L;
    else if (dx < -L / 2.0) dx += L;
}

void cotornoperiodicas(double &x, double L) {
    if (x < 0) x += L;
    if (x >= L) x -= L;
}

void inicializar(double x[N], double y[N], double vx[N], double vy[N], double ax[N], double ay[N],
                 double velocidad, bool casoEspecial = false) {
    mt19937 gen(random_device{}());

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
                if (sqrt(dx * dx + dy * dy) < 1.0) {
                    valid = false;
                    break;
                }
            }
        } while (!valid);

        x[i] = xi;
        y[i] = yi;

        if (!casoEspecial) {
            double theta = rand_double(0, 2 * M_PI, gen);
            vx[i] = velocidad * cos(theta);
            vy[i] = velocidad * sin(theta);
        } else {
            vx[i] = rand_double(0, 1, gen);
            vy[i] = 0.0;
        }

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
void simular(double velocidad_inicial, bool casoEspecial, ofstream &fout_total) {
    double x[N], y[N], vx[N], vy[N], ax[N], ay[N];
    inicializar(x, y, vx, vy, ax, ay, velocidad_inicial, casoEspecial);

    string nombre_archivo = "posiciones_vi";
    nombre_archivo += (casoEspecial ? "especial" : to_string(int(velocidad_inicial)));
    nombre_archivo += ".txt";

    ofstream fout_pos(nombre_archivo);

    double sumaP = 0.0, sumaT = 0.0;
    int conteo = 0;

     for (double k = 0.0; k < Pasos ; k += h) {
    

        for (int i = 0; i < N; ++i) {
            x_ant[i] = x[i];
            y_ant[i] = y[i];
            x[i] += vx[i] * h + 0.5 * ax[i] * h * h;
            y[i] += vy[i] * h + 0.5 * ay[i] * h * h;
            cotornoperiodicas(x[i], L);
            cotornoperiodicas(y[i], L);
        }

        double wx[N], wy[N];
        for (int i = 0; i < N; ++i) {
            wx[i] = vx[i] + 0.5 * ax[i] * h;
            wy[i] = vy[i] + 0.5 * ay[i] * h;
        }

        fuerza(x, y, ax, ay, L);

        for (int i = 0; i < N; ++i) {
            vx[i] = wx[i] + 0.5 * ax[i] * h;
            vy[i] = wy[i] + 0.5 * ay[i] * h;
        }

        if (k > 50.0) {
            double T = temperatura(vx, vy);
            double P = presion(x, y, vx, vy, L);
            sumaT += T;
            sumaP += P;
            conteo++;
        }

        // Guardar posiciones en el nuevo formato solicitado
        for (int i = 0; i < N; ++i) {
            fout_pos << x[i] << ", " << y[i] << "\n";
        }
        fout_pos << "\n"; // Línea en blanco entre instantes
    }

    double promT = sumaT / conteo;
    double promP = sumaP / conteo;

    fout_total << promT << ", " << promP;
    if (casoEspecial) fout_total << ", especial\n";
    else fout_total << ", vi=" << velocidad_inicial << "\n";

    fout_pos.close();
}


int main() {
    ofstream fout_total("presion_vs_temperatura.txt");
    fout_total << "# Temperatura promedio, Presion promedio, Condicion\n";

    vector<double> velocidades = {1.0, 2.0, 3.0, 4.0};
    for (double v : velocidades)
        simular(v, false, fout_total);

    simular(1.0, true, fout_total); // vx aleatoria, vy = 0

    fout_total.close();
    cout << "Simulación completa. Archivos generados:\n";
    cout << " - presion_vs_temperatura.txt\n";
    cout << " - posiciones_viX.csv por cada caso de velocidad\n";
    return 0;
}
