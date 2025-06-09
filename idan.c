#include <stdio.h>

#define N 400  // Jumlah langkah (simulasi 4 detik dengan dt = 0.01)
#define DT 0.01
#define M 1.0
#define K 10.0

// Data diskrit gaya
double t_data[] = {0, 1, 2, 3, 4};
double f_data[] = {0, 10, 15, 10, 0};
int f_len = 5;

// Interpolasi linear gaya F(t)
double interpolate_force(double t) {
    for (int i = 0; i < f_len - 1; i++) {
        if (t >= t_data[i] && t <= t_data[i+1]) {
            double t0 = t_data[i];
            double t1 = t_data[i+1];
            double f0 = f_data[i];
            double f1 = f_data[i+1];
            return f0 + (f1 - f0) * (t - t0) / (t1 - t0);
        }
    }
    return 0.0;
}

// Turunan pertama (dx/dt) dan kedua (dv/dt)
void derivatives(double t, double x, double v, double* dxdt, double* dvdt) {
    double F = interpolate_force(t);
    *dxdt = v;
    *dvdt = (F - K * x) / M;
}

int main() {
    double t = 0.0, x = 0.0, v = 0.0;
    double energy = 0.0;

    FILE *fout = fopen("result.csv", "w");
    fprintf(fout, "t,x,v,Ek\n");

    for (int i = 0; i < N; i++) {
        double k1x, k1v, k2x, k2v, k3x, k3v, k4x, k4v;

        derivatives(t, x, v, &k1x, &k1v);
        derivatives(t + DT/2, x + DT/2 * k1x, v + DT/2 * k1v, &k2x, &k2v);
        derivatives(t + DT/2, x + DT/2 * k2x, v + DT/2 * k2v, &k3x, &k3v);
        derivatives(t + DT, x + DT * k3x, v + DT * k3v, &k4x, &k4v);

        x += (DT/6)*(k1x + 2*k2x + 2*k3x + k4x);
        v += (DT/6)*(k1v + 2*k2v + 2*k3v + k4v);
        t += DT;

        double Ek = 0.5 * M * v * v;
        energy += Ek * DT;  // Trapezoidal approx (tanpa rata-rata: sederhana)

        fprintf(fout, "%.4f,%.6f,%.6f,%.6f\n", t, x, v, Ek);
    }

    fclose(fout);
    printf("Simulasi selesai. Total energi kinetik terintegrasi: %.4f J\n", energy);
    return 0;
}