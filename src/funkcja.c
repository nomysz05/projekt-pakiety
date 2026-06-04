#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void mnozenie_4x4(double* a, double* b, int n1, int m1, int m2, double* c) {
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < m2; j++) {
            int c_index = i * m2 + j;
            c[c_index] = 0.0;
            for (int k = 0; k < m1; k++) {
                c[c_index] += a[i * m1 + k] * b[k * m2 + j];
            }
        }
    }
}
int invert_4x4(double *m, double *inv) {
    double temp[16];
    for (int i = 0; i < 16; i++) {
        temp[i] = m[i];
        inv[i] = (i % 5 == 0) ? 1.0 : 0.0; 
    }

    for (int i = 0; i < 4; i++) {
        int pivot = i;
        for (int j = i + 1; j < 4; j++) {
            if (fabs(temp[j * 4 + i]) > fabs(temp[pivot * 4 + i])) pivot = j;
        }
        
        if (fabs(temp[pivot * 4 + i]) < 1e-12) return 0; 

        if (pivot != i) {
            for (int k = 0; k < 4; k++) {
                double t = temp[i * 4 + k]; temp[i * 4 + k] = temp[pivot * 4 + k]; temp[pivot * 4 + k] = t;
                t = inv[i * 4 + k]; inv[i * 4 + k] = inv[pivot * 4 + k]; inv[pivot * 4 + k] = t;
            }
        }

        double factor = temp[i * 4 + i];
        for (int k = 0; k < 4; k++) {
            temp[i * 4 + k] /= factor;
            inv[i * 4 + k] /= factor;
        }

        for (int j = 0; j < 4; j++) {
            if (j != i) {
                double f = temp[j * 4 + i];
                for (int k = 0; k < 4; k++) {
                    temp[j * 4 + k] -= f * temp[i * 4 + k];
                    inv[j * 4 + k] -= f * inv[i * 4 + k];
                }
            }
        }
    }
    return 1;
}

__attribute__ ((visibility ("default")))
SEXP kalman(SEXP r_time, SEXP r_lat, SEXP r_lon, SEXP r_speed, SEXP r_sigma_a, SEXP r_sigma_gps) {
    double *time = REAL(r_time);
    double *lat = REAL(r_lat);
    double *lon = REAL(r_lon);
    double *speed = REAL(r_speed);
    size_t n = XLENGTH(r_lat);
    
    double sigma_a = Rf_asReal(r_sigma_a);
    double sigma_gps = Rf_asReal(r_sigma_gps);

    SEXP r_new_lat = Rf_allocVector(REALSXP, n);
    SEXP r_new_lon = Rf_allocVector(REALSXP, n);
    PROTECT(r_new_lat);
    PROTECT(r_new_lon);
    
    double *new_lat = REAL(r_new_lat);
    double *new_lon = REAL(r_new_lon);

    new_lat[0] = lat[0];
    new_lon[0] = lon[0];
    
    double lat0 = lat[0];
    double lon0 = lon[0];

    double x_est[4] = { 0.0, 0.0, 0.0, 0.0 };

    double P_est[16] = {
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    };
    double A[16] = {
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    };
    double AT[16] = {
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    };

    for (size_t i = 1; i < n; i++) {
        double dt = time[i] - time[i-1];
        if (dt <= 0.0) dt = 1.0;

        double z_lat_m = (lat[i] - lat0) * 111111.0;
        double z_lon_m = (lon[i] - lon0) * 68000.0;
        double delta_lat = lat[i] - lat[i-1];
        double delta_lon = lon[i] - lon[i-1];
        double theta = atan2(delta_lon, delta_lat);
        
        double vx = speed[i] * cos(theta)/3.6;
        double vy = speed[i] * sin(theta)/3.6;
        
        if (i == 1) {
            x_est[2] = vx;
            x_est[3] = vy;
        }

        A[2]  = dt;
        A[7]  = dt; 
        AT[8] = dt; 
        AT[13] = dt; 
        
        double sigma_a_sq = sigma_a * sigma_a;
        double dt2 = dt * dt;
        double dt3 = dt2 * dt;
        double dt4 = dt3 * dt;
        
        double qt4 = (dt4 / 4.0) * sigma_a_sq;
        double qt3 = (dt3 / 2.0) * sigma_a_sq;
        double qt2 = dt2 * sigma_a_sq;
        
        double Q[16] = {
            qt4, 0.0, qt3, 0.0,
            0.0, qt4, 0.0, qt3,
            qt3, 0.0, qt2, 0.0,
            0.0, qt3, 0.0, qt2
        };
        
        double x_pred[4];
        x_pred[0] = x_est[0] + dt * x_est[2];
        x_pred[1] = x_est[1] + dt * x_est[3];
        x_pred[2] = x_est[2];                 
        x_pred[3] = x_est[3];       
        
        double P_pred[16];
        double temp[16]; 
        mnozenie_4x4(A, P_est, 4, 4, 4, temp);
        mnozenie_4x4(temp, AT, 4, 4, 4, P_pred);

        for (int j = 0; j < 16; j++) {
            P_pred[j] += Q[j];
        }     
        
        double z[4] = { z_lat_m, z_lon_m, vx, vy };
        double y[4];
        for (int j = 0; j < 4; j++) {
            y[j] = z[j] - x_pred[j];
        }
        
        double distance_error_sq = (y[0] * y[0]) + (y[1] * y[1]);
        
        if (i > 1 && distance_error_sq > (30.0 * 30.0)) {
            x_est[0] = x_pred[0];
            x_est[1] = x_pred[1];
            x_est[2] = x_pred[2];
            x_est[3] = x_pred[3];
            
            for (int j = 0; j < 16; j++) P_est[j] = P_pred[j];

            P_est[0]  += 25.0;
            P_est[5]  += 25.0;
            
            new_lat[i] = lat0 + (x_est[0] / 111111.0);
            new_lon[i] = lon0 + (x_est[1] / 68000.0);
            continue;
        }

        double r_pos = sigma_gps * sigma_gps;
        double r_vel = 0.5;
        double R[16] = {
            r_pos, 0.0,   0.0,   0.0,
            0.0,   r_pos, 0.0,   0.0,
            0.0,   0.0,   r_vel, 0.0,
            0.0,   0.0,   0.0,   r_vel
        };
        
        double S[16];
        for (int j = 0; j < 16; j++) S[j] = P_pred[j] + R[j];
        
        double S_inv[16];
        if (!invert_4x4(S, S_inv)) {
            for (int j = 0; j < 16; j++) S_inv[j] = (j % 5 == 0) ? 1.0 : 0.0; 
        }

        double K[16];
        mnozenie_4x4(P_pred, S_inv, 4, 4, 4, K);

        for (int j = 0; j < 4; j++) {
            double ky_sum = 0.0;
            for (int k = 0; k < 4; k++) {
                ky_sum += K[j * 4 + k] * y[k];
            }
            x_est[j] = x_pred[j] + ky_sum;
        }

        double I_minus_K[16];
        for (int j = 0; j < 16; j++) {
            double I_elem = (j % 5 == 0) ? 1.0 : 0.0;
            I_minus_K[j] = I_elem - K[j];
        }
        mnozenie_4x4(I_minus_K, P_pred, 4, 4, 4, P_est);

        new_lat[i] = lat0 + (x_est[0] / 111111.0);
        new_lon[i] = lon0 + (x_est[1] / 68000.0);
    }
    SEXP r_output_list = Rf_allocVector(VECSXP, 2);
    PROTECT(r_output_list);
    SET_VECTOR_ELT(r_output_list, 0, r_new_lat);
    SET_VECTOR_ELT(r_output_list, 1, r_new_lon);

    UNPROTECT(3); 
    return r_output_list;
}
