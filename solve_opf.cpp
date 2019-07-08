
#define SIZE 4
#define NUM_BUSES SIZE
#define NUM_LINES SIZE - 1

#include <stdio.h>
#include <math.h>

int adj_matrix[NUM_BUSES][NUM_BUSES] = {{0,-1,0,0},{1,0,-1,-1},{0,1,0,0},{0,1,0,0}};

float v_bus[NUM_BUSES];
float p_inj[NUM_BUSES] = {0,3.,2.,1};
float q_inj[NUM_BUSES] = {0,125.,23.,2.};

float P_line[NUM_LINES];
float Q_line[NUM_LINES];

float l_line[NUM_LINES];

float R_line[NUM_LINES];
float X_line[NUM_LINES];

struct state{
    float voltage;
    float p_inj;
    float q_inj;
    float P_line;
    float Q_line;
};

float x[5][NUM_BUSES];
float y[6][NUM_BUSES];

int main(){
    /* inicializacion */
    for (int i = 0; i < NUM_BUSES; i++){
        v_bus[i] = 1.0;
    }

    for (int i = 0; i < NUM_LINES; i++){
        P_line[i] = p_inj[i+1];
        Q_line[i] = q_inj[i+1];
        l_line[i] = (P_line[i]* P_line[i] + Q_line[i]*Q_line[i])/v_bus[i];
    }

    for (int i = 1; i < NUM_BUSES; i++){
        for (int j = 1; j < NUM_BUSES; j++){
            if (adj_matrix[i][j] == -1 ){
                float a  = P_line[i-1];
                P_line[i-1] += P_line[j - 1];
                printf("i = %d, j = %d, P_line[i] = %f, P_line[j-1] = %f, P_line_new[i] = %f\n", i,j,a,P_line[j-1],P_line[i-1]) ;
                Q_line[i-1] += Q_line[j - 1];
            }
        }
    }

    for (int i = 0; i < NUM_LINES; i++){
        l_line[i] = (P_line[i]* P_line[i] + Q_line[i]*Q_line[i])/v_bus[i];
    }

    for (int i = 0; i<NUM_LINES; i++){
        printf("Linea %d | P = %f , Q = %f, l = %f \n", i + 1, P_line[i],Q_line[i],l_line[i]);
    }

   return 0;
}