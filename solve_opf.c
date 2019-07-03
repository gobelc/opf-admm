
#define SIZE 4
#define NUM_BUSES SIZE
#define NUM_LINES SIZE - 1

#include<stdio.h>


int adj_matrix[NUM_BUSES][NUM_BUSES] = {{0,-1,0,0},{1,0,-1,-1},{0,1,0,0},{0,1,0,0}};

float v_bus[NUM_BUSES];
float p_inj[NUM_BUSES] = {0,15.,23.,2.};
float q_inj[NUM_BUSES] = {0,125.,23.,2.};

float P_line[NUM_LINES];
float Q_line[NUM_LINES];
float l_line[NUM_LINES];
float z_line[NUM_LINES];



int main(){
    /* inicializacion */
    for (int i = 0; i < NUM_BUSES; i++){
        v_bus[i] = 1.0;
    }

    for (int i = 0; i < NUM_LINES; i++){
        P_line[i] = p_inj[i+1];
        Q_line[i] = q_inj[i+1];
    }

    for (int i = 0; i < NUM_LINES; i++){
        for (int j = 0; j < NUM_LINES; j++){
            if (adj_matrix[i][j] == -1 ){
                P_line[i] += p_inj[j+1];
                Q_line[i] += q_inj[j+1];
            }
        }
    }


    for (int i =0; i<NUM_LINES; i++){
        if (i == 0) printf("Potencia activa en lineas:\n");
        printf("Linea %d | %f\n ", i+1, P_line[i]);
    }

    for (int i =0; i<NUM_LINES; i++){
        if (i == 0) printf("\nPotencia reactiva en lineas:\n");
        printf("Linea %d | %f\n ", i+1, Q_line[i]);
    }    



   return 0;
}