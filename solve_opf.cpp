
#define SIZE 4
#define NUM_BUSES SIZE
#define NUM_LINES SIZE - 1
#define MAX_NUMBER_CHILDREN 2
#define MAX_NUMBER_NEIGHBORS MAX_NUMBER_CHILDREN + 2

#define DEBUG 1

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <string.h>
#include <cmath>

using namespace std; 


int adj_matrix[NUM_BUSES][NUM_BUSES] = {{0,-1,0,0},{1,0,-1,-1},{0,1,0,0},{0,1,0,0}};
int neighbor_matrix[NUM_BUSES][MAX_NUMBER_NEIGHBORS];

float v_bus[NUM_BUSES] = {1.,1.,1.,1.};
float p_inj[NUM_BUSES] = {2.34,3.12,2.2,1.};
float q_inj[NUM_BUSES] = {.32,.13,.23,.35};

float P_line[NUM_LINES];
float Q_line[NUM_LINES];

float l_line[NUM_LINES];

float R_line[NUM_LINES]= {.11111,.22222,.33333};
float X_line[NUM_LINES]= {.12121,.21212,.31313};
 
struct state{
    float voltage;
    float p_inj;
    float q_inj;
    float P_line;
    float Q_line;
};

float x[NUM_BUSES][6+1+4*MAX_NUMBER_CHILDREN]; // xi = [vi_x,pi_x,qi_x,Pi_x,Qi_x,li_x]
float y[NUM_BUSES][6][MAX_NUMBER_NEIGHBORS]; // yij = [vij_y,pij_y,qij_y,Pij_y,Qij_y,lij_y] \ TRIDIMENSIONAL ARRAY: NODE * VARIABLE * NEIGHBOR.

// Create observation vectors containing floats

struct observation22{
    vector<int> y;
    int bus_number;
};

void print(vector <int> const &a) {
   cout << "The vector elements are : ";
   
   for(int i=0; i < a.size(); i++)
      cout << a.at(i) << ' ';
}


struct child_var{
    int node_ID;
    string type;
    float active_power_gen;
    float reactive_power_gen;
    float active_power;
    float reactive_power;
    float current;
};

struct node_var{
    float voltage_ancestor;
    float voltage;
    float active_power_gen;
    float reactive_power_gen;
    float active_power;
    float reactive_power;
    float current;
};


class Node
{
    public:
        Node(int n_childs,int ancestor_ID,vector<int> childrens_ID, float R, float X, string type);
        string type;
        float R;
        float X;
        int node_ID;
        int n_childs;
        int ancestor_ID;
        vector<int> childrens_ID;
        vector<child_var> children_measures;
        node_var node_measures;
        vector<float> state_vector;
        vector<vector<float>> matrix;

        void set_type(string type){
            this-> type = type;
        }

        void update_state(){
        }

        void update_observation(){

        }

        void update_multipliers(){

        }

        vector<float> generate_state_vector(){
            
            vector<float> state_vector;
            
            state_vector.push_back(this->node_measures.voltage_ancestor);
            state_vector.push_back(this->node_measures.voltage);
            state_vector.push_back(this->node_measures.active_power_gen);
            state_vector.push_back(this->node_measures.reactive_power_gen);
            state_vector.push_back(this->node_measures.active_power);
            state_vector.push_back(this->node_measures.reactive_power);
            state_vector.push_back(this->node_measures.current);
            
            for(int i=0;i<this->n_childs;i++){
                state_vector.push_back(this->children_measures[i].active_power);
                state_vector.push_back(this->children_measures[i].reactive_power);
                state_vector.push_back(this->children_measures[i].active_power_gen);
                state_vector.push_back(this->children_measures[i].reactive_power_gen);
                state_vector.push_back(this->children_measures[i].current);
            }

            return state_vector;
        }

};


Node::Node(int n_childs,int ancestor_ID,vector<int> childrens_ID, float R, float X, string type){
    // Constructor code

    // generate state vector
    this-> R = R;
    this-> X = X;
    this-> n_childs = n_childs;
    this-> ancestor_ID = ancestor_ID;
    this-> childrens_ID = childrens_ID;
    this-> node_measures.voltage_ancestor=1.; 
    this-> type = type;
    for(int i=0;i<n_childs;i++){
        child_var child_var_init;
        this->children_measures.push_back(child_var_init);
    }

    // generate A Matrix
    int CC = 7+5*n_childs;
    int RR = 3;
    vector<vector<float>> matrix;

    //cin>>CC; cin>>RR; already done
    for(int i = 0; i<RR; i++)
    {
        vector<float> myvector;
        for(int j = 0; j<CC; j++)
        {
            float tempVal = 0.;
            if (i==0 && j==0){
                tempVal = 1.; //A11
            }
            if (i==0 && j==1){
                tempVal = -1.; //A12
            }
            if (i==1 && j==2){
                tempVal = 1.; // A23
            }
            if (i==2 && j==3){
                tempVal = -1.; // A34
            }
            if (i==1 && j==4){
                tempVal = -1.; // A25
            }
            if (i==2 && j==5){
                tempVal = -1.; // A36
            }
            if (i==0 && j==4){
                tempVal = 2*R; // A15
            }
            if (i==0 && j==5){
                tempVal = 2*X; // A16
            }
            if (i==0 && j==6){
                tempVal = -(pow(X,2)+pow(R,2)) ; // A17
            }
            if (i==0 && j==6){
                tempVal = -(pow(X,2)+pow(R,2)) ; // A17
            }
            int child;
            if (j>6){
                for(int c=0;c<=n_childs;c++){
                    child = abs(j-7) % 5; 
                    if(i==1 && j==7+c*5){
                        tempVal  =  1.; 
                    }
                    if(i==1 && j==9+c*5){
                        tempVal  =  -R_line[childrens_ID[c]-1]; 
                    }
                    if(i==1 && j==10+c*5){
                        tempVal  =  1; 
                    }
                    if(i==2 && j==8+c*5){
                        tempVal  =  1.; 
                    }
                    if(i==2 && j==9+c*5){
                        tempVal  =  -X_line[childrens_ID[c]-1]; 
                    }
                    if(i==2 && j==11+c*5){
                        tempVal  =  1.; 
                    }
                }
            }

            myvector.push_back(tempVal);
        }
        matrix.push_back(myvector);
    }
    this-> matrix = matrix;
}


int main(){
    /* inicializacion */

    // Init 1
    for (int i = 0; i < NUM_BUSES; i++){
        v_bus[i] = 1.0;
    }

    // Init 2
    for (int i = 0; i < NUM_LINES; i++){
        P_line[i] = p_inj[i+1];
        Q_line[i] = q_inj[i+1];
    }

    // Init 3
    for (int i = 1; i < NUM_BUSES; i++){
        for (int j = 1; j < NUM_BUSES; j++){
            if (adj_matrix[i][j] == -1 ){
                float a  = P_line[i-1];
                P_line[i-1] += P_line[j - 1];
                if (DEBUG){
                    printf("i = %d, j = %d, P_line[i] = %f, P_line[j-1] = %f, P_line_new[i] = %f\n", i,j,a,P_line[j-1],P_line[i-1]) ;
                }
                Q_line[i-1] += Q_line[j - 1];
            }
        }
    }

    // Init 4
    for (int i = 0; i < NUM_LINES; i++){
        l_line[i] = (P_line[i]* P_line[i] + Q_line[i]*Q_line[i])/v_bus[i];
    }

    // Make neighbor matrix
    int c;
    neighbor_matrix[0][0] = -1; // tree root
    for (int i = 0;  i < NUM_BUSES; i++){
        c = 2;
            for (int k = 0;  k < NUM_BUSES; k++){ 
                //search for parent
                if (adj_matrix[i][k] > 0){
                    neighbor_matrix[i][0] = k;
                }
                //search for children
                if (adj_matrix[i][k] < 0){
                    neighbor_matrix[i][c] = k;
                    c = c + 1;
                }
            }
        neighbor_matrix[i][1] = i; // self
    }

    if (DEBUG){
        printf("\nNeighbors matrix:\n"); 
        for (int i = 0;  i < NUM_BUSES; i++){
            printf("Node %d: [%d,%d,%d,%d,%d]\n",i,neighbor_matrix[i][0],neighbor_matrix[i][1],neighbor_matrix[i][2],neighbor_matrix[i][3],neighbor_matrix[i][4]); 
        }
    }

    // Init state values

    // Root node
    x[0][0] = v_bus[0];
    x[0][1] = p_inj[0];
    x[0][2] = q_inj[0];
    x[0][3] = 0;
    x[0][4] = 0;
    x[0][5] = 0;
    
    // Non-root nodes
    for (int i = 1; i < NUM_BUSES; i++){
        //[vi_x,pi_x,qi_x,Pi_x,Qi_x,li_x]
        x[i][0] = v_bus[i];
        x[i][1] = p_inj[i];
        x[i][2] = q_inj[i];
        x[i][3] = P_line[i-1];
        x[i][4] = Q_line[i-1];
        x[i][5] = l_line[i-1];
        x[i][6] = P_line[i-1];
        x[i][7] = Q_line[i-1];
        x[i][8] = l_line[i-1];
        x[i][9] = l_line[i-1];

    }
    
    // Init 5
    int neighbor;
    for (int i = 0;  i < NUM_BUSES; i++){
        for (int m = 0; m < 6; m++){
            for (int j = 0; j < MAX_NUMBER_NEIGHBORS; j++){
                y[i][m][j] = 0; // init observation matrix as zeros
                neighbor = neighbor_matrix[i][j];
                if (j < 1) {
                    y[i][0][j]=x[neighbor][0]; // yij = vi_x
                }
    
                if (j > 1){ // children
                    y[i][3][j]=x[neighbor][3]; // yij_3 = Pi_x
                    y[i][4][j]=x[neighbor][4]; // yij_4 = Qi_x
                    y[i][5][j]=x[neighbor][5]; // yij_5 = li_x
                }
    
                if (j == 1){ // self
                    y[i][0][j]=x[neighbor][0]; // yij_0 = vi_x
                    y[i][1][j]=x[neighbor][1]; // yij_1 = pi_x
                    y[i][2][j]=x[neighbor][2]; // yij_2 = qi_x
                    y[i][3][j]=x[neighbor][3]; // yij_3 = Pi_x
                    y[i][4][j]=x[neighbor][4]; // yij_4 = Qi_x
                    y[i][5][j]=x[neighbor][5]; // yij_5 = li_x
                }
            }   
        }
    }


    for (int i = 0; i<NUM_LINES; i++){
        printf("Linea %d | P = %f , Q = %f, l = %f \n", i + 1, P_line[i],Q_line[i],l_line[i]);
    }

    if (DEBUG){
        printf("\nState matrix:\n"); 
        for (int i = 0;  i < NUM_BUSES; i++){
            printf("Node %d: [v = %f,p = %f,q = %f,P = %f,Q = %f,l = %f]\n",i,x[i][0],x[i][1],x[i][2],x[i][3],x[i][4],x[i][5]); 
            
        }
    }

    if (DEBUG){
        printf("\nObservation matrix:\n"); 
        for (int i = 0;  i < NUM_BUSES; i++){
            for (int j = 0; j < MAX_NUMBER_NEIGHBORS; j++){
                //printf("Node %d - Neighbor %d: [v = %f,p = %f,q = %f,P = %f,Q = %f,l = %f]\n",i,j,y[i][0][j],y[i][1][j],y[i][2][j],y[i][3][j],y[i][4][j],y[i][5][j]); 
                printf("y%d%d=[%f,%f,%f,%f,%f,%f]\n",i,j,y[i][0][j],y[i][1][j],y[i][2][j],y[i][3][j],y[i][4][j],y[i][5][j]); 
            }
        }
    }


    observation22 observations[NUM_BUSES];


    for(int i=0; i < NUM_BUSES;i++){
        observations[i].y = {i};
        observations[i].bus_number = i;
        print(observations[i].y);
    }

    vector<int> childrens_ID = { 4, 2, 3 ,3,1,4 };
    string type = "interior";
    Node nodo_test = Node(6,1,childrens_ID,.1,.01,type);
    vector<float> node_test_state = nodo_test.generate_state_vector();

    vector<Node> nodos;

    for (int i=0; i<SIZE; ++i){
        string type = "interior";
        Node newNode = Node (2,i,{ 2, 3},.1,.01,type);
        cout << "Node " << i << endl;
        int j = 1;

        for (int x : newNode.generate_state_vector()){
            cout << "Medida " << j << ": "  << x << " \n";
            j+=1;
        }

        int columna = 0;
        int fila = 0;
        for (vector<float> vector_temp : newNode.matrix){
            columna = 0;
            for (float x : vector_temp){
                cout << "Matriz " <<  fila << columna << ": " << x << " \n";
                columna+=1;
            }
            fila+=1;
        }

        
        nodos.push_back(newNode);
    }

    return 0;


}