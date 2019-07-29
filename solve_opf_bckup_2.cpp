
#define SIZE 4
#define NUM_BUSES SIZE
#define NUM_LINES SIZE - 1
#define MAX_NUMBER_CHILDREN 2
#define MAX_NUMBER_NEIGHBORS MAX_NUMBER_CHILDREN + 2

#define DEBUG 1

#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <string.h>
#include <cmath>
#include <fstream>
#include <stdlib.h>

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

class Operation{

    public:
        vector<float> substract(vector<float> a,vector<float> b ){
            vector<float> result;

            if (a.size() == b.size()){
                for(int i=0;i<a.size();i++){
                    result.push_back(a[i]-b[i]);
                }

            }

            return result;
        }

        vector<float> add(vector<float> a,vector<float> b ){
            vector<float> result;

            if (a.size() == b.size()){
                for(int i=0;i<a.size();i++){
                    result.push_back(a[i]+b[i]);
                }

            }

            return result;
        }

        vector<float> multiply_scalar(vector<float> a, float multiplier ){
            vector<float> result;

            for(int i=0;i<a.size();i++){
                result.push_back(multiplier*a[i]);
            }

            return result;
        }
};



class Node: public Operation
{
    public:
        Node(int node_rank, int node_ID, int n_childs,int ancestor_ID,vector<int> childrens_ID, float R, float X, string type);
        string type;
        float rho;
        float R;
        float X;
        int node_ID;
        int node_rank;
        int n_childs;
        int ancestor_ID;
        vector<int> childrens_ID;
        vector<child_var> children_measures;
        node_var node_measures;
        vector<float> state_vector;
        vector<float> observation_vector;
        vector<float> multipliers_vector;
        vector<vector<float>> matrix;

        void set_type(string type){
            this-> type = type;
        }

        void update_state(){
            string str = "sh m/update_x.sh ";
            system((str + to_string(this->node_ID)).c_str());
            std::ifstream ifs (to_string(this->node_ID)+"/.csv", std::ifstream::in);
            char c = ifs.get();
            int counter = 0;
            while (ifs.good()) {
                this->state_vector[counter] = std::atof(&c);
                c = ifs.get();
                counter+=1;
            }
            ifs.close();
        }

        void update_observation(){
            string str = "sh m/update_y.sh ";
            system((str + to_string(this->node_ID)).c_str());
            std::ifstream ifs (to_string(this->node_ID)+"/y.csv", std::ifstream::in);
            char c = ifs.get();
            int counter = 0;
            while (ifs.good()) {
                this->observation_vector[counter] = std::atof(&c);
                c = ifs.get();
                counter+=1;
            }
            ifs.close();
        }

        void write_state_vector(){
            ofstream out(to_string(this->node_ID)+"/x.csv");
            for (float x : this->state_vector){
                out << x << endl; ;    
            }
	        out.close();
        }

        void write_observation_vector(){
            ofstream out(to_string(this-> node_ID)+"/y.csv");
            for (float x : this->observation_vector){
                out << x << endl ;    
            }
	        out.close();
        }
        void write_multipliers_vector(){
            ofstream out(to_string(this-> node_ID)+"/mu.csv");
            for (float x : this->multipliers_vector){
                    out << x << endl ;    
            }
	        out.close();
        }
        void write_matrix(){
            ofstream out(to_string(this-> node_ID)+"/A.csv");
            int CC = 7+3*n_childs;
            int RR = 3;
            int counter;
            for (vector<float> vector_temp : this->matrix){
                counter = 1;
                for (float x : vector_temp){
                    if (counter<7+3*this->n_childs){
                        out << x << "," ;    
                    } else {
                        out << x << endl; 
                    }
                    counter+=1;
                }
            }
	        out.close();
        }

        void update_multipliers(){
            std::ifstream ifs (to_string(this-> node_ID)+"/mu.csv", std::ifstream::in);
            char c = ifs.get();
            int counter = 0;
            while (ifs.good()) {
                this->multipliers_vector[counter] = std::atof(&c);
                c = ifs.get();
                counter+=1;
            }
            ifs.close();
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
                state_vector.push_back(this->children_measures[i].current);
            }

            return state_vector;
        }

        vector<float> generate_observation_vector(){
            
            vector<float> observation_vector = this->state_vector;
        
            return observation_vector;
        }

        vector<float> generate_multipliers_vector(){
            
            //vector<float> multipliers_vector = substract(this->state_vector,this->observation_vector);
            vector<float> multipliers_vector(7+3*this->n_childs,.005);
        
            return multipliers_vector;
        }


};


Node::Node(int node_rank, int node_ID int n_childs,int ancestor_ID,vector<int> childrens_ID, float R, float X, string type){
    // Constructor code

    // generate state vector
    this-> rho = .1;
    this-> R = R;
    this-> X = X;
    this-> node_rank = node_rank;
    this-> node_ID = node_ID;
    this-> n_childs = n_childs;
    this-> ancestor_ID = ancestor_ID;
    this-> childrens_ID = childrens_ID;
    this-> node_measures.voltage_ancestor=1.; 
    this-> node_measures.voltage = 1.02;
    this-> node_measures.active_power = 100;
    this-> node_measures.active_power_gen = 100;
    this-> type = type;
    for(int i=0;i<n_childs;i++){
        child_var child_var_init;
        this->children_measures.push_back(child_var_init);
    }

    // generate A Matrix
    int CC = 7+3*n_childs;
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
                tempVal = 1.; // A34
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

            int child;
            if (j>6){
                for(int c=0;c<=n_childs;c++){
                    child = abs(j-7) % 3; 
                    if(i==1 && j==7+c*3){
                        tempVal  =  1.; 
                    }
                    if(i==1 && j==9+c*3){
                        tempVal  =  -R_line[childrens_ID[c]-1]; 
                    }

                    if(i==2 && j==8+c*3){
                        tempVal  =  1.; 
                    }
                    if(i==2 && j==9+c*3){
                        tempVal  =  -X_line[childrens_ID[c]-1]; 
                    }

                }
            }

            myvector.push_back(tempVal);
        }
        matrix.push_back(myvector);
    }
    this-> matrix = matrix;

    // generate state vector
    this-> state_vector = generate_state_vector();
    this-> observation_vector = generate_observation_vector();
    this-> multipliers_vector = generate_multipliers_vector();

}


int main(int argc,char *argv[]){


    /*
    MPI init
    */
    clock_t begin = clock();
    int numtasks, rank, sendcount, recvcount, source;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    
    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    
    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, rank, numtasks);

    /* inicializacion */
    observation22 observations[NUM_BUSES];


    for(int i=0; i < NUM_BUSES;i++){
        observations[i].y = {i};
        observations[i].bus_number = i;
        print(observations[i].y);
    }

    vector<int> childrens_ID = { 4, 2, 3 ,3,1,4 };
    string type = "interior";
    Node nodo_test = Node(rank,1,6,1,childrens_ID,.1,.01,type);
    vector<float> node_test_state = nodo_test.generate_state_vector();


    /*
    vector<Node> nodos;



    for (int i=0; i<SIZE; ++i){
        string type = "interior";
        Node newNode = Node (rank,i,2,i,{ 2, 3},.1,.01,type);
        newNode.write_matrix();
        newNode.write_state_vector();
        newNode.write_observation_vector();
        newNode.write_multipliers_vector();

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
    */
    // Init 1
    for (int i = 0; i < NUM_BUSES; i++){
        nodos[i].node_measures.voltage = 1.0;
    }

    // Init 2
    for (int i = 0; i < NUM_LINES; i++){
        nodos[i].node_measures.active_power = nodos[i].node_measures.active_power_gen;
        nodos[i].node_measures.reactive_power = nodos[i].node_measures.reactive_power_gen;
    }

    // Init 3
    for (int i = 1; i < NUM_BUSES; i++){
        for (int j = 1; j < NUM_BUSES; j++){
            if (adj_matrix[i][j] == -1 ){
                float a  = P_line[i-1];
                nodos[i].node_measures.active_power += nodos[j].node_measures.active_power;
                if (DEBUG){
                    printf("i = %d, j = %d, P_line[i] = %f, P_line[j-1] = %f, P_line_new[i] = %f\n", i,j,a,nodos[j].node_measures.active_power,nodos[i].node_measures.active_power) ;
                }
                nodos[i].node_measures.reactive_power += nodos[j].node_measures.reactive_power;
            }
        }
    }

    // Init 4
    for (int i = 0; i < NUM_LINES; i++){
        nodos[i].node_measures.current = (nodos[i].node_measures.active_power* nodos[i].node_measures.active_power + nodos[i].node_measures.reactive_power*nodos[i].node_measures.reactive_power) / nodos[i].node_measures.voltage;
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


    for (int i = 0; i<NUM_LINES; i++){
        printf("Linea %d | P = %f , Q = %f, l = %f \n", i + 1, nodos[i].node_measures.active_power,nodos[i].node_measures.reactive_power,nodos[i].node_measures.current);
    }
    

    int iters = 0;
    while(iters<20){       
        std::cout << "\nIteracion: " << iters << " | Nodo: " << rank << std::endl;
        std::cout << "\nx-update... " << iters << " | Nodo: " << rank << std::endl;
        nodos[rank].update_state();
        std::cout << "\ny-update... " << iters << " | Nodo: " << rank << std::endl;
        nodos[rank].update_observation();
        std::cout << "\nmultipliers-update... " << iters << " | Nodo: " << rank << std::endl;
        nodos[rank].update_multipliers();
        iters+=1;
    }
    
    MPI_Finalize();

    return 0;


}