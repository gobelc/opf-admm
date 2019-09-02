#include <unistd.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <string.h>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include "NetworkModel_8bus.hpp"

#define DEBUG 1
#define SYNC 0

/// ADMM Parameters
#define RHO 1.
#define MAX_ITER 100


/// NETWORK (THIS MUST BE A COPY OF THE NETWORK MODEL PARAMATERS)
#define SIZE 8
#define NUM_BUSES SIZE
#define NUM_LINES SIZE - 1
#define MAX_NUMBER_CHILDREN 2
#define MAX_NUMBER_NEIGHBORS MAX_NUMBER_CHILDREN + 2

using namespace std;


extern string path;
extern int adj_matrix[NUM_BUSES][NUM_BUSES];
extern float v_bus[NUM_BUSES];
extern float p_inj[NUM_BUSES];
extern float q_inj[NUM_BUSES];
extern float R_line[NUM_BUSES-1];
extern float X_line[NUM_BUSES-1];
extern float S_base;

float P_line[NUM_LINES];
float Q_line[NUM_LINES];
float l_line[NUM_LINES];

char buffer_logFile [100];


int neighbor_matrix[NUM_BUSES][MAX_NUMBER_NEIGHBORS];
//float rho = RHO;
// Create observation vectors containing floats

/*struct observation22{
    vector<int> y;
    int bus_number;
};
*/

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

struct observation_var{
    float voltage_ancestor;
    float voltage;
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
        string filePath;
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
        observation_var node_observations;

        vector<float> state_vector;
        vector<float> observation_vector;
        vector<float> multipliers_vector;
        vector<vector<float>> matrix;

        void set_type(string type){
            this-> type = type;
        }

        void update_state(){
            
            std::ifstream ifs1 (to_string(this->node_ID)+"/rho.csv", std::ifstream::in);
            char c_rho[10];
            ifs1.getline(c_rho,10);
            ifs1.close();
            rho = std::atof(c_rho);
            
            cout << fixed << setprecision(5) << rho <<endl;
            string str = "sh m/update_x.sh ";
            system((str + to_string(this->node_ID) + " " + to_string(rho) + " " + path).c_str());
            std::ifstream ifs (to_string(this->node_ID)+"/x.csv", std::ifstream::in);
            char c[10];
            ifs.getline(c,10);
            int counter = 0;
            float measure;
            while (ifs.good()) {
                measure = std::atof(c);
                this->state_vector[counter] = measure;
                if (counter==0){ this->node_measures.voltage_ancestor = measure;}
                if (counter==1){ this->node_measures.voltage = measure;}
                if (counter==2){ this->node_measures.active_power_gen = measure;}
                if (counter==3){ this->node_measures.reactive_power_gen = measure;}
                if (counter==4){ this->node_measures.active_power = measure;}
                if (counter==5){ this->node_measures.reactive_power = measure;}
                if (counter==6){ this->node_measures.current = measure;}
                if (counter==7){ this->children_measures[0].active_power = measure;}
                if (counter==8){ this->children_measures[0].reactive_power = measure;}
                if (counter==9){ this->children_measures[0].current = measure;}
                if (counter==10){ this->children_measures[1].active_power = measure;}
                if (counter==11){ this->children_measures[1].reactive_power = measure;}
                if (counter==12){ this->children_measures[1].current = measure;}
                ifs.getline(c,10);
                counter+=1;
            }
            ifs.close();
        }

        void update_observation(){
            
            std::ifstream ifs2 (to_string(this->node_ID)+"/rho.csv", std::ifstream::in);
            char c_rho[10];
            ifs2.getline(c_rho,10);
            ifs2.close();
            rho = std::atof(c_rho);
            
            string str = "sh m/update_y.sh ";
            system((str + to_string(this->node_ID) + " " + to_string(rho) + " " + path).c_str());
            std::ifstream ifs (to_string(this->node_ID)+"/y.csv", std::ifstream::in);
            char c[10];
            ifs.getline(c,10);
            int counter = 0;
            float measure;
            n_childs = this->n_childs;
            while (ifs.good()) {
                measure = std::atof(c);
                this->observation_vector[counter] = measure;
                if (counter==0){ this->node_observations.voltage_ancestor = measure;}
                if (counter==1){ this->node_observations.voltage = measure;}
                if (counter==2){ this->node_observations.active_power_gen = measure;}
                if (counter==3){ this->node_observations.reactive_power_gen = measure;}
                if (counter==4){ this->node_observations.active_power = measure;}
                if (counter==5){ this->node_observations.reactive_power = measure;}
                if (counter==6){ this->node_observations.current = measure;}
                if (counter==7){ this->children_measures[0].active_power = measure;}
                if (counter==8){ this->children_measures[0].reactive_power = measure;}
                if (counter==9){ this->children_measures[0].current = measure;}
                if (counter==10){ this->children_measures[1].active_power = measure;}
                if (counter==11){ this->children_measures[1].reactive_power = measure;}
                if (counter==12){ this->children_measures[1].current = measure;}
                
                ifs.getline(c,10);
                counter+=1;
            }
            ifs.close();
        }

        void write_state_vector(){
            ofstream out(to_string(this->node_ID)+"/x.csv");
            for (float x : this->state_vector){
                out << x << endl;    
            }
	        out.close();
        }

        void write_observation_vector(){
            ofstream out(to_string(this-> node_ID)+"/y.csv");
            for (float x : this->observation_vector){
                out <<  x << endl;    
            }
	        out.close();
        }
        void write_multipliers_vector(){
            ofstream out(to_string(this-> node_ID)+"/mu.csv");
            for (float x : this->multipliers_vector){
                    out << fixed << setprecision(5) << x << endl ;    
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
            char c[10];
            ifs.getline(c,10);
            int counter = 0;
            while (ifs.good()) {
                //cout << std::atof(c) << endl;
                this->multipliers_vector[counter] = std::atof(c);
                ifs.getline(c,10);
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


        vector<float> init_observation_vector(){
            
            vector<float> observation_vector = this->state_vector;

            return observation_vector;
        }

        vector<float> generate_observation_vector(){
            
            vector<float> observation_vector;

            observation_vector.push_back(this->node_observations.voltage_ancestor);
            observation_vector.push_back(this->node_observations.voltage);
            observation_vector.push_back(this->node_observations.active_power_gen);
            observation_vector.push_back(this->node_observations.reactive_power_gen);
            observation_vector.push_back(this->node_observations.active_power);
            observation_vector.push_back(this->node_observations.reactive_power);
            observation_vector.push_back(this->node_observations.current);
            
            for(int i=0;i<this->n_childs;i++){
                observation_vector.push_back(this->children_measures[i].active_power);
                observation_vector.push_back(this->children_measures[i].reactive_power);
                observation_vector.push_back(this->children_measures[i].current);
            }    

            return observation_vector;
        }

        vector<float> generate_multipliers_vector(){
            
            //vector<float> multipliers_vector = substract(this->state_vector,this->observation_vector);
            vector<float> multipliers_vector(7+3*this->n_childs,.005);
        
            return multipliers_vector;
        }

        void log_screen(std::string message){
            string now = getCurrentDateTime("now");
            std::cout << now << '\t' << "| Nodo " << this-> node_rank << "> " << message <<  std::endl;
        }

        inline void log_file(std::string logMsg){
            string now = getCurrentDateTime("now");
            ofstream ofs(this->filePath.c_str(), std::ios_base::out | std::ios_base::app );
            ofs << now << '\t' << logMsg << '\n';
            ofs.close();
        }

        inline string getCurrentDateTime( string s ){
            time_t now = time(0);
            struct tm  tstruct;
            char  buf[80];
            tstruct = *localtime(&now);
            if(s=="now")
                strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
            else if(s=="date")
                strftime(buf, sizeof(buf), "%Y-%m-%d", &tstruct);
            return string(buf);
        };

};


Node::Node(int node_rank, int node_ID, int n_childs,int ancestor_ID,vector<int> childrens_ID, float R, float X, string type){
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
    this-> filePath = path + "/" + to_string(this->node_ID) + "/" +  getCurrentDateTime("date")+".txt";
    
    child_var child_var_init;
    child_var_init.current = 0.;
    child_var_init.active_power = 0.;
    child_var_init.reactive_power = 0.;


    int root = 0;
    for(int i=0;i<n_childs;i++){
        this->children_measures.push_back(child_var_init);
    }

    // generate A Matrix
    int CC = 7+3*n_childs;
    int RR = 3;
    vector<vector<float>> matrix;
    float R_anc, X_anc, g, b = 0;

    /*
    if (node_ID != root){
        R_anc = R_line[ancestor_ID];
        X_anc = X_line[ancestor_ID];
        g=R_anc/(pow(X_anc,2)+pow(R_anc,2));
        b=X_anc/(pow(X_anc,2)+pow(R_anc,2));
    } // shunt impedance
    */
    

    //cin>>CC; cin>>RR; already done
    for(int i = 0; i<RR; i++)
    {
        vector<float> myvector;
        for(int j = 0; j<CC; j++)
        {
            float tempVal = 0.;
            if (i==0 && j==0){
                tempVal = 1.; //A11
                if (node_ID == root) tempVal=0;
            }
            if (i==0 && j==1){
                tempVal = -1.; //A12
            }
            if (i==1 && j==0){
                tempVal = -g; //A21
                if (node_ID == root) tempVal=0;
            }
            if (i==1 && j==2){
                tempVal = 1.; // A23
            }
            if (i==2 && j==0){
                tempVal = -b; //A21
                if (node_ID == root) tempVal=0;

            }
            if (i==2 && j==3){
                tempVal = 1.; // A34
            }
            if (i==1 && j==4){
                tempVal = -1.; // A25
                if (node_ID == root) tempVal=0;

            }
            if (i==2 && j==5){
                tempVal = -1.; // A36
                if (node_ID == root) tempVal=0;

            }
            if (i==0 && j==4){
                tempVal = 2*R; // A15
            }
            if (i==0 && j==5){
                tempVal = 2*X; // A16
            }
            if (i==0 && j==6){
                tempVal = -(pow(X,2)+pow(R,2)) ; // A17
                if (node_ID == root) tempVal=0;
            }

            if (j>6){
                int c = abs(j-7) / 3 + 1; 
                cout << "Nodo: " << this->node_ID <<"Hijo: " <<  c << ", i= " << i << ", j= " << j <<endl;
                if(i==1 && j==4+c*3){
                    tempVal  =  1.; 
                }
                if(i==1 && j==6+c*3){
                    tempVal  =  -R_line[childrens_ID[c-1]-1]; 
                    cout << "Impedancia R linea hijo:" << tempVal <<endl;
                }

                if(i==2 && j==5+c*3){
                    tempVal  =  1.; 
                }
                if(i==2 && j==6+c*3){
                    tempVal  = -X_line[childrens_ID[c-1]-1]; 
                    cout << "Impedancia Q linea hijo:" << tempVal<<endl;

                }
            }

            myvector.push_back(tempVal);
        }
        matrix.push_back(myvector);
    }
    this-> matrix = matrix;

    // generate state vector
    this-> state_vector = generate_state_vector();
    this-> observation_vector = init_observation_vector();
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

    /* inicializacion 
    observation22 observations[NUM_BUSES];


    for(int i=0; i < NUM_BUSES;i++){
        observations[i].y = {i};
        observations[i].bus_number = i;
        print(observations[i].y);
    }
    */
    
    // Make neighbor matrix
    int c;
    neighbor_matrix[0][0] = -1; // tree root
    //for (int i = 0;  i < NUM_BUSES; i++){
    c = 2;
    for (int k = 0;  k < NUM_BUSES; k++){ 
        //search for parent
        if (adj_matrix[rank][k] > 0){
            neighbor_matrix[rank][0] = k;
        }
        //search for children
        if (adj_matrix[rank][k] < 0){
            neighbor_matrix[rank][c] = k;
            c = c + 1;
        }
    }
    neighbor_matrix[rank][1] = rank; // self

    // Node creation
    string type = "interior";
    int i = rank;

    int n_childs = 0;
    vector<int> childrens_ID;
    for(int i=2;i<MAX_NUMBER_NEIGHBORS;i++){
        if (neighbor_matrix[rank][i] > 0){
            n_childs +=1;
            childrens_ID.push_back(neighbor_matrix[rank][i]);
        }
    }
    int ancestor_ID = neighbor_matrix[rank][0];
    Node nodo = Node(rank,rank,n_childs,ancestor_ID,childrens_ID, i>0? R_line[i-1]:.0,i>0? X_line[i-1]:.0,type);

    nodo.write_matrix();
    nodo.write_state_vector();
    nodo.write_observation_vector();
    nodo.write_multipliers_vector();

    cout << "Node " << i << endl;
    int j = 1;

    for (int x : nodo.generate_state_vector()){
        cout << "Medida " << j << ": "  << x << " \n";
        j+=1;
    }

    int columna = 0;
    int fila = 0;
    for (vector<float> vector_temp : nodo.matrix){
        columna = 0;
        for (float x : vector_temp){
            cout << "Matriz " <<  fila << columna << ": " << x << " \n";
            columna+=1;
        }
        fila+=1;
    }

    nodo.node_measures.voltage = 1.0;
    nodo.node_measures.active_power = nodo.node_measures.active_power_gen;
    nodo.node_measures.reactive_power = nodo.node_measures.reactive_power_gen;
    nodo.node_measures.current = (nodo.node_measures.active_power*nodo.node_measures.active_power + nodo.node_measures.reactive_power*nodo.node_measures.reactive_power) / nodo.node_measures.voltage;


    if (DEBUG){
        printf("\nNeighbors matrix:\n"); 
        printf("Node %d: [%d,%d,%d,%d,%d]\n",rank,neighbor_matrix[rank][0],neighbor_matrix[rank][1],neighbor_matrix[rank][2],neighbor_matrix[rank][3],neighbor_matrix[rank][4]); 
    }

    printf("Linea %d | P = %f , Q = %f, l = %f \n", rank + 1, nodo.node_measures.active_power,nodo.node_measures.reactive_power,nodo.node_measures.current);
    
    int iters = 0;

    while(iters<MAX_ITER){

        // MESSAGE PASSING
        struct children_observations{
            int nodeID;
            float current_obs;
            float active_power_obs;
            float reactive_power_obs;
        };


        float voltage_obs = nodo.node_observations.voltage;
        float current_obs = nodo.node_observations.current;
        float active_power_obs = nodo.node_observations.active_power;
        float reactive_power_obs = nodo.node_observations.reactive_power;

        int res = MPI_Barrier(MPI_COMM_WORLD);
        if (SYNC){
            int res = MPI_Barrier(MPI_COMM_WORLD);
            if (res != MPI_SUCCESS){
                fprintf (stderr, "MPI_Barrier failed\n");
                exit (0);
            } else {
                cout << "Node " << nodo.node_ID << " attained the sync barrier.\n" << endl;
            }
        }

        nodo.update_state();
        nodo.update_observation();

        MPI_Request req_v;
        MPI_Status status_v;

        MPI_Request req_current[SIZE];
        MPI_Status status_current[SIZE];
        
        MPI_Request req_power[SIZE];
        MPI_Status status_power[SIZE];
        
        MPI_Request req_reactive[SIZE];
        MPI_Status status_reactive[SIZE];

        if (SYNC && nodo.node_ID ==0){
            int res = MPI_Barrier(MPI_COMM_WORLD);
            if (res != MPI_SUCCESS){
                fprintf (stderr, "MPI_Barrier failed\n");
            } 
        }

        cout << "Node " << nodo.node_ID << " passed the sync barrier after update.\n" << endl;
        
        //Send measures to children
        for(int i=0; i<n_childs;i++){
            sprintf(buffer_logFile, "Node %d is sending voltage to children %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.childrens_ID[i],voltage_obs, iters);
            nodo.log_file(buffer_logFile);
            //cout << std::fixed << std::setprecision(3) << "Node " << nodo.node_ID << " is sending voltage to children " << nodo.childrens_ID[i]  << ". Measure: " << voltage_obs << endl;
            if (SYNC){
                MPI_Send(&voltage_obs, 1, MPI_FLOAT, nodo.childrens_ID[i], nodo.node_ID, MPI_COMM_WORLD);
            } else {
                MPI_Isend(&voltage_obs, 1, MPI_FLOAT, nodo.childrens_ID[i], nodo.node_ID, MPI_COMM_WORLD,&req_v);
            }
            
        }

        // Send measures to ancestor
        if (rank>0){
            if (SYNC){
                MPI_Send(&current_obs, 1, MPI_FLOAT, nodo.ancestor_ID, nodo.node_ID+100, MPI_COMM_WORLD);
                MPI_Send(&active_power_obs, 1, MPI_FLOAT, nodo.ancestor_ID, nodo.node_ID+200, MPI_COMM_WORLD);
                MPI_Send(&reactive_power_obs, 1, MPI_FLOAT, nodo.ancestor_ID, nodo.node_ID+300, MPI_COMM_WORLD);
                sprintf(buffer_logFile, "Node %d is sending current to ancestor %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.ancestor_ID,current_obs, iters);
                nodo.log_file(buffer_logFile);
                sprintf(buffer_logFile, "Node %d is sending active power to ancestor %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.ancestor_ID,active_power_obs, iters);
                nodo.log_file(buffer_logFile); 
                sprintf(buffer_logFile, "Node %d is sending reactive power to ancestor %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.ancestor_ID,reactive_power_obs, iters);
                nodo.log_file(buffer_logFile);             
            } else {
                MPI_Isend(&current_obs, 1, MPI_FLOAT, nodo.ancestor_ID, nodo.node_ID+100, MPI_COMM_WORLD,&req_current[nodo.node_ID]);
                MPI_Isend(&active_power_obs, 1, MPI_FLOAT, nodo.ancestor_ID, nodo.node_ID+200, MPI_COMM_WORLD,&req_power[nodo.node_ID]);
                MPI_Isend(&reactive_power_obs, 1, MPI_FLOAT, nodo.ancestor_ID, nodo.node_ID+300, MPI_COMM_WORLD,&req_reactive[nodo.node_ID]);
                sprintf(buffer_logFile, "Node %d is sending current to ancestor %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.ancestor_ID,current_obs, iters);
                nodo.log_file(buffer_logFile);
                sprintf(buffer_logFile, "Node %d is sending active power to ancestor %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.ancestor_ID,active_power_obs, iters);
                nodo.log_file(buffer_logFile); 
                sprintf(buffer_logFile, "Node %d is sending reactive power to ancestor %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.ancestor_ID,reactive_power_obs, iters);
                nodo.log_file(buffer_logFile);             
            }
        }

        //Receive voltage from ancestor
        float voltage_ancestor_obs;
        if (rank>0){
            int request_complete = 0;

            //cout << "Node " << nodo.node_ID << " is receiving voltage from ancestor " <<  nodo.ancestor_ID  << endl;
            if (SYNC){
                MPI_Recv(&voltage_ancestor_obs, 1, MPI_FLOAT, nodo.ancestor_ID, nodo.ancestor_ID, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                nodo.node_measures.voltage_ancestor = voltage_ancestor_obs;
                nodo.node_observations.voltage_ancestor  = voltage_ancestor_obs;
                //cout << std::fixed << std::setprecision(3) << "Node " << nodo.node_ID << " received voltage from ancestor " << nodo.ancestor_ID << ". Measure: " << voltage_ancestor_obs << endl;
                sprintf(buffer_logFile, "Node %d is receiving voltage from ancestor %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.ancestor_ID,voltage_obs, iters);
                nodo.log_file(buffer_logFile);                
            } else {
                MPI_Irecv(&voltage_ancestor_obs, 1, MPI_FLOAT, nodo.ancestor_ID, nodo.ancestor_ID, MPI_COMM_WORLD,
                        &req_v);
                MPI_Test(&req_v,&request_complete,&status_v);
                if (request_complete){
                    nodo.node_measures.voltage_ancestor = voltage_ancestor_obs;
                    nodo.node_observations.voltage_ancestor  = voltage_ancestor_obs;
                    //cout << std::fixed << std::setprecision(3) << "Node " << nodo.node_ID << " received voltage from ancestor " << nodo.ancestor_ID << ". Measure: " << voltage_ancestor_obs << endl;
                    sprintf(buffer_logFile, "Node %d is receiving voltage from ancestor %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.ancestor_ID,voltage_obs, iters);
                    nodo.log_file(buffer_logFile);                
                }
            }

        }

        //Receive measures from children
        float current_rec_1, current_rec_2;
        float active_power_rec_1, active_power_rec_2;
        float reactive_power_rec_1,reactive_power_rec_2;

        if (n_childs > 0){
            int request_complete_current = 0;
            int request_complete_power = 0;
            int request_complete_reactive = 0;

            if (SYNC){
                MPI_Recv(&current_rec_1, 1, MPI_FLOAT, nodo.childrens_ID[0], nodo.childrens_ID[0]+100, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(&active_power_rec_1, 1, MPI_FLOAT, nodo.childrens_ID[0], nodo.childrens_ID[0]+200, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(&reactive_power_rec_1, 1, MPI_FLOAT, nodo.childrens_ID[0], nodo.childrens_ID[0]+300, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                sprintf(buffer_logFile, "Node %d received current from children %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.childrens_ID[0],current_rec_1, iters);
                nodo.log_file(buffer_logFile);
                sprintf(buffer_logFile, "Node %d received active power from children %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.childrens_ID[0],active_power_rec_1, iters);
                nodo.log_file(buffer_logFile);
                sprintf(buffer_logFile, "Node %d received reactive power from children %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.childrens_ID[0],reactive_power_rec_1, iters);
                nodo.log_file(buffer_logFile);                
                nodo.children_measures[0].current = current_rec_1;
                nodo.children_measures[0].active_power = active_power_rec_1;
                nodo.children_measures[0].reactive_power = reactive_power_rec_1;
            } else {
                MPI_Irecv(&current_rec_1, 1, MPI_FLOAT, nodo.childrens_ID[0], nodo.childrens_ID[0]+100, MPI_COMM_WORLD,
                        &req_current[nodo.childrens_ID[0]]);
                MPI_Irecv(&active_power_rec_1, 1, MPI_FLOAT, nodo.childrens_ID[0], nodo.childrens_ID[0]+200, MPI_COMM_WORLD,
                        &req_power[nodo.childrens_ID[0]]);
                MPI_Irecv(&reactive_power_rec_1, 1, MPI_FLOAT, nodo.childrens_ID[0], nodo.childrens_ID[0]+300, MPI_COMM_WORLD,
                        &req_reactive[nodo.childrens_ID[0]]);
                
                //nodo.children_measures[0].current = .5*(nodo.children_measures[0].current + current_rec_1);
                //nodo.children_measures[0].active_power = .5*(nodo.children_measures[0].active_power + active_power_rec_1);
                //nodo.children_measures[0].reactive_power = .5*(nodo.children_measures[0].reactive_power + reactive_power_rec_1);
                MPI_Test(&req_current[nodo.childrens_ID[0]],&request_complete_current,&status_current[nodo.childrens_ID[0]]);
                MPI_Test(&req_power[nodo.childrens_ID[0]],&request_complete_power,&status_power[nodo.childrens_ID[0]]);
                MPI_Test(&req_reactive[nodo.childrens_ID[0]],&request_complete_reactive,&status_reactive[nodo.childrens_ID[0]]);

                if (request_complete_current){
                    nodo.children_measures[0].current = current_rec_1;
                    //cout << std::fixed << std::setprecision(3) << "Node " << nodo.node_ID << " received current from children " << nodo.childrens_ID[0] << ". Measure: " << current_rec_1 << endl;
                    sprintf(buffer_logFile, "Node %d received current from children %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.childrens_ID[0],current_rec_1, iters);
                    nodo.log_file(buffer_logFile);
                }
                if (request_complete_power){
                    nodo.children_measures[0].active_power = active_power_rec_1;
                    //cout << std::fixed << std::setprecision(3) << "Node " << nodo.node_ID << " received active power from children " << nodo.childrens_ID[0] << ". Measure: " << active_power_rec_1 << endl;
                    sprintf(buffer_logFile, "Node %d received active power from children %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.childrens_ID[0],active_power_rec_1, iters);
                    nodo.log_file(buffer_logFile);                
                }
                if (request_complete_reactive){
                    nodo.children_measures[0].reactive_power = reactive_power_rec_1;
                    //cout << std::fixed << std::setprecision(3) << "Node " << nodo.node_ID << " received reactive power from children " << nodo.childrens_ID[0] << ". Measure: " << reactive_power_rec_1 << endl;
                    sprintf(buffer_logFile, "Node %d received reactive power from children %d. Measure: %f. Iter: %d",nodo.node_ID,nodo.childrens_ID[0],reactive_power_rec_1, iters);
                    nodo.log_file(buffer_logFile);                
                }
            }
        }

        if (n_childs > 1){
            int request_complete_current = 0;
            int request_complete_power = 0;
            int request_complete_reactive = 0;
            if (SYNC){
                MPI_Recv(&current_rec_2, 1, MPI_FLOAT, nodo.childrens_ID[1], nodo.childrens_ID[1]+100, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(&active_power_rec_2, 1, MPI_FLOAT, nodo.childrens_ID[1], nodo.childrens_ID[1]+200, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(&reactive_power_rec_2, 1, MPI_FLOAT, nodo.childrens_ID[1], nodo.childrens_ID[1]+300, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                nodo.children_measures[1].current = current_rec_2;
                nodo.children_measures[1].active_power = active_power_rec_2;
                nodo.children_measures[1].reactive_power = reactive_power_rec_2;
            } else {
                MPI_Irecv(&current_rec_2, 1, MPI_FLOAT, nodo.childrens_ID[1], nodo.childrens_ID[1]+100, MPI_COMM_WORLD,
                        &req_current[nodo.childrens_ID[1]]);
                MPI_Irecv(&active_power_rec_2, 1, MPI_FLOAT, nodo.childrens_ID[1], nodo.childrens_ID[1]+200, MPI_COMM_WORLD,
                        &req_power[nodo.childrens_ID[1]]);
                MPI_Irecv(&reactive_power_rec_2, 1, MPI_FLOAT, nodo.childrens_ID[1], nodo.childrens_ID[1]+300, MPI_COMM_WORLD,
                        &req_reactive[nodo.childrens_ID[1]]);

                //nodo.children_measures[1].current = .5*(nodo.children_measures[1].current + current_rec_2);
                //nodo.children_measures[1].active_power = .5*(nodo.children_measures[1].active_power + active_power_rec_2);
                //nodo.children_measures[1].reactive_power = .5*(nodo.children_measures[1].reactive_power + reactive_power_rec_2);
                MPI_Test(&req_current[nodo.childrens_ID[1]],&request_complete_current,&status_current[nodo.childrens_ID[1]]);
                MPI_Test(&req_power[nodo.childrens_ID[1]],&request_complete_power,&status_power[nodo.childrens_ID[1]]);
                MPI_Test(&req_reactive[nodo.childrens_ID[1]],&request_complete_reactive,&status_reactive[nodo.childrens_ID[1]]);

                
                if (request_complete_current){
                    nodo.children_measures[1].current = current_rec_2;
                    cout << std::fixed << std::setprecision(3) << "Node " << nodo.node_ID << " received current from children " << nodo.childrens_ID[1] << ". Measure: " << current_rec_2 << endl;

                }
                if (request_complete_power){
                    nodo.children_measures[1].active_power = active_power_rec_2;
                    cout << std::fixed << std::setprecision(3) << "Node " << nodo.node_ID << " received active power from children " << nodo.childrens_ID[1] << ". Measure: " << active_power_rec_2 << endl;
                }
                if (request_complete_reactive){
                    nodo.children_measures[1].reactive_power = reactive_power_rec_2;
                    cout << std::fixed << std::setprecision(3) << "Node " << nodo.node_ID << " received reactive power from children " << nodo.childrens_ID[2] << ". Measure: " << reactive_power_rec_2 << endl;

                }
            }
        }
        
        //nodo.state_vector = nodo.generate_state_vector();
        //nodo.write_state_vector();
        nodo.observation_vector = nodo.generate_observation_vector();
        nodo.write_observation_vector();

        iters+=1;

        char buffer [15];
        sprintf(buffer, "Iteracion: %d", iters);
        nodo.log_screen(buffer);
    }

    MPI_Finalize();

    return 0;


}