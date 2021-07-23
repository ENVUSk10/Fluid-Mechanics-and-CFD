#include<bits/stdc++.h>
using namespace std;
int main(){
float dom_length = 1;
int node_points = 101, iterations = 0;
float node = node_points;
float dx = dom_length/(node - 1);
float U = 0.02, alpha = 0.01;
float pe[5] = {0, 0.5, 1, 2.5, 10};
float T_record[5][node_points], T[node_points], T_new[node_points];
double error_req = 1e-6;


for(int num = 0; num<5 ; num++){
    for( int i = 0; i<node_points; i++){
         T[i] = 0;
         T_new[i] = 0;

    }

    T_new[node_points - 1] = 1;
    T[node_points - 1] = 1;
    float Pe = pe[num];
    float ap = Pe + (2/dx);
    float ae = 1/dx;
    float aw = Pe + (1/dx);

    double error_mag = 1;
    while (error_mag > error_req){
        for( int i = 1; i<node_points - 1; i++)T_new[i]=(aw*T[i-1] + ae*T[i+1])/(ap);
        error_mag = 0;
        for( int i = 0; i<node_points; i++)error_mag = error_mag + abs(T[i] - T_new[i]);
        iterations = iterations + 1;
        for( int i = 0; i<node_points; i++)T[i] = T_new[i];
    }
for(int i = 0; i<node_points; i++)T_record[num][i] = T[i];
}
for( int i = 0; i<node_points; i++)
    cout << T_record[4][i] << " " << i <<endl;
ofstream outfile;
outfile.open("velocity");
if(outfile.is_open()){
    for(int i = 0; i<5; i++){
        for(int j = 0; j<node_points; j++){
            outfile << T_record[i][j] << " " << j << endl;
        }
        outfile << endl;
    }
    outfile.close();
}
else cout << "ERROR";

}
