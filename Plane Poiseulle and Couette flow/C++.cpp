// Same as MATLAB
#include<bits/stdc++.h>
using namespace std;
int main()
{
float domain_size, Node_points, dy,  dpdx = 1, mu = 2, velocity;
int iteration = 0, node_points;

double error_mag = 1, error_req = 1e-6;
double error_record[10000];

cout << "Enter the domain size and the node points" << endl;
cin >> domain_size >> node_points;
Node_points = node_points;
vector<float> u(node_points);
vector<float> u_new(node_points);

dy = domain_size/(Node_points - 1);

string couette;
cout << "Do you want Couette condition (y/n)" << endl;
cin >> couette;


if(couette == "y" ||couette == "Y"){
    cout << "Enter the Velocity" << endl;
    cin >> velocity;
    u[node_points-1] = velocity;
    u_new[node_points-1] = velocity;
}

while (error_mag > error_req)
{

    for(int i = 1; i<u.size()-1; i++)u_new[i] = (u[i-1] + u[i+1] + (dy*dy*dpdx/mu))/2;


    error_mag = 0;
    iteration++;
    for(int i = 1; i<u.size()-1; i++)error_mag=error_mag+abs(u[i] - u_new[i]);
    error_record[iteration] = error_mag;

    if(!(iteration%100))cout << error_record[iteration] <<" " <<iteration << endl;

    for(int i = 1; i<u.size()-1; i++)u[i] = u_new[i];

}
cout << "Solution Convereged" << endl;

// writing the velocity and error in a file
ofstream outfile;
ofstream outfile2;
outfile.open("velocity");
if(outfile.is_open()){
    for(int i = 0; i<u.size(); i++) outfile << i << " " << u[i] << endl;
    outfile.close();
}
outfile2.open("Residual error");
if(outfile2.is_open()){
    for(int i = 0; i<iteration; i++) outfile2 << i <<" " << error_record[i] << endl;
    outfile2.close();
}
else cout << "ERROR";

}