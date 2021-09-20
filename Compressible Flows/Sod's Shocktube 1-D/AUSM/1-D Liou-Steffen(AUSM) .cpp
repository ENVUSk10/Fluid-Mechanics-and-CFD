#include<bits/stdc++.h>
using namespace std;
// Function to find Minimum and Maximum
double max(double x, double y)
{
  if (x>y) return x;
  else return y;
}
double min(double x, double y)
{
  if (x>y) return y;
  else return x;
}

int main()
{
double domain_length = 1.0, nodes = 3000.0;
int node_points = 3000;
double dx = domain_length/(nodes - 1);

double Tinit = 0.0, Tfinal = 0.25, dt = 1e-4;
double FinalTime = (Tfinal - Tinit)/dt;
double gamma = 1.4;

double RhoL, RhoR,cL, cR, etL, etR, RhouL, RhouR, RhoetL, RhoetR, uL, uR, pL, pR, hL, hR;
double MachL, MachR, M1, M2,P1,P2, M, P;

double U[3][node_points+1], F[3][node_points], x[node_points], U_final[3][node_points];

// Left side Variables
RhoL = 1.0;
pL = 1.0;
uL = 0.0;

// Right side variables
RhoR = 0.125;
pR = 0.1;
uR = 0.0;

// Velocity of Sound
cL = pow(((gamma*pL)/RhoL),0.5);
cR = pow(((gamma*pR)/RhoR), 0.5);

// Density*velocity
RhouL = RhoL*uL;
RhouR = RhoR*uR;

etL = pL/(gamma - 1) + RhouL*0.5*uL;
etR = pR/(gamma - 1) + RhouR*0.5*uR;

for(int i = 0; i<node_points+1; i++)
{

    if((dx*i) < 0.5)
    {
		U[0][i] = RhoL;
        U[1][i] = RhouL;
        U[2][i] = etL;

	}

     if((dx*i) >= 0.5)
	{
        U[0][i] = RhoR;
        U[1][i] = RhouR;
        U[2][i] = etR;

	}

}

for(int time = 1; time<= FinalTime; time++)
{
	for(int i = 0; i<node_points; i++)
    {
	   RhoL = U[0][i];
       RhoR = U[0][i+1];

       RhouL = U[1][i];
       RhouR = U[1][i+1];

       RhoetL = U[2][i];
       RhoetR = U[2][i+1];

       uL = U[1][i]/RhoL;
       uR = U[1][i+1]/RhoR;

       pL = (RhoetL - RhouL*uL*0.5)*(gamma - 1);
       pR = (RhoetR - RhouR*uR*0.5)*(gamma - 1);

       hL = (RhoetL + pL)/RhoL;
       hR = (RhoetR + pR)/RhoR;

       cL = pow(((gamma*pL)/RhoL),0.5);
       cR = pow(((gamma*pR)/RhoR), 0.5);

       MachL = uL/cL;
       MachR = uR/cR;

	   if(MachL <= -1.0)
	   {
			M1 = 0.0;
            P1 = 0.0;
       }
        else if(MachL < 1.0)
        {
			M1 = 0.25*(MachL + 1)*(MachL + 1);
            P1 = 0.25*(2 - MachL)*pow((MachL + 1),2)*pL;
		}
        else
        {
			M1 = MachL;
            P1 = pL;
        }


        if(MachR <= -1.0)
        {
			M2 = MachR;
            P2 = pR;
        }
		else if(MachR < 1.0)
        {
			M2 = -0.25*(MachR - 1)*(MachR - 1.0);
            P2 = 0.25*(2.0 + MachR)*pow((MachR - 1.0),2)*pR;
        }
		else
        {
			M2 = 0.0;
            P2 = 0.0;
        }

		M = M1 + M2;
        P = P1 + P2;
		
        F[0][i] = max(0.0, M)*RhoL*cL + min(0.0,M)*RhoR*cR;
        F[1][i] = max(0.0, M)*RhouL*cL + min(0.0,M)*RhouR*cR + P;
        F[2][i] = max(0.0, M)*RhoL*hL*cL + min(0.0,M)*RhoR*hR*cR;

	}

	for(int i = 1; i<node_points; i++)
    {
		U[0][i] = U[0][i] - (dt/dx)*(F[0][i] - F[0][i-1]);
        U[1][i] = U[1][i] - (dt/dx)*(F[1][i] - F[1][i-1]);
        U[2][i] = U[2][i] - (dt/dx)*(F[2][i] - F[2][i-1]);

    }

	U[0][0] = U[0][1];
    U[1][0] = U[1][1];
    U[2][0] = U[2][1];
}
for(int i =0; i<=2; i++)
{
	for(int j = 0; j< node_points; j++)
	{
		U_final[i][j] = U[i][j];
	}
}

for(int i =0; i<node_points; i++)
{
	x[i] = i*dx;
}

// Output Variables Density, Velocity, Pressure
ofstream outfile, outfile1, outfile2;
outfile.open("Density.dat");
if(outfile.is_open()){
    outfile << endl;
	for(int i = 0; i<node_points; i++){

        outfile << x[i] << " " <<  U_final[0][i] ;
        outfile << endl;
    }
    outfile.close();
}
else cout << "ERROR";

outfile.open("Velocity.dat");
if(outfile.is_open()){
    outfile << endl;
	for(int i = 0; i<node_points; i++){

        outfile << x[i] << " " <<  U_final[1][i]/U_final[0][i];
        outfile << endl;
    }
    outfile.close();
}
else cout << "ERROR";

outfile.open("Pressure.dat");
if(outfile.is_open()){
    outfile << endl;
	for(int i = 0; i<node_points; i++){

        outfile << x[i] << " " <<  (gamma-1.) * ( U_final[2][i] - 0.5 * U_final[1][i] * U_final[1][i] / U_final[0][i] );
        outfile << endl;
    }
    outfile.close();
}
else cout << "ERROR";
}
