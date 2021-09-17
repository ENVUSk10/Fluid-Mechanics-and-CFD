#ifndef UTILITY_FUNCS_H_INCLUDED
#define UTILITY_FUNCS_H_INCLUDED
#include<iostream>
// Equating 2-D arrays
int (*Equate( int a[][grid_points], int b[][grid_points]) ) [grid_points]
{
    for (int i = 0; i<grid_points;i++)
    {
        for (int j = 0; j<grid_points;j++)
        {
            b[i][j] = a[i][j];
        }
    }
}

// Print 2-D arrays
void Display( int a[][grid_points])
{
    for(int i = 0; i<grid_points ; i++)
    {
        for(int j = 0; j<grid_points ; j++)
    {
        std :: cout << a[i][j] << " ";
    }
        std :: cout << std :: endl;
    }

}

int (*Dirichlet( int a[][grid_points], char *boundary, int value))[grid_points]
{
    if(boundary == "Top")
    {
         for(int i = 0; i<grid_points ; i++)
            {
                a[0][i] = value;
            }
    }

    else if(boundary == "Bottom")
    {
         for(int i = 0; i<grid_points ; i++)
            {
                a[grid_points-1][i] = value;
            }
    }

    else if(boundary == "Left")
    {
         for(int i = 0; i<grid_points ; i++)
            {
               a[i][0] = value;
            }
    }

    else
    {
         for(int i = 0; i<grid_points ; i++)
            {
               a[i][grid_points-1] = value;
            }
    }

}


int (*Initialize( int a[][grid_points], int value))[grid_points]
{
    for(int i = 0; i<grid_points ; i++)
    {
         for(int j = 0; j<grid_points ; j++)
    {

            a[i][j] = value;
    }


    }

}
#endif // UTILITY_FUNCS_H_INCLUDED
