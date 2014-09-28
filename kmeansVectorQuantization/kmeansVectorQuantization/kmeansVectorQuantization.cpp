// kmeansVectorQuantization.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <string.h>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <iostream>

using namespace std;

void CalculateAutoCorCoeff(long long* s, long long* acc, int N, int p)
{
	
	
	for (int k=0; k<=p ; k++)
	{
			acc[k]=0;
			for(int i=0; i<N-k; i++)
			{
				acc[k] += s[i]*s[i+k];
				
				
			}
			//cout << acc[k] << endl;
			
	}

	
}


void CalculateAlphas(long long* acc, double* alpha, int p) //Implementing Durbin's Algorithm
{
	double* E = new double [p+1];
	double* k = new double [p+1];
	
	double** alphaDurbin = new double* [p+1];
	for(int i=0; i<=p; i++)
	{
		alphaDurbin[p] = new double [p+1];
	}

	E[0] = acc[0];
	

	for(int i=1; i<=p; i++)
	{

		double sum =0;
		for(int j=1; j<i; j++)
		{
			sum += alphaDurbin[j][i-1]*acc[i-j];
		}

		k[i] = (acc[i] - sum)  /E[i-1];

		
		E[i] = (1 - k[i]*k[i])*E[i-1];

		for(int j=1; j<i; j++)
		{
			alphaDurbin[j][i] = alphaDurbin[j][i-1] - k[i]*(alphaDurbin[i-j][i-1]);
		}
	
		alphaDurbin[i][i] = k[i];

	}

	

	for(int i=1; i<=p; i++)
	{

		alpha[i] = alphaDurbin[i][p];
	

	}

	




}






int _tmain(int argc, _TCHAR* argv[])
{
	int p = 12;
	string input = "..\\..\\ceptest\\input.txt";

	ifstream inp;
	inp.open(input.c_str());
	int a;
	ofstream out;
	out.open("output.txt");
	int count=0;
	long long *s = new long long[320];

	while(!inp.eof())
	{
		inp >> s[count];
		//cout << a << endl;
		count++;
	}

	cout << count << endl;

	long long *acc = new long long[13];
	CalculateAutoCorCoeff(s, acc, count, p);
	for(int i=0; i<=p; i++)
	{
		out << acc[i] << endl;
	}
	cin >> a;
	return 0;
}


