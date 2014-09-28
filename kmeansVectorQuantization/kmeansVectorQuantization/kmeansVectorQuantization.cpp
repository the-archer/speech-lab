// kmeansVectorQuantization.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <string.h>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <math.h>       

#define PI 3.14159265

const int p = 12;
const int NormMax = 5000;
const int fr_size = 320;
const int fr_shift = 80;

using namespace std;

void CalculateAutoCorCoeff(double* s, double* acc, int N, int p)
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


void CalculateAlphas(double* acc, double* alpha, int p) //Implementing Durbin's Algorithm
{
	double* E = new double [p+1];
	double* k = new double [p+1];
	
	double** alphaDurbin = new double* [p+1];
	for(int i=0; i<=p; i++)
	{
		alphaDurbin[i] = new double [p+1];
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


void CalculateCepCoeff(double* alpha, double* cep, int p, int q)
{
	for(int i=1; i<=q; i++)
	{
		double cur=0;
		for(int k=1; k<i; k++)
		{
			cur += (k/double(i))*cep[k]*alpha[i-k];
		}
		if(i<=p)
		{
			cur += alpha[i];
		}
		
		cep[i] = cur;

	}


}

void ApplyHammingWindow(double* s, int N)
{
	for(int i=0; i<N; i++)
	{
		s[i] *= (0.54 - double(0.46)*cos((2*PI*i)/(N-1)));
		
	}


}

void ApplyParameterWeighting(double* cep, int p)
{
	for(int i=1; i<=p; i++)
	{
		cep[i] *= (1 + (p/2)*sin((PI*i)/p));
	}

}


void Normalise(double* sample, int len)
{
	
    double maximum = 0;

    for(int j=0; j<len; j++)
    {

        if(sample[j]>0 && sample[j]>maximum)
        {
            maximum = sample[j];

        }
        else if(sample[j]<0 && (-1*sample[j])>maximum)
        {
            maximum = -1*sample[j];

        }




    }

    double adjust = NormMax/maximum;

    for(int i=0; i<len; i++)
    {
        sample[i] *= adjust;


    }
    return;


}

void GetFeatureVector(double* s, double* cep, int N, int p)
{

	ApplyHammingWindow(s, N);
	double *acc = new double[p+1];
	CalculateAutoCorCoeff(s, acc, N, p);

	double* alpha = new double[p+1];

	CalculateAlphas(acc, alpha, p);

	CalculateCepCoeff(alpha, cep, p, p);
	ApplyParameterWeighting(cep, p);



}


double CalculateDistance(double* vec1, double* vec2, int p)
{
	double dist=0;
	for(int i=1; i<=p; i++)
	{
		dist += pow((vec1[i] - vec2[i]), 2);
	}
	return dist;
}




void BuildTrainingSet()
{
	ofstream out;
	out.open("output.txt");
	
	fstream vec;
	vec.open("..\\..\\speech-samples\\training\\trainingset.txt", fstream::app);

	fstream map;
	map.open("..\\..\\speech-samples\\training\\trainingmap.txt", fstream::app);


	for(int i=1; i<=20; i++)
	{
		string input = "..\\..\\speech-samples\\raw\\txt\\";
		string cur = "u";
		if(i<10)
		{
			cur += "0";
		}
		cur +=  to_string((_ULonglong)i);
		out << cur << endl;
		input += cur + ".txt";
		ifstream inp;
		inp.open(input.c_str());
		int len=0;
		double* sample = new double[160000];
		while(!inp.eof())
		{
			inp >> sample[len];
			len++;
		
		}
		len -= 1; //Correcting for an extra empty line at the end of every speech file
		inp.close();
		Normalise(sample, len);

		for(int i=0; (i+fr_size)<=len; )
		{
			out << i << endl;
			double* s = new double[fr_size];
			for(int j=0; j<fr_size; j++)
			{
				s[j] = sample[i+j];
			}
			int N = fr_size;
			double* cep = new double[p+1];
			GetFeatureVector(s, cep, N, p);
			for(int j=1; j<p; j++)
			{
				vec << cep[j] << " ";
			}
			vec << cep[p] << endl;
			map << cur << endl;

			if((i+fr_size+fr_shift)>len && (i+fr_size)<len)
			{
				i += (len - i - fr_size);
				
			}
			else
			{
				i += fr_shift;
			}


		}


	}


	

	

	
	out.close();
	vec.close();
	return;



}




int _tmain(int argc, _TCHAR* argv[])
{
		
	

		


	return 0;
}


