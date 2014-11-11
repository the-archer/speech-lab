// hmm.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include<iostream>

using namespace std;

int N=5;
int M=32;



void ForwardProcedure(vector<vector<long double>>& a, vector<vector<long double>>& b, vector<long double>& start, vector<int>& obs, int T, vector<vector<long double>>& alpha)
{
	
	vector<long double> init;
	init.push_back(0);
	for(int i=1;i<=N; i++)
	{
		init.push_back(start[i]*b[i][obs[i]]);
	}
	alpha.push_back(init);
	alpha.push_back(init);

	for(int i=2; i<=T; i++)
	{
		vector<long double> temp;
		temp.push_back(0);
		for(int j=1; j<=N; j++)
		{
			
			long double sum=0;
			for(int k=1; k<=N; k++)
			{
				sum+=alpha[i-1][k]*a[k][j];
			}
			temp.push_back(sum*b[j][obs[i]]);
			
		}
		alpha.push_back(temp);

	}

	long double prob=0;
	for(int i=1; i<=N; i++)
	{
		prob+=alpha[T][i];
	}

	return;

}


void BackwardProcedure(vector<vector<long double>>& a, vector<vector<long double>>& b, vector<long double>& start, vector<int>& obs, int T, vector<vector<long double>>& beta)
{
	for(int i=0; i<=T; i++)
	{
		vector<long double> temp(N+1, 0);
		
		beta.push_back(temp);
	}
	for(int i=1; i<=N; i++)
	{
		beta[T][i]=1;
	}
	for(int t=T-1; t>0; t--)
	{
		for(int i=1; i<=N; i++)
		{
			long double sum=0;
			for(int j=1; j<=N; j++)
			{
				sum+=(a[i][j]*b[j][obs[t+1]]*beta[t+1][j]);
			}
			beta[t][i]=sum;

		}

	}

	return;


}


void CheckZeroEntriesInB(vector<vector<long double>>& b)
{
	for(int i=1; i<=N; i++)
	{
		long double max=0;
		int zeroes=0;
		int ind=0;
		for(int j=1; j<=M; j++)
		{
			if(b[i][j]>max)
			{
				max=b[i][j];
				ind =j;
			}
			if(b[i][j]==0)
			{
				zeroes+=1;
			}

		}

		if
	}

}


void ParameterReestimation(vector<vector<long double>>& a, vector<vector<long double>>& b, vector<long double>& start, vector<int>& obs, int T)
{
	vector<vector<long double> > alpha;
	vector<vector<long double> > beta;
	ForwardProcedure(a, b, start, obs, T, alpha);
	BackwardProcedure(a, b, start, obs, T, beta);
	vector<vector<vector<long double> > > prob;
	for(int i=0; i<=T; i++)
	{
		vector< vector<long double> > temp;
		for(int j=0; j<=N; j++)
		{
			vector<long double> v(N+1, 0);
			temp.push_back(v);
		}
		prob.push_back(temp);
	}
	for(int t=1; t<T; t++)
	{
		long double denom=0;
		for(int i=1; i<=N; i++)
		{
			for(int j=1; j<=N; j++)
			{
				denom+=(alpha[t][i]*a[i][j]*b[j][obs[t+1]]*beta[t+1][j]);
			}

		}

		for(int i=1; i<=N; i++)
		{
			for(int j=1; j<=N; j++)
			{
				prob[t][i][j]=(alpha[t][i]*a[i][j]*b[j][obs[t+1]]*beta[t+1][j])/denom;


			}

		}

	}

	vector<vector<long double> > gamma;
	for(int i=0; i<=T; i++)
	{
		vector<long double> temp(N+1, 0);
		gamma.push_back(temp);
	}
	
	for(int t=1; t<=T; t++)
	{
		for(int i=1; i<=N; i++)
		{
			long double sum=0;
			for(int j=1; j<=N; j++)
			{
				sum += prob[t][i][j];
			}
		}

	}

	for(int i=1; i<=N; i++)
	{
		start[i]=gamma[1][i];
	}

	for(int i=1; i<=N; i++)
	{
		for(int j=1; j<=N; j++)
		{
			long double num=0;
			long double den=0;
			for(int t=1; t<T; t++)
			{
				num+=prob[t][i][j];
				den+=gamma[t][i];
			}
			a[i][j]=num/den;
		}
	}
	for(int j=1; j<=N; j++)
	{
		for(int k=1; k<=M; k++)
		{
			long double num=0;
			long double den=0;
			
			for(int t=1; t<=T; t++)
			{
				if(obs[t]==k)
				{
					num+=gamma[t][j];
				}
				den+=gamma[t][j];
			}

			b[j][k]=num/den;


		}
	}

	CheckZeroEntriesinB(b);

	return;
}




long double ViterbiAlgo(vector<vector<long double>>& a, vector<vector<long double>>& b, vector<long double>& start, vector<int>& obs, int T, vector<int>& states)
{
	vector<vector<long double> > score;
	vector<long double> init;
	vector<vector<int> > maxm;
	init.push_back(0);
	vector<int> minit;
	minit.push_back(0);
	for(int i=1; i<=N; i++)
	{
		init.push_back(start[i]*b[i][obs[1]]);
		minit.push_back(0);
	}
	score.push_back(init);
	score.push_back(init);
	maxm.push_back(minit);
	maxm.push_back(minit);

	for(int t=2; t<=T; t++)
	{
		vector<int> temp;
		vector<long double> temp2;
		temp2.push_back(0);

		temp.push_back(0);
		for(int j=1; j<=N; j++)
		{
			long double best=0;
			int bstate=0;
			for(int i=1; i<=N; i++)
			{
				long double cur=score[t-1][i];
				if(cur>best)
				{
					best=cur;
					temp.push_back(i);

					
				}
			}

			temp2.push_back(best*b[j][obs[t]]);
			

		}
		maxm.push_back(temp);
		score.push_back(temp2);

	}
	long double bestscore=0;
	int st=-1;

	for(int i=1; i<=N; i++)
	{
		if(score[T][i]>bestscore)
		{
			bestscore=score[T][i];
			st=i;
		}
	}
	for(int i=0; i<=T; i++)
	{
		states.push_back(0);
	}
	states[T]=st;
	for(int t=T-1; T>0; t--)
	{
		states[t]=maxm[t+1][states[t+1]];
	}

	return bestscore;

}

void Initialise(vector<vector<long double>>& a, vector<vector<long double>>& b, vector<long double>& start)
{
	for(int i=0; i<=N; i++)
	{
		start.push_back(0);
		vector<long double> temp(N+1, 0);
		a.push_back(temp);
		vector<long double> temp2(M+1, 1/long double(M));
		b.push_back(temp2);
	}
	start[1]=1;
	for(int i=1; i<N; i++)
	{
		a[i][i]=0.8;
		a[i][i+1]=0.2;
	}
	a[N][N]=1;

	return;
}





int _tmain(int argc, _TCHAR* argv[])
{
	//LoadObservation(); yet to implement
	


	return 0;
}

