// hmm.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <string.h>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <math.h>   
#include <vector>
#include <time.h>
#include <float.h>
#include <windows.h>

using namespace std;

int N=5;
int M=32;


void GetAllFiles(vector<string>& fname, string dirname)
{
	HANDLE hFind;
	WIN32_FIND_DATAA data;

	hFind = FindFirstFileA(dirname.c_str(), &data);
	if (hFind != INVALID_HANDLE_VALUE) {
	  do {
		//printf("%s\n", data.cFileName);
		fname.push_back(data.cFileName);
	  } while (FindNextFileA(hFind, &data));
	  FindClose(hFind);
	}

}


long double ForwardProcedure(vector<vector<long double>>& a, vector<vector<long double>>& b, vector<long double>& start, vector<int>& obs, int T, vector<vector<long double>>& alpha)
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

	return prob;

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
		long double a = 1E-20;
		
		if(a*zeroes>max)
		{
			cout << "Problem in checking zeroes" << endl;
		}
		else
		{
			for(int j=1; j<=M; j++)
			{
				if(b[i][j]==0)
				{
					b[i][j]+=a;
				}
			}
			b[i][ind]-=a*zeroes;
		}

	
	}

	return;

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
			gamma[t][i]=sum;
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

	CheckZeroEntriesInB(b);

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
	for(int t=T-1; t>0; t--)
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

void LoadModel(string init, vector<vector<long double>>& a, vector<vector<long double>>& b, vector<long double>& start)
{
	ifstream inp;
	inp.open(init.c_str());
	
	
		long double temp;
		start.push_back(0);
		for(int i=1; i<=N; i++)
		{	
			inp >> temp;
			start.push_back(temp);
		}
		vector<long double> tp(N+1);
		a.push_back(tp);
		for(int i=1; i<=N; i++)
		{
			vector<long double> tp;
			tp.push_back(0);
			for(int j=1; j<=N; j++)
			{
				inp >> temp;
				tp.push_back(temp);
			}
			a.push_back(tp);
		}
		vector<long double> tp2(M+1);
		b.push_back(tp2);
		for(int i=1; i<=N; i++)
		{
			vector<long double> tp;
			tp.push_back(0);
			for(int j=1; j<=M; j++)
			{
				inp >> temp;
				tp.push_back(temp);
			}
			b.push_back(tp);
		}

	
	return;
}

void CopyVector(vector<long double>& x, vector<long double>& y)
{
	for(int i=0; i<x.size(); i++)
	{
		y.push_back(x[i]);
	}
	return;
}

void CopyVector2(vector<vector<long double>>& x, vector<vector<long double>>& y)
{
	for(int i=0; i<x.size(); i++)
	{
		y.push_back(x[i]);
	}
	return;
}

void AddModel(vector<vector<long double>>& a1, vector<vector<long double>>& b1, vector<long double>& start1, vector<vector<long double>>& a, vector<vector<long double>>& b, vector<long double>& start)
{
	for(int i=0; i<=N; i++)
	{
		for(int j=0; j<=N; j++)
		{
			a1[i][j] += a[i][j];
		}
		for(int j=0; j<=M; j++)
		{
			b1[i][j] += b[i][j];
		}
		start1[i] += start[i];
	}
	return;

}

void DivideModel(vector<vector<long double>>& a, vector<vector<long double>>& b, vector<long double>& start, int k)
{
	for(int i=0; i<=N; i++)
	{
		for(int j=0; j<=N; j++)
		{
			a[i][j] /= k;
		}
		for(int j=0; j<=M; j++)
		{
			b[i][j] /= k;
		}
		start[i] /= k ;
	}
}

void Train()
{
	string init="init.txt";
	//cout << "File name with initial config:" ;
	//cin >> init;
	int k=10;
	//cout << "Enter number of different variatons of a word to train:";
	//cin >> k;
	string models="zero.txt";
	//cout << "Enter file name of file with lists of inputs:";
	//cin >> models;

	vector<vector<long double> > a_orig;
	vector<vector<long double> > b_orig;
	vector<long double> start_orig;
	LoadModel(init, a_orig, b_orig, start_orig);
	vector<vector<long double> > a_new;
	vector<vector<long double> > b_new;
	vector<long double> start_new;
	for(int i=0; i<=N; i++)
	{
		vector<long double> temp(N+1, 0);
		vector<long double> temp2(M+1, 0);
		a_new.push_back(temp);
		b_new.push_back(temp2);
		start_new.push_back(0);

	}
	
	
	ifstream inp;
	inp.open(models.c_str());
	string f;
	for(int i=0; i<k; i++)
	{
		inp >> f;
		cout << f << endl;
		ifstream ob;
		ob.open(f.c_str());
		vector <int> obs;
		int temp;
		int T=0;
		obs.push_back(0);
		while(!ob.eof())
		{
			ob >> temp;
			obs.push_back(temp);
			T+=1;
			
		}
		ob.close();
		vector<vector<long double> > a;
		vector<vector<long double> > b;
		vector<long double> start;
		CopyVector2(a_orig, a);
		CopyVector2(b_orig, b);

		CopyVector(start_orig, start);
		vector <vector<long double> > alpha;
		//cout << "here" << endl;
		long double old_prob=0;
		long double prob = ForwardProcedure(a, b, start, obs, T, alpha);
		cout << prob << endl;
		int count=0;
		while(prob-old_prob>(1E-150) || count==0)
		{
			cout << count << endl;
			count++;
			old_prob = prob;
			vector<int> states;
			ViterbiAlgo(a_orig, b_orig, start_orig, obs, T, states);
			ParameterReestimation(a, b, start, obs, T);
			alpha.clear();
			prob = ForwardProcedure(a, b, start, obs, T, alpha);
			
			cout << prob << endl;
			if(count>100)
				break;
			
		}
		
		/*for(int i=1; i<=T; i++)
		{
			cout << states[i] << " ";
		}
		cout << endl;*/
		
		
		
		
	
		AddModel(a_new, b_new, start_new, a, b, start);


	}
	inp.close();

	DivideModel(a_new, b_new, start_new, k);

	ofstream nmod;
	nmod.open("nmodel.txt");
	for(int i=1; i<=N; i++)
	{
		nmod << start_new[i] << " " ;
	}
	nmod << endl;

	for(int i=1; i<=N; i++)
	{
		for(int j=1; j<=N; j++)
		{

			nmod << a_new[i][j] << " " ;
		}

		nmod << endl;
	}

	for(int i=1; i<=N; i++)
	{
		for(int j=1; j<=M; j++)
		{
			nmod << b_new[i][j] << " " ;

		}
		nmod << endl;
	}
	nmod.close();

}

string IdentifyWord(vector<int>& obs, int T)
{
	string dirname="models\\*.txt";
	vector<string> fname;
	GetAllFiles(fname, dirname);
	long double best=0;
	int ind=-1;
	for(int i=0; i<fname.size(); i++)
	{
		//ifstream inp;
		string path="models\\" + fname[i];
		vector<vector<long double>> a;
		vector<vector<long double>> b;
		vector<long double> start;
		LoadModel(path, a, b, start);


		vector <vector<long double> > alpha;
		long double prob = ForwardProcedure(a, b, start, obs, T, alpha);
		if(prob>best)
		{
			best=prob;
			ind = i;
		}


	}
	return fname[ind].substr(0, fname[ind].size()-4);

}



void Test()
{
	string dirname="test\\*.txt";
	vector<string> fname;
	GetAllFiles(fname, dirname);
	int total=fname.size();
	int correct=0;
	for(int i=0; i<fname.size(); i++)
	{
		cout << fname[i] << endl;
		ifstream inp;
		string path="test\\" + fname[i];
		inp.open(path.c_str());
		int temp;
		vector <int> obs;
		obs.push_back(0);
		int T=0;
		while(!inp.eof())
		{
			inp >> temp;
			if(temp>32)
				temp=32;
			obs.push_back(temp);
			T+=1;

		}
		string word = IdentifyWord(obs, T);
		cout << word << endl;
		size_t found=fname[i].find(word);
		if (found!=std::string::npos)
		{
			cout << "Correct" << endl;
			correct+=1;
		}
		else
		{
			cout << "Incorrect" << endl;
		}

	}
	cout << "Accuracy: " << double(correct)/total*100 << endl;
}


int _tmain(int argc, _TCHAR* argv[])
{
	//Train();
	Test();
	string init;
	cout << "Done" << endl;
	cin >> init;

	return 0;
}

