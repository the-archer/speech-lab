// kmeansVectorQuantization.cpp : Defines the entry point for the console application.
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

#define PI 3.14159265

const int p = 12;
const int NormMax = 5000;
const int fr_size = 320;
const int fr_shift = 80;
const int cb_size = 32;
const double epsilon = 0.02;
const double delta = 0.000000000001;
using namespace std;


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

void CalculateAutoCorCoeff(double* s, double* acc, int N, int p)
{
	
	
	for (int k=0; k<=p ; k++)
	{
			acc[k]=0;
			for(int i=0; i<N-k; i++)
			{
				acc[k] += s[i]*s[i+k];
				
				
			}
		
			
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


double CalculateEuclDistance(vector<double>& vec1, vector<double>& vec2)
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
	vector<string> fname;
	string dirname = "..\\..\\speech-samples\\classify\\*.txt";
	GetAllFiles(fname, dirname);

	for(int i=0; i<fname.size(); i++)
	{

		string path = "..\\..\\speech-samples\\classify\\";
		cout << fname[i] << endl;
		ifstream inp;
		path += fname[i];
		inp.open(path.c_str());
		//string input = "..\\..\\speech-samples\\raw\\txt\\";
		//string cur = "u";
		
		
		
	
		
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
			map << 0 << endl;

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

void LoadTrainingSet(vector<vector<double>>& cep)
{
	ifstream vec;
	vec.open("..\\..\\speech-samples\\training\\trainingset.txt");

	ifstream map;
	map.open("..\\..\\speech-samples\\training\\trainingmap.txt");

	string vowel;

	int count=0;

	while(!vec.eof())
	{
		map >> vowel;

		vector<double> cur(p+1);
		
		switch(vowel[0])
		{
		case 'a':
			cur[0]=0;
			break;
		case 'e':
			cur[0]=1;
			break;
		case 'i':
			cur[0]=2;
			break;
		case 'o':
			cur[0]=3;
			break;
		case 'u':
			cur[0]=4;
			break;
	

		}
		
		for(int i=1; i<=p; i++)
		{
			vec >> cur[i];

		}

		cep.push_back(cur);

	}

	vec.close();
	map.close();

	return;

}


void ChooseInitialCodeVec(vector<vector<double>>& fvec, vector<vector<double>>& codevec, vector<vector<int>>& cluster)
{
	int tsetsize = fvec.size();
	bool* gen = new bool[tsetsize];
	
	for(int i=0; i<tsetsize; i++)
	{
		gen[i] = false;
		
	}
	srand (time(NULL));

	for(int i=0; i<cb_size; i++)
	{
		int r = rand()%tsetsize;
		if(gen[r])
		{
			i--;
			continue;
		}
		gen[r]=true;

		codevec.push_back(fvec[r]);
		vector<int> emp;
		cluster.push_back(emp);

	}

	return;

}

void AssignToCluster(vector<vector<double>>& fvec, vector<vector<double>>& codevec, vector<vector<int>>& cluster, int m)
{
	int tsetsize = fvec.size();
	for(int i=0; i<cb_size; i++)
	{

		cluster[i].clear();
	}
	fstream log;
	log.open("logfile.txt", fstream::app);
	for(int i=0; i<tsetsize; i++)
	{
		double min_dist = DBL_MAX;
		int low=-1;
		for(int j=0; j<m; j++)
		{
			double dist = CalculateEuclDistance(fvec[i], codevec[j]); 
			if(dist < min_dist)
			{
				min_dist = dist;
				low = j;
			}

		}
		
		cluster[low].push_back(i);

	}
	for(int i=0; i<m; i++)
	{
		log << "Cluster No.: " << i << " No. of vectors assigned: " << cluster[i].size() << endl;
	}

	log.close();

	return;
}

void UpdateCodeVec(vector<vector<double>>& fvec, vector<vector<double>>& codevec, vector<vector<int>>& cluster, int m)
{
	for(int i=0; i<m; i++)
	{
		if(cluster[i].size()==0) //If empty cluster, choose another random codevector
		{
			cout << "Empty Cluster!" << endl;
			srand (time(NULL));
			int r = rand()%(fvec.size());

			codevec[i] = fvec[r];
			continue;

		}
		
		vector<double> centroid(p+1, 0);
	
		int count[5]={0};
		for(int j=0; j<cluster[i].size(); j++)
		{
			count[(int)fvec[cluster[i][j]][0]]++;

			for(int k=1; k<=p; k++)
			{
				centroid[k] += (fvec[cluster[i][j]][k]/cluster[i].size());
			}

			


		}
		codevec[i]  = centroid;
		

		int max_count=0;
		int vow=-1;
		for(int j=0; j<5; j++)
		{
			if(count[j]>max_count)
			{
				max_count=count[j];
				vow = j;
			}
		}
		codevec[i][0]=vow;
	}

	return;
}

double CalculateTotalDistortion(vector<vector<double>>& fvec, vector<vector<double>>& codevec, vector<vector<int>>& cluster, int m)
{
	double dist = 0;
	
	for(int i=0; i<m; i++)
	{
		for(int j=0; j<cluster[i].size(); j++)
		{
			dist += CalculateEuclDistance(fvec[cluster[i][j]], codevec[i]);
			
		}
	}
	
	return dist;
}

void PrintCluster(vector<vector<int>>& cluster)
{
	ofstream cl;
	cl.open("cluster.txt");
	for(int i=0; i<cb_size; i++)
	{
		cl << i << ": " ;
		for(int j=0; j<cluster[i].size(); j++)
		{
			cl << cluster[i][j] << " ";
		}
		cl << endl;
	}
}

void PrintCodeVec(vector<vector<double>>& codevec)
{
	ofstream cv;
	cv.open("codevector.txt");
	for(int i=0; i<cb_size; i++)
	{
		for(int j=0; j<=p; j++)
		{
			cv << codevec[i][j] << " " ;
		}
		cv <<  endl;
	}
}


void RunKMeansAlgo(vector<vector<double>>& fvec, vector<vector<double>>& codevec, int MaxIter)
{
	vector<vector<int>> cluster;
	ChooseInitialCodeVec(fvec, codevec, cluster);
	
	ofstream iter;
	iter.open("itervsdist.txt");
	fstream log;
	log.open("logfile.txt", fstream::app);
	for(int i=1; i<=MaxIter; i++)
	{
		AssignToCluster(fvec, codevec, cluster, cb_size);
		//PrintCluster(cluster);
		UpdateCodeVec(fvec, codevec, cluster, cb_size);
		double dist = CalculateTotalDistortion(fvec, codevec, cluster, cb_size);
		dist = dist/fvec.size();
		iter << i << " " << dist << endl;
		
		log << i << " Distortion: " << dist << endl;
		cout << i << " Distortion: " << dist << endl;
		

	}
	iter.close();
	log.close();
	return;

}


void GetCentroid(vector<vector<double>>& fvec, vector<double>& centroid)
{
	int size = fvec.size();
	for(int i=0; i<size; i++)
	{
		for(int j=1; j<=p; j++)
		{
			centroid[j]+=(fvec[i][j]/size);
		}

	}

	return;

}


void SplitCodebook(vector<vector<double>>& codevec, int m)
{
	for(int i=0; i<m/2; i++)
	{
		vector<double> newcentroid(p+1, 0);
		for(int j=1; j<=p; j++)
		{
			newcentroid[j]=codevec[i][j]*(1+epsilon);
			codevec[i][j]=codevec[i][j]*(1-epsilon);
		}
		codevec.push_back(newcentroid);
	}
}

void RunModifiedKMeansAlgo(vector<vector<double>>& fvec, vector<vector<double>>& codevec, int MaxIter)
{
	vector<vector<int>> cluster(cb_size, 0);
	
	
	vector<double> centroid(p+1, 0);
	GetCentroid(fvec, centroid);
	codevec.push_back(centroid);
	
	double last=DBL_MAX;
	ofstream iter;
	iter.open("cbsizevsdist.txt");
	fstream log;
	log.open("logfile.txt", fstream::app);
	for(int m=2; m<=cb_size; m*=2)
	{
		SplitCodebook(codevec, m);
		for(int i=1; i<=MaxIter; i++)
		{
			



			AssignToCluster(fvec, codevec, cluster, m);
			
			//PrintCluster(cluster);
			UpdateCodeVec(fvec, codevec, cluster, m);
			double dist = CalculateTotalDistortion(fvec, codevec, cluster, m);
			dist = dist/fvec.size();
			if(last-dist<delta)
			{
				break;
			}
			last=dist;
			iter << m << " " << i << " Distortion: " << dist << endl;
			cout << m << " " << i << " Distortion: " << dist << endl;
			iter << m << " " << i << " " << dist << endl;
			
		

		}

	}
	log.close();
	iter.close();
	return;

}

int AssignCodeVector(vector<double>& cep, vector<vector<double>>& codevec)
{
	double dist = DBL_MAX;
	int lab=-1;
	for(int i=0; i<codevec.size(); i++)
	{
		double cur = CalculateEuclDistance(cep, codevec[i]);
		if(cur<dist)
		{
			dist = cur;
			lab = i+1;
		}
	}

	return lab;
}


void GenerateCodeBookUsingKMeans(bool modifiedKMeans)
{
	vector <vector<double>> fvec;
	vector <vector<double>> codevec;

	
	LoadTrainingSet(fvec);
	if(!modifiedKMeans)
	RunKMeansAlgo(fvec, codevec, 50);
	else
	RunModifiedKMeansAlgo(fvec, codevec, 100);
	
	ofstream cb;
	cb.open("codebook.txt");
	for(int i=0; i<cb_size; i++)
	{
		for(int j=0; j<p; j++)
		{
			cb << codevec[i][j] << " " ;
		}
		cb << codevec[i][p] << endl;

	}

	cb.close();
}

void LoadCodeBook(vector<vector<double>>& codevec)
{
	ifstream cb;
	cb.open("codebook.txt");

	while(!cb.eof())
	{
		vector<double> cur(p+1);
		for(int j=0; j<=p; j++)
		{
			cb >> cur[j];
		}
		codevec.push_back(cur);

	}
	return;
}

void ClassifyVectors()
{
	vector<vector<double>> codevec;

	LoadCodeBook(codevec);

	vector<string> fname;
	string dirname = "..\\..\\speech-samples\\classify\\*.txt";
	GetAllFiles(fname, dirname);

	for(int i=0; i<fname.size(); i++)
	{
		string path = "..\\..\\speech-samples\\classify\\";
		cout << fname[i] << endl;
		ifstream inp;
		path += fname[i];
		inp.open(path.c_str());
		vector<int> quant;
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
		int labcount[5]={0};
		int framecount=0;
		for(int k=0; (k+fr_size)<=len; )
		{
			framecount++;
			//out << i << endl;
			double* s = new double[fr_size];
			for(int j=0; j<fr_size; j++)
			{
				s[j] = sample[k+j];
			}
			int N = fr_size;
			double* cep = new double[p+1];
			GetFeatureVector(s, cep, N, p);
			vector<double> cepvec(p+1);
			for(int m=0; m<=p; m++)
			{
				cepvec[m] = cep[m];
			}
			int lab = AssignCodeVector(cepvec, codevec);
			
			//labcount[lab]++;
			quant.push_back(lab);

			if((k+fr_size+fr_shift)>len && (k+fr_size)<len)
			{
				k += (len - k - fr_size);
				
			}
			else
			{
				k += fr_shift;
			}
		
	
		}
		delete(sample);
		//double* sample = new double[160000];
		string path2 = "..\\..\\speech-samples\\classifyout\\";
		ofstream oop;
		path2 += fname[i];
		oop.open(path2.c_str());
		for(int i=0; i<quant.size(); i++)
		{
			oop << quant[i] << endl;
		}
		oop.close();

	}

	return;

}





int AssignLabel  (vector<double>& cep, vector<vector<double>>& codevec)
{
	double dist = DBL_MAX;
	int lab=-1;
	for(int i=0; i<codevec.size(); i++)
	{
		double cur = CalculateEuclDistance(cep, codevec[i]);
		if(cur<dist)
		{
			dist = cur;
			lab = codevec[i][0];
		}
	}

	return lab;
}




void TestCodeBook()
{
	vector<vector<double>> codevec;

	LoadCodeBook(codevec);
	vector<string> fname;
	string dirname = "..\\..\\speech-samples\\test\\*.txt";
	GetAllFiles(fname, dirname);

	int correct=0;
	int total=fname.size();

	for(int i=0; i<fname.size(); i++)
	{
		string path = "..\\..\\speech-samples\\test\\";
		cout << fname[i] << endl;
		ifstream inp;
		path += fname[i];
		inp.open(path.c_str());
		int actual = -1;
		switch(fname[i][0])
		{
		case 'a':
			actual=0;
			break;
		case 'e':
			actual=1;
			break;
		case 'i':
			actual=2;
			break;
		case 'o':
			actual=3;
			break;
		case 'u':
			actual=4;
			break;
		}
		double* sample = new double[160000];
		int len=0;
		while(!inp.eof())
		{
			inp >> sample[len];
			len++;
		}

		len -= 1; //Correcting for an extra empty line at the end of every speech file
		inp.close();
		Normalise(sample, len);
		int labcount[5]={0};
		int framecount=0;
		for(int k=0; (k+fr_size)<=len; )
		{
			framecount++;
			//out << i << endl;
			double* s = new double[fr_size];
			for(int j=0; j<fr_size; j++)
			{
				s[j] = sample[k+j];
			}
			int N = fr_size;
			double* cep = new double[p+1];
			GetFeatureVector(s, cep, N, p);
			vector<double> cepvec(p+1);
			for(int m=0; m<=p; m++)
			{
				cepvec[m] = cep[m];
			}
			int lab = AssignLabel(cepvec, codevec);
			
			labcount[lab]++;

			if((k+fr_size+fr_shift)>len && (k+fr_size)<len)
			{
				k += (len - k - fr_size);
				
			}
			else
			{
				k += fr_shift;
			}


		}

		int max=0;
		int val=-1;

		for(int m=0; m<5; m++)
		{
			if(labcount[m]>max)
			{
				max = labcount[m];
				val = m;
			}
		}
		cout << val << endl;
		cout << "Certainty: " << labcount[val]/double(framecount)*100 << endl;
		if(actual==val)
		{
			correct++;
			cout << "Correct: /" << fname[i][0] << "/" << endl;
		}


	}
	cout << "Accuracy: " << double(correct)/total*100 << endl;

}

int _tmain(int argc, _TCHAR* argv[])
{
	
	int choice=0;
	cout << "1. Build Codebook using K Means \n2. Build Codebook using Modified K-Means\n3. Test Codebook\n 4.Classify feature vectors\n 5.Build Training Set\nChoose (1/2/3/4/5):";
	cin >> choice;
	
	choice-=1;
	if(choice<2)
	GenerateCodeBookUsingKMeans(choice);
	else if(choice==2)
	TestCodeBook();
	else if(choice==3)
	ClassifyVectors();
	else
	BuildTrainingSet();
	
	
	
	
	
	
	
	cout << "Done" << endl;
	
	int a;
	cin >> a;


	return 0;
}


