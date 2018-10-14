#include <iostream>
#include<fstream>
#include <time.h>
#include<vector>
#include<algorithm>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>

using namespace std;

#define MAX 107
#define INF 1000000007
#define EPS (1e-12)

ofstream fout("1005081_ratio.csv");
ofstream fout2("1005081_DP.csv");
ofstream fout3("1005081_AP.csv");

void Pivot( long m,long n,double A[MAX+7][MAX+7],long *B,long *N,long r,long c )
{
    long i,j;
    swap( N[c],B[r] );
    A[r][c] = 1/A[r][c];
    for( j=0;j<=n;j++ )
        if( j!=c )
            A[r][j] *= A[r][c];
    for( i=0;i<=m;i++ )
    {
        if( i!=r )
        {
            for( j=0;j<=n;j++ )
                if( j!=c )
                    A[i][j] -= A[i][c]*A[r][j];
            A[i][c] = -A[i][c]*A[r][c];
        }
    }
}

long Feasible( long m,long n,double A[MAX+7][MAX+7],long *B,long *N )
{
    long r,c,i;
    double p,v;
    while( 1 )
    {
        for( p=INF,i=0;i<m;i++ )
            if( A[i][n]<p )
                p = A[r=i][n];
        if( p > -EPS ) return 1;
        for( p=0,i=0;i<n;i++ )
            if( A[r][i]<p )
                p = A[r][c=i];
        if( p > -EPS ) return 0;
        p = A[r][n]/A[r][c];
        for( i=r+1;i<m;i++ )
        {
            if( A[i][c] > EPS )
            {
                v = A[i][n]/A[i][c];
                if( v<p )
                {
                    r=i,p=v;
                }
            }
        }
        Pivot( m,n,A,B,N,r,c );
    }
}

long Simplex( long m,long n,double A[MAX+7][MAX+7],double *b,double &Ret )
{
   long B[MAX+7],N[MAX+7],r,c,i;
   double p,v;
   for( i=0;i<n;i++ ) N[i] = i;
   for( i=0;i<m;i++ ) B[i] = n+i;
   if( !Feasible( m,n,A,B,N ) ) return 0;
   while( 1 )
   {
      for( p=0,i=0;i<n;i++ ) if( A[m][i] > p ) p = A[m][c=i];
      if( p<EPS )
      {
           for( i=0;i<n;i++ ) if( N[i]<n ) b[N[i]] = 0;
           for( i=0;i<m;i++ ) if( B[i]<n ) b[B[i]] = A[i][n];
           Ret = -A[m][n];
           return 1;
       }
       for( p=INF,i=0;i<m;i++ )
       {
           if( A[i][c] > EPS )
           {
               v = A[i][n]/A[i][c];
               if( v<p ) p = v,r = i;
           }
       }
       if( p==INF ) return -1;
       Pivot( m,n,A,B,N,r,c );
    }

}

double  AP(int N,int M,int weight[],int* k,int E[][100])
{
    double A[MAX+7][MAX+7];
    double b[M];
    double Ret;
    int i,j;

    for(i=0;i<N+1;i++)
        for(j=0;j<M+1;j++)
            A[i][j]=0;
    for(i=0;i<M;i++) A[N][i]=-weight[i];
    for(i=0;i<N;i++) A[i][M]=-1;
    for(i=0;i<M;i++)
    {
        for(j=0;j<k[i];j++)
        {
            A[E[i][j]][i]=-1;
        }
    }
    cout<<"LP Matrix : "<<endl;
    for(i=0;i<N+1;i++)
    {
     for(j=0;j<M+1;j++)
        {
            cout<<A[i][j]<<"  ";
        }cout<<endl;
    }
    Simplex(N,M,A,b,Ret);
    int f[N];
    for(int i=0;i<N;i++)
    {
        f[i]=0;
    }
    for(int i=0;i<M;i++)
    {
        for(int j=0;j<k[i];j++)
        {
            f[E[i][j]]++;
        }
    }
    double maximum =0.0;
    for(int i=0;i<N;i++)
    {
        if(f[i] > maximum) maximum =f[i];
    }
    double cost=0.0;
    for(int i=0;i<M;i++)
    {
        if(b[i]>=1/maximum )
        {
            b[i]=1;
            cost = cost+weight[i];
        }
        else b[i]=0;
    }
    cout<<"IDs of Subsets : ";
    for(int i=0;i<M;i++) if(b[i]==1) cout<<i<<"  ";
    cout<<endl;
    return cost;
}

double task2(int N,int M,int W[],int* K,int E[][100])
{
    clock_t start, stop;
    cout<<"Approximation LP approach \n";
    start = clock();
    double approx=AP(N,M,W,K,E);
    cout<<"cost for AP: "<<approx;
    stop = clock();
    printf("\nfor %d elements and %d subsets execution time is %lf seconds\n", N,M,(((double)stop)-start)/CLOCKS_PER_SEC);
    fout3 << N << ", " << (((double)stop)-start)/CLOCKS_PER_SEC << endl;
    cout<<endl;
    return approx;
}


int mask(int nowConsiderIndex,int K[],int E[][100])
{
    int m= (0 << nowConsiderIndex) ;
    for(int i=0;i<K[nowConsiderIndex];i++)
    {
        m= m | (1<<E[nowConsiderIndex][i]);
    }
    return m;
}


int  DP(int coveredMask,int nowConsiderIndex,int N,int M, int weight[], vector<vector<int> > &f,int k[],int E[][100],int printSolution[][100])
{
    f.resize(coveredMask+1);
    f[coveredMask].resize(nowConsiderIndex+1);
    if(nowConsiderIndex==M && coveredMask!=((1<<N)-1))
    {
        f[coveredMask][nowConsiderIndex]=999;
    }
    else if(nowConsiderIndex==M && coveredMask==((1<<N)-1))
    {
        f[coveredMask][nowConsiderIndex]=0;
    }
    else
    {
        int valChoose=DP(coveredMask | mask(nowConsiderIndex,k,E), nowConsiderIndex+1,N,M,weight,f,k,E,printSolution) + weight[nowConsiderIndex];
        int valNotChoose= DP(coveredMask, nowConsiderIndex+1,N,M,weight,f,k,E,printSolution);
        f[coveredMask][nowConsiderIndex] = min ( valChoose, valNotChoose );
        if(valChoose<valNotChoose)
        {
           printSolution[coveredMask][nowConsiderIndex]=1;
        }
        else
        {
            printSolution[coveredMask][nowConsiderIndex]=0;
        }

    }
    return f[coveredMask][nowConsiderIndex];
}

double task1(int N,int M,int W[],int* K,int E[][100])
{
    clock_t start, stop;
    cout<<"Bitmask DP approach \n";
    start = clock();
    int msk = (0 << M);
    vector<vector<int> > F;
    int printSolution[(1<<N) - 1][100];
    int cost=DP(msk,0,N,M,W,F,K,E,printSolution);
    cout<<"IDs of Subsets : ";
    for(int x=0;x<((1<<N) - 1);x++)
    {
        for(int y=0;y<M;y++)
        {
            if(printSolution[x][y]==1)
            {
                cout<<y<<"  ";
                x=x | mask(y,K,E);
            }
        }
    }
    cout<<endl;
    cout<<"cost for DP: "<<cost;
    stop = clock();
    printf("\nfor %d elements and %d subsets execution time is %lf seconds\n", N,M,(((double)stop)-start)/CLOCKS_PER_SEC*1000);
    fout2 << N << ", " << (((double)stop)-start)/CLOCKS_PER_SEC << endl;
    cout<<endl;
    return cost;
}

int main()
{
    ifstream fin ("1005081_input.txt");
    fout2 << "Element, Execution Time" << endl;
    fout3 << "Element, Execution Time" << endl;
    fout << "DP, AP, Ratio" << endl;
    int caseno =1;
    int T,N,M;
    fin >> T;
    for(int i=0;i<T;i++)
    {
        cout<<"\n---------------------------CASE "<<caseno<<" ----------------------------\n\n";
        caseno++;
        fin >> N;
        fin >> M;
        int W[M],K[M],E[M][100];
        for(int j=0;j<M;j++)
        {
            fin>>W[j];
            fin>>K[j];
            for(int k=0;k<K[j];k++)
            {
                fin>>E[j][k];
            }
        }
        double cost=task1(N,M,W,K,E);
        double approx=task2(N,M,W,K,E);
        fout <<cost<<",";
        fout <<approx<<",";
        fout <<approx/cost<<",";
        fout<<endl;
        cout<<endl;
    }
    fin.close();
    return 0;
}
