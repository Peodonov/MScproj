#include<iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include<cmath>
#include<math.h>
#include <stdlib.h>
#include <ctime>
#include <random>
#include <time.h>
#include <numeric>



using namespace std;

template <class T>
class Field
{
  private:
  // vector of type T which will contain the data at each point in the lattice
    std::vector<T>		data_; 
  // vector which contains dimensions of the lattice
    std::vector<int> 		  lat_;	
  // function which inputs the location in the lattice and outputs the data at that point 
    int index_(int x, int y, int z) const 
    {
        // Periodic boundary conditions implemented using modulo operator
        int elem = (x%lat_[0] + lat_[0])%lat_[0] 
        + ((y%lat_[1] + lat_[1])%lat_[1])*(lat_[0]) 
        + ((z%lat_[2] + lat_[2])%lat_[2])*(lat_[0]*lat_[1]);
        return elem;
    }

  public: 
    Field(std::vector<int> latsize, int a) : data_(a), lat_(latsize)
    {
      // Initialize each element to be a vector filled with zeros
      std::vector<double> start(latsize.size(), 0.0);
      for (int i = 0; i < a; i++)
      {
        data_[i] = 0.0;
      }
    }

    void save_lattice() 
    {
      std::ofstream file1;
      file1.open("lattice_data.bin", ios::out | ios::binary);
      file1.write ((char*)&data_, data_.size()*sizeof(double));
      file1.close();
      cout << "Lattice saved" << endl;
    }

    void load_lattice()
    {  
      std::vector<T> vec;
      std::ifstream file2;
      file2.open("lattice_data.bin", ios::in | ios::binary);
      file2.read((char*)&vec, data_.size()*sizeof(double));
      data_ = vec;
      file2.close();
      cout << "Lattice loaded" << endl;
    }

    // Returns number of sites in the lattice
    int get_size() {return data_.size();}
    // Returns the dimension of the lattice
    int get_dim() {return lat_.size();}
    // returns the size of each side of the lattice in vector form
    std::vector<int> get_lat_dim() {return lat_;}
    T  operator() (int x, int y, int z) const { return data_[index_(x,y,z)]; }
    T& operator() (int x, int y, int z)       { return data_[index_(x,y,z)]; }
}; 

class Action
{
    private:
        double ma_;
    public:
        Action(double ma): ma_(ma){}

        double local_S(Field<double>& s, int i, int j, int k)
        {
            return 0.5*(6.0 + ma_*ma_)*s(i,j,k)*s(i,j,k)
                    - (s(i+1,j,k) + s(i-1,j,k) 
                    + s(i,j+1,k) + s(i,j-1,k) 
                    + s(i,j,k+1) + s(i,j,k-1))*s(i,j,k);
        }
        double global_S(Field<double>& s)
        {
            double g_S = 0.0;

            for(int i = 0; i < s.get_lat_dim()[0]; i++)
            {
                for(int j = 0; j < s.get_lat_dim()[1]; j++)
                {
                    for(int k = 0; k < s.get_lat_dim()[2]; k++)
                    {
                        g_S = g_S + 0.5*((6.0 + ma_*ma_)*s(i,j,k)*s(i,j,k)
                        - (s(i+1,j,k) + s(i-1,j,k) 
                        + s(i,j + 1,k) + s(i,j-1,k) 
                        + s(i,j,k + 1) + s(i,j,k - 1))*s(i,j,k));
                    }
                }
            }
            return g_S;
        }

};


template <typename T>
void update(Field<T>& sigma, Action& S, double delta)
{
  for(int i=0;i<sigma.get_lat_dim()[0];i++)
  {
    for(int j=0;j<sigma.get_lat_dim()[1];j++)
    {
      for(int k=0;k<sigma.get_lat_dim()[2];k++)
      {
        // Compute the current action at lattice site n
        double s_old = S.local_S(sigma,i,j,k);
    
        double sigma_old = sigma(i,j,k);
        
        // Generate new value of x at lattice site n
        sigma(i,j,k) = sigma(i,j,k) + ((double) rand()/RAND_MAX - 0.5)*delta*2.0;

        // Compute difference between old action and new action
        double ds = S.local_S(sigma,i,j,k) - s_old;

        // Compute probability
        double prob = std::min(1.0, exp(- double(ds)));
        // Test probability
        if (prob < (double) rand()/RAND_MAX)
        {
          sigma(i,j,k) = sigma_old;
        }
      }
    }
  }
}

double twopoint(Field<double>& sigma, int t)
{
  // initialize the values which will store the field summed over spatial coordinates
  double phi_t = 0;

  // For loops summing over position coordinates
  for(int x = 0; x < sigma.get_lat_dim()[0]; x++)
  {
    for(int y = 0; y < sigma.get_lat_dim()[1]; y++)
    {
      phi_t = phi_t + sigma(x,y,t);
    }
  }
  return phi_t*sigma(0,0,0);
}

double onepoint(Field<double>& sigma, int t)
{
  double phit = 0;
  for(int x = 0; x < sigma.get_lat_dim()[0]; x++)
  {
    for(int y = 0; y < sigma.get_lat_dim()[1]; y++)
    {
      phit = phit + sigma(x,y,t);
    }
  }
  return phit;
}


int main(){

    srand(time(NULL));
    int n = 12;
    // The size and dimension of the lattice is given by this vector
    std::vector<int> latsize{n,n,n};
    int size = n*n*n;
    // Initialize the lattice
    Field<double> lattice(latsize, size);

    // set simulation variables
    double m = 0.1;      // mass
    double delta = 2.0;  // configuration space range

    int thermN = 1000;    // iterations for thermalization
    int N = 1000;         // iterations for average
    Action S(m);


    
    for(int l =0; l< 100; l++)
    {
        update(lattice, S, delta);
    }

    double g1 = S.global_S(lattice);
    double l1 = S.local_S(lattice, 1,1,1);
    lattice(1,1,1) = 1.0;
    cout << S.global_S(lattice) - g1 << endl;
    cout << S.local_S(lattice, 1,1,1) - l1 << endl;


    // stores 2 point function at different times
    std::vector<double> correlation2(n,0.0);

    // stores 1 point function at different times
    std::vector<double> correlation1(n,0.0);

    ofstream file1 ("onept.txt"); //store data in txt file
    ofstream file2 ("twopt.txt");

    for(int l =0; l< thermN; l++)
    {
      update(lattice, S, delta);
      for(int t = 0; t < lattice.get_lat_dim()[2]; t++)
      {
        correlation2[t] = correlation2[t] + twopoint(lattice,t);
        correlation1[t] = correlation1[t] + onepoint(lattice,t);
      }
      file1 << l << ", " << correlation1[0]/double(l + 1)<< "\n";
      file2 << l << ", " << correlation2[0]/double(l + 1)<< "\n";
    }
    file1.close();
    file2.close();

    
    // Stores correlation functions at different times
    std::vector<double> cor1(12,0.0);
    std::vector<double> cor2(12,0.0);
    double onept = 0.0;
    for(int p = 0; p < N; p++)
    {
      for(int i = 0; i<5;i++)
      {
        update(lattice, S, delta);
      }

      for(int t = 0; t < lattice.get_lat_dim()[2]; t++)
      {        
        cor1[t] = cor1[t] + onepoint(lattice,t);
        cor2[t] = cor2[t] + twopoint(lattice,t);
      }
    }

    ofstream file3 ("cor1.txt");
    ofstream file4 ("cor2.txt");

    for(int t = 0; t < n; t++)
    {
        file3 << double(t) << ", " << cor1[t]/double(N) << "\n";
        file4 << double(t) << ", " << cor2[t]/double(N) << "\n";
    }
    file3.close();
    file4.close();
    
    
    return 0;
}