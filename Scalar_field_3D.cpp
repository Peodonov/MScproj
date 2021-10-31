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

class Data {
        int    key;
        double value;
    };

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
      file1.write ((char*)&data_, sizeof(data_));
      file1.close();
      cout << "Lattice saved" << endl;
    }

    void load_lattice()
    {  
      std::vector<T> vec;
      std::ifstream file2;
      file2.open("lattice_data.bin", ios::in | ios::binary);
      file2.read((char*)&vec, sizeof(data_));
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


double action(Field<double> sigma, int i, int j, int k, double m, double a)
{
  return ((6.0/pow(a,2.0) + pow(m,2.0))*pow(sigma(i,j,k),2.0)*0.5*pow(a,3.0) 
      - (sigma(i+1,j,k) + sigma(i-1,j,k) 
      + sigma(i,j + 1,k) + sigma(i,j-1,k) 
      + sigma(i,j,k + 1) + sigma(i,j,k - 1))*(a*sigma(i,j,k)*0.5));
}


template <typename T>
void update(double a, double m, double delta, Field<T>& sigma)
{
  for(int i=0;i<sigma.get_lat_dim()[0];i++)
  {
    for(int j=0;j<sigma.get_lat_dim()[1];j++)
    {
      for(int k=0;k<sigma.get_lat_dim()[2];k++)
      {
        // Compute the current action at lattice site n
        double s_old = action(sigma, i, j, k, m, a);
    
        double sigma_old = sigma(i,j,k);
        
        // Generate new value of x at lattice site n
        sigma(i,j,k) = sigma(i,j,k) + ((double) rand()/RAND_MAX - 0.5)*delta*2.0;

        // Compute difference between old action and new action
        double ds = action(sigma, i, j, k, m, a) - s_old;

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
  int n = 2;
  // The size and dimension of the lattice is given by this vector
  std::vector<int> latsize{n,n,n};
  int size = int(pow(n,3));
  // Initialize the lattice
  Field<double> lattice(latsize, size);





  // set simulation variables
  double a = 1.0;      // lattice spacing
  double m = 0.5;      // mass
  double delta = 3.75;  // configuration space range

  int thermN = 300;    // iterations for thermalization

  for(int l =0; l< thermN; l++)
  {
    update(a, m, delta, lattice);
  }

  lattice.save_lattice();

  Field<double> lattice1(latsize, size);
  lattice1.load_lattice();
  cout << lattice1(1,1,1) << endl;
  
  return 0;
}