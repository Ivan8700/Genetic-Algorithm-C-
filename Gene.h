#pragma once
#include <random>
#include <iostream>
#include <list>;
#include <vector>;
using std::cout; using std::cin; using std::endl; using std::list;
class Gene
{
private:
	int length;
	int machines;
	int makespan;
	int* gene;
	int* sum_array;
	double value;
	double probability;
	double lowerProb;
	double upperProb;
public:
	Gene(int n, int m);
	Gene(const Gene& g);
	void randomGene(const int& m, const int& n, const int arr[]);
	~Gene();
	bool validGene(int m);
	void updateGene(const int& n, const int& m, const int arr[]);
	int getMakespan();
	double getLowerProb();
	double getUpperProb();
	void setLowerProb(const double& k);
	void setUpperProb(const double& k);
	int getGeneAtIndex(int i);
	void setGeneAtIndex(const int& i, const int& k);
	int* getGene();
	void printGene(const int& lengthOfGene);
	double getValue();
	double getProbability();
	void setProbabiltyOfGene(double q);
	void shortcut(const int& a, const int arr[], const int& n, const int& sum, const int& m, const int& opt);
	void valueGeneByDeviationFromOPTfunction(const int arr[], const int& length, const int& sum, const int& m, const int& opt);
	void valueGeneByDeviationFromOPTsquaredFunction(const int arr[], const int& n, const int& sum, const int& m, const int& opt);
	void valueGeneByMakespanMinusMinimalFunction(const int arr[], const int& length, const int& sum, const int& m, const int& opt);
	void valueGeneByMakespanMinusOPTplusOne(const int arr[], const int& length, const int& sum, const int& m);
	void changeGene(Gene* g, const int& n, const int arr[], const int& print);
	void decodeAndPrint(int n, int m, const int arr[]);
	void mutate(const int& length, const int& m, const int arr[], const int& print, const int& iteration_counter, const int& print_every_x_iterations);
	void correctGene(const int& n, const int& m, const int arr[]);
};

