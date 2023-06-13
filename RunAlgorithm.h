#pragma once
#include <iostream>
#include <fstream>
#include <iomanip>
#include <list>
#include <vector>
#include "Gene.h"
using std::cout; using std::cin; using std::endl; using std::list; using std::vector;
class RunAlgorithm
{
private:
	int machine_number;
	int number_of_jobs;
	int population;
	int mutation;
	int sum_of_jobs;
	int opt;
	int print;
	int print_every_x_iterations;
	int function_number;
	int converge;
	int* input_arr;

public:
	RunAlgorithm(const int& option_to_read_input, const int& print, const int& population = 100, const int& mutation = 2, const int& amount_of_different_genes_to_converge=5, const int& function_number=1);
	~RunAlgorithm();
	void startAlgo(const int& amount_of_different_genes);
	void fix_population(vector<Gene*>* arrG);
	void printArray(int arr[]);
	void printSolutions(vector<Gene*>& arrG);
	void setLowerUpperProbabilities(vector<Gene*>* arrG, const double& sum_probabilities, const int& iteration_counter);
	int binarySearchForGene(vector<Gene*>* arrG, const double& probability_randomed);
	void initializationOfPopulation(vector<Gene*>* array_of_Genes, list<Gene*>* temp_list_of_genes);
	int amountOfDifferentGenes(const vector<Gene*>*arrG);
};

