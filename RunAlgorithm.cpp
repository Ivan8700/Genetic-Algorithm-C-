#include "RunAlgorithm.h"
RunAlgorithm::RunAlgorithm(const int& option_to_read_input, const int& print, const int& population, const int& mutation, const int& amount_of_different_genes_to_converge, const int& function_number) : population(population), mutation(mutation), print(print), function_number(function_number)
{
	converge = amount_of_different_genes_to_converge;
	sum_of_jobs = 0;
	number_of_jobs = 0;
	std::setprecision(30);
	int cntODD = 0;
	if (option_to_read_input == 0) //Read input from .txt file
	{
		std::ifstream infile;
		infile.open("inputs/input1.txt");
		infile >> machine_number;
		cout << "Enter number of machines : ";
		cout << machine_number << "\n"; 
		infile >> number_of_jobs;
		cout << "Enter Amount of numbers : ";
		cout << number_of_jobs << "\n"; 
		input_arr = new int[number_of_jobs];
		if (infile)
		{
			int i = 0;
			while (infile >> input_arr[i] && i < number_of_jobs)
			{
				sum_of_jobs += input_arr[i];
				if (input_arr[i] % 2 == 1)
					cntODD++;
				i++;
			}

			infile.close();
		}
		else
		{
			cout << "File was not found" << endl;
		}
	}
	else
	{
		if (option_to_read_input == 1) //Manual insertion input
		{
			cout << "Enter amount of machines : "; cin >> machine_number;
			cout << "/n Enter amount of jobs : "; cin >> number_of_jobs; cout << "/n";
			input_arr = new int[number_of_jobs];
			for (int i = 0; i < number_of_jobs; i++)
			{
				cout << "Enter a job : ";
				cin >> input_arr[i]; cout << "/n";
				this->sum_of_jobs += input_arr[i];
				if (input_arr[i] % 2 == 1)
					cntODD++;
			}
		}
		else //Random input
		{
			srand((unsigned int)time(NULL));
			cout << "Enter amount of machines : "; cin >> machine_number;
			cout << "/n Enter amount of jobs : "; cin >> number_of_jobs; cout << "/n";
			input_arr = new int[number_of_jobs];
			for (int i = 0; i < number_of_jobs; i++)
			{
				input_arr[i] = (rand() % 40) + 21; // Size of jobs [21,60]
				this->sum_of_jobs += input_arr[i];
				if (input_arr[i] % 2 == 1)
					cntODD++;
			}
		}
	}
	opt = ((double)(sum_of_jobs) / machine_number);
	opt = (opt % 2 == 1) ? (opt + 1) : opt;
	cout << "==================================================\n";
	cout << "The input is : \n";
	printArray(input_arr);
	if (print >0)
	{
		cout << "\nThe sum is " << sum_of_jobs << "\n";
		cout << "OPT is bound by " << opt << " \n\n";
	}
	cout << "==================================================\n";
	startAlgo(amount_of_different_genes_to_converge);
}
RunAlgorithm::~RunAlgorithm()
{
	delete[] input_arr;
}
/*Randomization of genes for the initial population*/
void RunAlgorithm::initializationOfPopulation(vector<Gene*>* array_of_Genes, list<Gene*>* temp_list_of_genes)
{
	for (int i = 0; i <population; i++) //create 'population' genes one by one 
	{
		Gene* G = new Gene(number_of_jobs,machine_number);
		G->randomGene(number_of_jobs,machine_number, input_arr); //call a function to random the encoding of the gene.
		array_of_Genes->insert(array_of_Genes->begin(), G);
	}
}
/*Check how many different genes in the generation (used for discoverying how much did we converge) */
int RunAlgorithm::amountOfDifferentGenes(const vector<Gene*>* arrG)
{
	list<Gene*>* l = new list<Gene*>;
	Gene* g1G = new Gene(*arrG->at(0));
	l->push_front(g1G);
	int cnt = 1;
	bool found, flag;
	for (int i = 1; i < population; i++) //Traverse the popultion
	{
		flag = false;
		for (list<Gene*>::iterator it = l->begin(); it != l->end(); ++it)
		{
			found = true;
			for (int j = 0; j < number_of_jobs; j++)
			{
				if (arrG->at(i)->getGeneAtIndex(j) != (**it).getGeneAtIndex(j))
				{
					found = false;
					break;
				}
			}
			if (found) //if found = true, the gene is alreay there, don't add it and go to the element.
			{
				flag = true;
				break;
			}
		}
		if (!flag) //flag= false means the gene is different from all the genes, add this gene.
		{
			Gene* g2G = new Gene(*arrG->at(i));
			l->push_front(g2G);
			cnt++;
		}
	}

	for (list<Gene*>::iterator it = l->begin(); it != l->end(); ++it)
	{
		delete (*it);
	}
	delete l;
	return cnt;
}
/*if the schedule is invalid (odd summation in some machine, fix it)*/
void RunAlgorithm::fix_population(vector<Gene*>* arrG)
{
	for (int i = 0; i < population; i++)
		if (!(arrG->at(i)->validGene(machine_number)))
			arrG->at(i)->correctGene(number_of_jobs,machine_number,input_arr);
}
/*Set interval of probabilities for each gene for binary search purposes O(logn) instead linear search.*/
void RunAlgorithm::setLowerUpperProbabilities(vector <Gene*>* arrG, const double& sum_probabilities, const int& iteration_counter)
{
	double run_probability = 0;
	for (int j = 0; j < population; j++)
	{
		arrG->at(j)->setProbabiltyOfGene(arrG->at(j)->getValue() / sum_probabilities);
		arrG->at(j)->setLowerProb(run_probability);
		run_probability += arrG->at(j)->getProbability();
		arrG->at(j)->setUpperProb(run_probability);
		if (print == 1 || (print > 1 && (iteration_counter % print == 0)))
			cout << "Evaluation fo the " << j + 1 << "th gene is " << arrG->at(j)->getValue() << " with probability " << arrG->at(j)->getProbability() << "\n";
	}
}
/*Binary search on interval of probabilities to find the needed gene out of an array*/
int RunAlgorithm::binarySearchForGene(vector<Gene*>* arrG, const double& probablity_randomed)
{
	int lower_bound = 0, upper_bound = population, search_value = (population / 2);
	while (probablity_randomed < arrG->at(search_value)->getLowerProb() || probablity_randomed >arrG->at(search_value)->getUpperProb())
	{
		if (probablity_randomed < arrG->at(search_value)->getLowerProb())
		{
			upper_bound = search_value;
			search_value = ((search_value - lower_bound) / 2);
		}
		else
		{
			lower_bound = search_value;
			search_value += ((upper_bound - search_value) / 2);
		}
	}
	return search_value;
}
/*Main Algorithm, works as follows :
1) evaluate genes based on an evaluation function
2) set probabilities for each gene based on it's evaluation
3) Switching Procedure - create almost all the next population using switching procedure - random 2 probabilities, random 2 indecies and switch between the genes in the
left hand side of the first randomed index and the right hand side of the second randomed index.
4) Mutate Procedure - random probability and take the gene, random some indecies and change the value there again by randomization */
void RunAlgorithm::startAlgo(const int& amount_of_different_genes_to_converge)
{
	vector<Gene*>* array_of_Genes = new vector<Gene*>; //will contain the population of the current generation through out the algorithm
	list<Gene*>* temp_list_of_genes = new list<Gene*>; //will contain the population of the next population through out the algorithm
	initializationOfPopulation(array_of_Genes,temp_list_of_genes);
	int chosenGene = 1;
	size_t chosenGeneG1, chosenGeneG2;
	double  sumProb;
	Gene* g1 = NULL; Gene* g2 = NULL;
	double run_probability, random_first_probability, random_second_probability;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> ds(0.0, 1.0);
	Gene* gBest = NULL; //Best Solution
	int iteration_counter = 0;
	int amount_of_different_genes = amountOfDifferentGenes(array_of_Genes);
	while (amount_of_different_genes > this->converge)//Run till amount of genes will be less than or equal to 5, seeking convergence.
	{
		fix_population(array_of_Genes); //Fixing invalid solutions to valid ones
		//print the population always or every 'x' iterations
		if ((print > 1 && iteration_counter % print == 0) || print == 1)
		{
			cout << "#################################################################################################### \n";
			cout << "Starting Iteration number " << (iteration_counter + 1) << "\n";
			cout << "\n";
			printSolutions(*array_of_Genes);
			cout << "\n";
		}
		//================== Evaluating genes ========================//
		sumProb = 0;
		for (int j = 0; j < population; j++)
		{
			array_of_Genes->at(j)->shortcut(function_number, input_arr, number_of_jobs, sum_of_jobs, machine_number, opt); //Call a function by index 'functionNumber' instead of typing the function.
			if (array_of_Genes->at(j)->getMakespan() == opt) //if found optimal solution
			{
				iteration_counter = -1;
				gBest = array_of_Genes->at(j);
			}
			sumProb += array_of_Genes->at(j)->getValue();
		}
		if (iteration_counter == -1)
			break;

		//================ Setting probabilities =====================//
		if ((print > 1 && (iteration_counter % print == 0)) || print == 1)
			cout << "summation of evaluations is " << sumProb << "\n"; //To check if the algorithm calculates properly
		setLowerUpperProbabilities(array_of_Genes, sumProb, iteration_counter);

		//================= Building next generation =================//
		if (print == 1 || (print > 1 && (iteration_counter % print == 0)))
		{
			cout << "\nStarting switching procedure \n";
			cout << "======================================================================================================= \n";
		}
		//============ Swiching Procedure ===============
		for (int j = 0; j < ((population - mutation) / 2); j++)
		{
			g1 = NULL; g2 = NULL;
			random_first_probability = ds(gen);
			chosenGeneG1 = (binarySearchForGene(array_of_Genes, random_first_probability) + 1); //point to the found gene
			g1 = array_of_Genes->at(chosenGeneG1 - 1);
			random_second_probability = ds(gen);
			while (random_second_probability > g1->getLowerProb() && random_second_probability < g1->getUpperProb())
				random_second_probability = ds(gen);
			chosenGeneG2 = (binarySearchForGene(array_of_Genes, random_second_probability) + 1);
			g2 = array_of_Genes->at(chosenGeneG2 - 1);
			if (print == 1 || (print > 1 && iteration_counter % print == 0))
			{
				cout << "Switching gene number " << chosenGeneG1 << " and " << chosenGeneG2 << "\n";
				cout << "Switching between "; g1->printGene(number_of_jobs); cout << " and "; g2->printGene(number_of_jobs); cout << "\n";
			}
			Gene* g1G = new Gene(*g1);
			Gene* g2G = new Gene(*g2);
			g1G->changeGene(g2G, number_of_jobs,input_arr, print);
			temp_list_of_genes->push_front(g1G);
			temp_list_of_genes->push_front(g2G);
			if ((print > 1 && iteration_counter % print == 0) || print == 1) 
			{
				cout << "Result of switching :\n";
				cout << "Gene number " << chosenGeneG1 << " is "; g1G->printGene(number_of_jobs); cout << "\nGene number " << chosenGeneG2 << " is "; g2G->printGene(number_of_jobs); cout << "\n\n\n";
			}

		}
		//============ Mutation Procedure ===============
		if ((print > 1 && iteration_counter % print == 0) || print == 1) 
		{
			cout << "\n Starting Mutation procedure " << endl;
			cout << "======================================================================================================= \n";
		}
		for (int t = 0; t < mutation; t++)
		{
			random_first_probability = ds(gen);
			chosenGene = (binarySearchForGene(array_of_Genes, random_first_probability) + 1);
			g1 = array_of_Genes->at(chosenGene - 1);
			if (print == 1 || (print > 1 && iteration_counter % print == 0))
			{
				cout << "Chosen gene is number " << (chosenGene + 1) << "\n";
				cout << "Mutating the gene "; (*g1).printGene(number_of_jobs); cout << "\n";
			}
			Gene* g1G = new Gene(*g1);
			g1G->mutate(number_of_jobs, machine_number, input_arr, print, iteration_counter,print);
			temp_list_of_genes->push_front(g1G);
		}
		if (print == 1 || (print > 1 && iteration_counter % print == 0))
			cout << endl << endl;

		//====================== Switching between the current generation and the next one ===================================
		for(int j=0; j<array_of_Genes->size(); j++)
			delete (array_of_Genes->at(j));
		array_of_Genes->clear();
		for (list<Gene*>::iterator it = temp_list_of_genes->begin(); it != temp_list_of_genes->end(); ++it)
		{
			Gene* g1G = new Gene(**it);
			array_of_Genes->insert(array_of_Genes->begin(), g1G);
			delete (*it);
		}
		temp_list_of_genes->clear();
		amount_of_different_genes = amountOfDifferentGenes(array_of_Genes);
		iteration_counter++;
	}

	//fix the solutions in case they are invalid. Find best solution to return to the user.
	fix_population(array_of_Genes);
	if (print > 0)
		cout << "\n\n--------------------------------------------------------------------------- \n\n";
	double max1 = 0;
	if (gBest == NULL) //Didn't find optimal solution through the run.
		for (int j = 0; j < population; j++) //Find best solution by evaluation
		{
			array_of_Genes->at(j)->shortcut(function_number,input_arr,number_of_jobs,sum_of_jobs,machine_number,opt);
			if (array_of_Genes->at(j)->getValue() > max1)
			{
				max1 = array_of_Genes->at(j)->getValue();
				gBest = array_of_Genes->at(j);
			}
		}
	gBest->decodeAndPrint(number_of_jobs,machine_number,input_arr);
	cout << "The makespan is : " << gBest->getMakespan() << endl;
	cout << "OPT is bound by " << opt << endl;
	amount_of_different_genes = amountOfDifferentGenes(array_of_Genes);
	cout << "Amount of different genes at the end is " << amount_of_different_genes << "\n";
	cout << "Amount of iterations done " << iteration_counter << "\n";
	cout << "Deviation from opt is " << (gBest->getMakespan() - opt) << "\n";

	//================== Deletion ===========================//
	for (int j = 0; j < array_of_Genes->size(); j++)
		delete (array_of_Genes->at(j));
	array_of_Genes->clear();
	delete temp_list_of_genes;
	delete array_of_Genes;
}
void RunAlgorithm::printArray(int arr[])
{
	for (int i = 0; i < number_of_jobs; i++)
		cout << "(" << arr[i] << "," << i + 1 << ")\t";
	cout << endl << endl;
}
void RunAlgorithm::printSolutions(vector<Gene*>& arrG)
{
	for (int i = 0; i < population; i++)
	{
		cout << "Number of gene is " << i + 1 << " :   ";
		arrG.at(i)->printGene(number_of_jobs);
		cout << "\n";
	}
}