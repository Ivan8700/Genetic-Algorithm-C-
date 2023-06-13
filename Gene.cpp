#include "Gene.h"

Gene::Gene(int n = 0, int m = 0)
{
	this->length = n;
	this->machines = m;
	this->gene = new int[n];
	this->sum_array = new int[m];
	for (int i = 0; i < m; i++)
		sum_array[i] = 0;
	makespan = 0;
	value = 0;
	probability = 0;
	lowerProb = 0;
	upperProb = 0;
}
Gene::Gene(const Gene& g)
{
	this->length = g.length;
	this->machines = g.machines;
	this->gene = new int[g.length];
	this->sum_array = new int[g.machines];
	for (int i = 0; i < g.length; i++)
		this->gene[i] = g.gene[i];
	for (int i = 0; i < g.machines; i++)
		this->sum_array[i] = g.sum_array[i];
	this->makespan = g.makespan;
	this->value = g.value;
	this->probability = g.probability;
	this->lowerProb = g.lowerProb;
	this->upperProb = g.upperProb;
}
Gene::~Gene()
{
	delete[] this->gene;
	delete[] this->sum_array;
}

int Gene::getMakespan() { return this->makespan; }
double Gene::getLowerProb() { return this->lowerProb; }
double Gene::getUpperProb() { return this->upperProb; }
void Gene::setLowerProb(const double& k) { this->lowerProb = k; }
void Gene::setUpperProb(const double& k) { this->upperProb = k; }
int Gene::getGeneAtIndex(int i) { return this->gene[i]; }
/*Randomization of the gene (used only in initialization)*/
void Gene::randomGene(const int& n, const int& m, const int arr[])
{
	int k;
	for (int i = 0; i < n; i++)
	{
		k = rand() % m;
		this->setGeneAtIndex(i, k);
		this->sum_array[k] += arr[i];
	}
}
/*Check if the schedule of the gene is invalid*/
bool Gene::validGene(int m)
{
	bool ok = true;
	for (int i = 0; i < m; i++)
		if (sum_array[i] % 2 == 1)
			ok = false;
	return ok;
}
/*update the sum_array of the gene*/
void Gene::updateGene(const int& n, const int& m, const int arr[])
{
	for (int i = 0; i < m; i++)
		this->sum_array[i] = 0;
	for (int i = 0; i < n; i++)
	{
		this->sum_array[this->getGeneAtIndex(i)] += arr[i];
	}
}

void Gene::setGeneAtIndex(const int& index, const int& value_to_set)
{
	this->gene[index] = value_to_set;
}

int* Gene::getGene() { return this->gene; }
void Gene::printGene(const int& lengthOfGene)
{
	for (int i = 0; i < lengthOfGene; i++)
		cout << this->gene[i];
}

double Gene::getValue() { return this->value; }
double Gene::getProbability() { return this->probability; }
void Gene::setProbabiltyOfGene(double q) { this->probability = q; }
/*running the required function, Y= makespan X=opt bound.
The functions are :
1) 1/(Y-X+1)
2) 1/sum_on_all_machine_of(load on each machine - opt)
3) Same as function number 2 but (load on each machine - opt)^(0.5)
4) 1/(Y-M) where M is the minimal machine load.*/
void Gene::shortcut(const int& function_number_to_run, const int arr[], const int& n, const int& sum, const int& m, const int& opt)
{
	if (function_number_to_run == 1)
		valueGeneByMakespanMinusOPTplusOne(arr, n, sum, m);
	if (function_number_to_run == 2)
		valueGeneByDeviationFromOPTfunction(arr, n, sum, m, opt);
	if (function_number_to_run == 3)
		valueGeneByDeviationFromOPTsquaredFunction(arr, n, sum, m, opt);
	if (function_number_to_run == 4)
		valueGeneByMakespanMinusMinimalFunction(arr, n, sum, m, opt);
}

void Gene::valueGeneByDeviationFromOPTfunction(const int arr[], const int& length, const int& sum, const int& m, const int& opt)
{
	this->makespan = 0;
	int sump = 0;
	this->updateGene(length, m, arr);
	for (int i = 0; i < m; i++)
	{
		if (this->sum_array[i] > this->makespan)
			this->makespan = this->sum_array[i];
		sump += abs(this->sum_array[i] - opt);
	}
	if (sump == 0)
		this->value = INT_MAX;
	else
		this->value = (1.0) / (sump);
}

void Gene::valueGeneByDeviationFromOPTsquaredFunction(const int arr[], const int& n, const int& sum, const int& m, const int& opt)
{
	{
		this->makespan = 0;
		int sump = 0;
		this->updateGene(length, m, arr);
		for (int i = 0; i < m; i++)
		{
			if (this->sum_array[i] > this->makespan)
				this->makespan = this->sum_array[i];
			sump += pow(abs((this->sum_array[i] - opt)), 0.5);
		}
		if (sump == 0)
			this->value = INT_MAX;
		else
			this->value = (1.0) / (sump);
	}
}

void Gene::valueGeneByMakespanMinusMinimalFunction(const int arr[], const int& length, const int& sum, const int& m, const int& opt)
{
	this->makespan = 0;
	int minimal_load_machine = INT_MAX;
	this->updateGene(length, m, arr);
	for (int i = 0; i < m; i++)
	{
		if (this->sum_array[i] > this->makespan)
			this->makespan = this->sum_array[i];
		if (this->sum_array[i] < minimal_load_machine)
			minimal_load_machine = this->sum_array[i];
	}
	if ((this->makespan - minimal_load_machine) == 0)
		this->value = INT_MAX;
	else
		this->value = (1.0) / (this->makespan - minimal_load_machine);
}

void Gene::valueGeneByMakespanMinusOPTplusOne(const int arr[], const int& length, const int& sum, const int& m)
{
	this->makespan = 0;
	this->updateGene(length, m, arr);
	for (int i = 0; i < m; i++)
	{
		if (this->sum_array[i] > this->makespan)
			this->makespan = this->sum_array[i];
	}
	int opt = ceil((double)sum / m);
	opt = (opt % 2 == 1) ? opt + 1 : opt;
	this->value = (1.0) / ((this->makespan - opt) + 1);
}
/*Switching function of 2 genes.*/
void Gene::changeGene(Gene* g, const int& n, const int arr[], const int& print)
{
	int first_index_randomed = rand() % ((n / 2) - 1);
	int second_index_randomed = (rand() % (n / 2)) + (n / 2) + 1;
	if (second_index_randomed >= n)
		second_index_randomed = (n - 1);
	if (print == 1)
		cout << "Changing the first " << first_index_randomed + 1 << " coordinates and the last " << (n - second_index_randomed) << " coordinates\n";
	int* temp = new int[n];
	for (int i = 0; i < n; i++)
	{
		if (i <= first_index_randomed || i >= second_index_randomed)
		{
			temp[i] = g->getGeneAtIndex(i);
			this->sum_array[temp[i]] += arr[i];
			this->sum_array[this->getGeneAtIndex(i)] -= arr[i];
			g->sum_array[temp[i]] -= arr[i];
			g->sum_array[this->getGeneAtIndex(i)] += arr[i];
			g->setGeneAtIndex(i, this->gene[i]);
			this->gene[i] = temp[i];
		}
	}
	delete[] temp;
}
/*Decode the gene and print the schedule*/
void Gene::decodeAndPrint(int n, int m, const int arr[])
{
	list<std::pair<int, int>>* schedule = new list<std::pair<int, int>>[m];
	for (int i = 0; i < n; i++)
	{
		schedule[this->getGeneAtIndex(i)].push_front(std::make_pair(arr[i], i + 1));
	}
	for (int i = 0; i < m; i++)
	{
		cout << "\nMachine number " << (i + 1) << "\n";
		cout << "{";
		for (list<std::pair<int, int>>::iterator it = schedule[i].begin(); it != schedule[i].end(); ++it)
		{
			cout << "(" << (*it).first << "," << (*it).second << ") ";
		}
		cout << "} \n";
		cout << "->Sum of The Machine is : " << this->sum_array[i] << "\n \n";
	}
	delete[] schedule;
}
/*Mutation procedure, random 2% of the length of the gene and change the values in the indecies.*/
void Gene::mutate(const int& length, const int& m, const int arr[], const int& print, const int& iteration_counter, const int& print_every_x_iterations)
{
	int index_randomed;
	int randomed_machine_index;
	int* placeHolder = new int[length / 50]; //mutate only if the amount of jobs is bigger than 50, otherwise do nothing.
	bool ok = true;
	for (int i = 0; i < length / 50; i++)
	{
		index_randomed = (rand() % length);
		placeHolder[i] = index_randomed;
		// Make sure the randomed index is different.
		for (int j = 0; j < i; j++)
			if (placeHolder[j] == index_randomed)
				ok = false;
		while (!ok)
		{
			ok = true;
			index_randomed = (rand() % length);
			for (int j = 0; j < i; j++)
				if (placeHolder[j] == index_randomed)
					ok = false;
		}

		randomed_machine_index = (rand() % m);
		while (randomed_machine_index == this->getGeneAtIndex(index_randomed)) //make sure we change the value in the randomed index.
			randomed_machine_index = rand() % m;
		this->sum_array[this->getGeneAtIndex(index_randomed)] -= arr[index_randomed];
		this->sum_array[randomed_machine_index] += arr[index_randomed];
		if (print == 1 || (print > 1 && iteration_counter % print_every_x_iterations == 0))
		{
			cout << "at index " << index_randomed + 1 << " to " << randomed_machine_index << "\n";
		}
		this->setGeneAtIndex(index_randomed, randomed_machine_index);
	}
	if (print == 1 || (print > 1 && iteration_counter % print == 0))
	{
		cout << "Result after the mutation is "; this->printGene(length); cout << "\n\n";
	}
	delete[] placeHolder;
}
/*Correcting invalid solutions (odd summation at some machine) */
void Gene::correctGene(const int& n, const int& m, const int arr[])
{
	int index_min_odd, index_max_odd, minimum_odd_job_at_maximum_machine, index_of_job, maximal_machine_odd, minimal_machine_odd;
	this->makespan = 0;
	//this->updateGene(n, m, arr);
	while (!this->validGene(m)) //while solution is invalid
	{
		index_min_odd = 0; index_max_odd = 0, maximal_machine_odd = 0, minimal_machine_odd = 999999;
		for (int i = 0; i < m; i++)
		{
			if (sum_array[i] < minimal_machine_odd && sum_array[i] % 2 == 1)
			{
				index_min_odd = i;
				minimal_machine_odd = sum_array[i];
			}
		}
		for (int i = 0; i < m; i++)
		{
			if (sum_array[i] > maximal_machine_odd&& sum_array[i] % 2 == 1)
			{
				index_max_odd = i;
				maximal_machine_odd = sum_array[i];
			}
		}
		if (index_min_odd == index_max_odd)
		{
			for (int i = 0; i < m; i++)
				if (i != index_max_odd && sum_array[i] == sum_array[index_min_odd])
				{
					index_min_odd = i;
					break;
				}
		}
		minimum_odd_job_at_maximum_machine = INT_MAX;
		for (int i = 0; i < n; i++)
		{
			if (this->getGeneAtIndex(i) == index_max_odd && arr[i] % 2 == 1 && arr[i] < minimum_odd_job_at_maximum_machine)
			{
				minimum_odd_job_at_maximum_machine = arr[i];
				index_of_job = i;
			}
		}
		if (index_min_odd == index_max_odd)
		{

			for (int i = 0; i < n; i++)
			{
				if (this->getGeneAtIndex(index_min_odd) && arr[i] % 2 == 1 && arr[i] < minimum_odd_job_at_maximum_machine)
				{
					minimum_odd_job_at_maximum_machine = arr[i];
				}
			}
		}
		sum_array[index_max_odd] -= minimum_odd_job_at_maximum_machine;
		sum_array[index_min_odd] += minimum_odd_job_at_maximum_machine;
		this->setGeneAtIndex(index_of_job, index_min_odd);
	}

};