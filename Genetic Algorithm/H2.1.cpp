#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <thread>

using namespace std;
unsigned seed = chrono::system_clock::now().time_since_epoch().count();
mt19937_64 generator(seed);
mt19937_64 g_randomGenerator;
const int m = 100;
const double m_p = 0.001, c_p = 0.7;
int n, d = 5, N, l, L, f;
double a, b, minimum_global;

void initializare()
{
	N = (b - a) * pow(10, d);
	l = ceil(log2(N));
	L = l * n;
}
double DeJong(vector<double> v)
{
	double sum = 0;
	int j=0;
	for (int i = 1; i <= n; i++,j++)
	{
		sum = sum + v[j] * v[j];
	}
	return sum;
}
double Schwefel(vector<double> v)
{
	double sum = 0;
	int j=0;
	for (int i = 1; i <= n; i++,j++)
	{
		if (v[j] < 0)
		{
			v[j] = v[j] * (-1);
		}
		sum = sum + v[j] * sin(sqrt(v[j]));
	}
	return 418.9829 * n - sum;
}
double Rastrigin(vector<double> v)
{
	double sum = 0;
	int j=0;
	for (int i = 1; i <= n; i++,j++)
	{
		sum = sum + v[j] * v[j] - 10 * cos(2 * 3.14 * v[j]);
		//atan(1)=arctg(1)=pi/4
	}
	return 10 * n + sum;
}
double Michalewicz(vector<double> v)
{
	double sum = 0;
	int j=0;
	for (int i = 1; i <= n; i++,j++)
	{
		sum = sum + sin(v[j]) * pow(sin((i * v[j] * v[j]) / 3.14), 20);
	}
	return (-1) * sum;
}
int numar_random(int x, int y) {
	std::random_device rd;  
	std::mt19937 gen(rd()); 
	std::uniform_int_distribution<> distrib(x, y);

	return distrib(gen);
}
/// Functia de pe randul 73 a fost preluata de pe chatgpt pentru a genera un numar random intr-un interval [x,y]
vector<double> valoare_pop_linie(vector<bool> pop_linie)
{
	vector<double> pop_linie_rezultat(L);

	for (int i = 0; i < n; i++)
	{
		double rezultat = 0;
		for (int j = l * i; j < l * i + l; j++)
		{
			if (pop_linie[j] == 1) 
			{
				rezultat = rezultat + pow(2, l - 1 - j + l * i); 
			}
		}
		rezultat = a + rezultat * (b - a) / (pow(2, l) - 1);
		pop_linie_rezultat[i] = rezultat;
	}

	return pop_linie_rezultat;
}
double rezultat_linie(vector<bool> pop_linie)
{
	vector<double> pop_linie_rezultat(L);
	pop_linie_rezultat = valoare_pop_linie(pop_linie);

	if (f == 1)
	{
		return DeJong(pop_linie_rezultat);
	}
	if (f == 2)
	{
		return Schwefel(pop_linie_rezultat);
	}
	if (f == 3)
	{
		return Rastrigin(pop_linie_rezultat);
	}
	if (f == 4)
	{
		return Michalewicz(pop_linie_rezultat);
	}
}
vector<vector<bool>> generator_populatie()
{
	vector<vector<bool>> pop(m);
	for (int i = 0; i < m; i++)
		pop[i] = vector<bool>(L);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < L; j++)
			pop[i][j] = 0;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < L; j++)
			pop[i][j] = numar_random(0, 1);
	}

	return pop;
}
vector<double> evaluate(vector<vector<bool>> pop)
{
	vector<double> pop_rezultat(m);

	for (int i = 0; i < m; i++)
	{
		pop_rezultat[i] = rezultat_linie(pop[i]);
	}

	return pop_rezultat;
}
vector<vector<bool>> mutate(vector<vector<bool>> pop)
{
	double p;
	for (int i = 0; i < m; i++) 
	{
		for (int j = 0; j < L; j++)
		{
			p = (double)numar_random(1, 99999) / 100000;
			if (p < m_p)
			{
				if (pop[i][j] == 1)
				{
					pop[i][j] = 0;
				}
				else
				{
					pop[i][j] = 1;
				}
			}
		}
	}
	return pop;
}
vector<vector<bool>> cross_over(vector<vector<bool>> pop)
{
	vector<vector<bool>> pop2(m);
	for (int i = 0; i < m; i++)
		pop2[i] = vector<bool>(L);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < L; j++)
			pop2[i][j] = 0;
	vector<bool> parinte1(L), parinte2(L), copil1(L), copil2(L);
	vector<double> p(m);
	int line = 0, pp, nr = 0, size2 = 0;
	int rand, taietura1, taietura2;
	double rezultat_parinte1, rezultat_parinte2, rezultat_copil1, rezultat_copil2;

	for (int i = 0; i < m; i++)
	{
		p[i] = (double)numar_random(1, 99999) / 100000;
	}

	while (p[line] < c_p && line < m)
	{
		line++;
	}
	if (line % 2 != 0) 
	{
		pp = (double)numar_random(1, 99999) / 100000;
		if (pp < 0.5) 
		{
			line++; 
		}
		else
		{
			line--; 
		}
	}

	while (nr < line)
	{
		parinte1= pop[nr];
		parinte2 = pop[nr + 1];
		nr = nr + 2;
		rand = numar_random(1, L - 2);
		taietura1 = numar_random(1, rand);//taietura va fi la stanga indexului dat
		taietura2 = numar_random(taietura1 + 1, L - 1); //taietura va fi la stanga indexului dat

		for (int i = 0; i < taietura1; i++) //face incrucisarea (in zigzag)
		{
			copil1[i] = parinte1[i];
			copil2[i] = parinte2[i];
		}
		for (int i = taietura1; i < taietura2; i++)
		{
			copil1[i] = parinte2[i];
			copil2[i] = parinte1[i];
		}
		for (int i = taietura2; i < L; i++)
		{
			copil1[i] = parinte1[i];
			copil2[i] = parinte2[i];
		}

		rezultat_parinte1 = rezultat_linie(parinte1);
		rezultat_parinte2 = rezultat_linie(parinte2);
		rezultat_copil1 = rezultat_linie(copil1);
		rezultat_copil2 = rezultat_linie(copil2);

		if (rezultat_parinte1 >= rezultat_copil1) 
		{
			pop2[size2] = copil1;
			size2++;
		}
		else
		{
			pop2[size2] = parinte1;
			size2++;
		}

		if (rezultat_parinte2 >= rezultat_copil2) 
		{
			pop2[size2] = copil2;
			size2++;
		}
		else
		{
			pop2[size2] = parinte2;
			size2++;
		}
	}

	vector<double> pop_rezultat(m), pop_rezultat_sortat(m);
	pop_rezultat_sortat = pop_rezultat;
	sort(pop_rezultat_sortat.begin(), pop_rezultat_sortat.end());

	vector<int> ver(m);
	for (int i = 0; i < m; i++)
	{
		ver[i] = 0; 
	}
	for (int i = size2; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (ver[j] != 0)
			{
				continue;
			}

			if (pop_rezultat_sortat[i] == pop_rezultat[j])
			{
				ver[j] = 1;
				pop2[i] = pop[j];
				break;
			}
		}
	}

	return pop2;
}
vector<vector<bool>> select(vector<vector<bool>> pop)
{
	vector<vector<bool>> pop2(m);
	for (int i = 0; i < m; i++)
		pop2[i] = vector<bool>(L);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < L; j++)
			pop2[i][j] = 0;

	vector<double> pop_rezultat(m), pop_p(m), pop_ap(m);
	pop_rezultat = evaluate(pop);
	double total = 0, p;
	int k;

	vector<double> pop_rezultat_sortat(m);
	pop_rezultat_sortat = pop_rezultat;
	sort(pop_rezultat_sortat.begin(), pop_rezultat_sortat.end()); 
	for (int i = 0; i < m; i++)
	{
		total += pop_rezultat[i];
	}
	for (int i = 0; i < m; i++)
	{
		pop_p[i] = pop_rezultat[i] / total;
	}
	pop_ap[0] = 0;
	for (int i = 0; i < m - 1; i++)
	{
		pop_ap[i + 1] = pop_ap[i] + pop_p[i];
	}
	vector<int> ver(m);
	for (int i = 0; i < m; i++)
	{
		ver[i] = 0; 
	}

	for (int i = 0; i < 4; i++) 
	{
		for (int j = 0; j < m; j++)
		{
			if (ver[j] != 0)
			{
				continue;
			}

			if (pop_rezultat_sortat[i] == pop_rezultat[j])
			{
				ver[j] = 1;
				pop2[i] = pop[j];
				break;
			}
		}
	}
	///Roata norocului
	for (int i = 4; i < m; i++)
	{
		p = (double)numar_random(1, 99999) / 100000;
		k = 0;
		for (k = 0; k < m - 1; k++)
		{
			if (pop_ap[k] < p && p <= pop_ap[k + 1])
			{
				break;
			}
		}

		ver[k] = 1; 
		pop2[i] = pop[k];
	}

	return pop2;
}
double alegere_rezultat(vector<double> pop_rezultat)
{
	double rez = pop_rezultat[0], minim = 1000, t;
	vector<double> pop_rezultat_sortat(m);
	pop_rezultat_sortat = pop_rezultat;
	sort(pop_rezultat_sortat.begin(), pop_rezultat_sortat.end());

	for (int i = 0; i < pop_rezultat.size(); i++)
	{
		t = pop_rezultat[i] - minimum_global;
		if (t < 0)
		{
			t = t * (-1);
		}
		if (minim > t)
		{
			minim = t;
		}
	}

	for (int i = 0; i < pop_rezultat.size(); i++)
	{
		t = pop_rezultat_sortat[i] - minimum_global;
		if (t < 0)
		{
			t = t * (-1);
		}
		if (minim == t)
		{
			rez = pop_rezultat_sortat[i];
			break;
		}
	}

	return rez;
}
void algoritm_genetic()
{
	int nr = 0;
	double rezultat;

	vector<vector<bool>> pop(m);
	for (int i = 0; i < m; i++)
		pop[i] = vector<bool>(L);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < L; j++)
			pop[i][j] = 0;
	pop = generator_populatie();

	vector<double> pop_rezultat(m);

	pop_rezultat = evaluate(pop);

	while (nr < 1000)
	{
		nr++;

		pop = mutate(pop);
		pop = cross_over(pop);
		pop = select(pop);
		pop_rezultat = evaluate(pop);
	}
	rezultat = alegere_rezultat(pop_rezultat);
	cout << "Rezultat: " << rezultat << endl;
}
int main()
{
	cout << fixed << setprecision(d);
	cout << "Pentru functia DeJong:1, pentru functia Schwefel:2, pentru functia Rastrigin:3, pentru functia Michalewicz:4 :";
	cin >> f;
	cout << '\n' << "Numarul de dimensiuni: ";
	cin >> n;
	if (f == 1)
	{
		a = -5.12;
		b = 5.12;
		minimum_global = 0;

		initializare();

		for (int i = 1; i < 31; i++)
		{
			cout << i << " : ";
			algoritm_genetic();
		}
	}
	else
		if (f == 2)
		{
			a = -500;
			b = 500;
			minimum_global = 0;

			initializare();

			for (int i = 1; i < 31; i++)
			{
				cout << i << " : ";
				algoritm_genetic();
			}
		}
	if (f == 3)
	{
		a = -5.12;
		b = 5.12;
		minimum_global = 0;

		initializare();

		for (int i = 1; i < 31; i++)
		{
			cout << i << " : ";
			algoritm_genetic();
		}
	}
	else
		if (f == 4)
		{
			a = 0;
			b = atan(1) * 4;
			if (n == 5)
			{
				minimum_global = -4.68765;
			}
			if (n == 10)
			{
				minimum_global = -9.66015;
			}
			if (n == 30)
			{
				minimum_global = -29.63088;
			}
			initializare();

			for (int i = 1; i < 31; i++)
			{
				cout << i << " : ";
				algoritm_genetic();
			}
		}
		else
			cout << "Nu ati ales nici o functie!";
	/*
	cout << "Dimensiunea" << n << endl;
	cout << "DeJong" << '\n';
	a = -5.12;
	b = 5.12;
	f = 1;
	minimum_global = 0;

	initializare();

	for (int i = 1; i < 31; i++)
	{
		cout << i << " : ";
		algoritm_genetic();
	}
	cout << "Schwefel" << '\n';
	a = -500;
	b = 500;
	f = 2;
	minimum_global = 0;

	initializare();

	for (int i = 1; i < 31; i++)
	{
		cout << i << " : ";
		algoritm_genetic();
	}
	cout << "Rastrigin" << '\n';
	a = -5.12;
	b = 5.12;
	f = 3;
	minimum_global = 0;

	initializare();

	for (int i = 1; i < 31; i++)
	{
		cout << i << " : ";
		algoritm_genetic();
	}
	cout << "Michalewicz" << '\n';
	a = 0;
	b = atan(1) * 4;
	f = 4;
	if (n == 5)
	{
		minimum_global = -4.68765;
	}
	if (n == 10)
	{
		minimum_global = -9.66015;
	}
	if (n == 30)
	{
		minimum_global = -29.63088;
	}
	initializare();

	for (int i = 1; i < 31; i++)
	{
		cout << i << " : ";
		algoritm_genetic();
	}
	*/
	return 0;
}

