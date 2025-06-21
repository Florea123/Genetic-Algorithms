#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <map>
#include <cmath>  
#include <vector>

const int popSize = 200;

using namespace std;

vector<vector<int>> createGraph(int n, int m) {
    vector<vector<int>> graph(n, vector<int>(n, 0));  // Initializează un graf cu n noduri, fără muchii

    for (int i = 0; i < m; ++i) {
        int u, v;
        cout << "Enter edge (u v): ";
        cin >> u >> v;
        --u; --v;  // Dacă nodurile sunt indexate de la 1, le scădem cu 1 pentru a lucra pe indexare 0
        graph[u][v] = 1;
        graph[v][u] = 1;  // Adaugă muchia și în sens invers pentru graf neorientat
    }

    return graph;
}

int getMaxColors(const vector<vector<int>>& graph, int n) {
    int maxNumColors = 0;
    for (int i = 0; i < n; ++i) {
        int currMax = 0;
        for (int j = 0; j < n; ++j) {
            currMax += graph[i][j];
        }
        maxNumColors = max(maxNumColors, currMax);
    }
    return maxNumColors;
}

vector<int> createChromosome(int n, int maxNumColors) {
    vector<int> chromosome(n);
    for (int i = 0; i < n; ++i) {
        chromosome[i] = rand() % maxNumColors + 1;
    }
    return chromosome;
}

vector<vector<int>> createPopulation(int n, int maxNumColors) {
    vector<vector<int>> population(popSize);
    for (int i = 0; i < popSize; ++i) {
        population[i] = createChromosome(n, maxNumColors);
    }
    return population;
}

int calcFitness(const vector<vector<int>>& graph, const vector<int>& chromosome, int n) {
    int penalty = 0;
    for (int vertex1 = 0; vertex1 < n; ++vertex1) {
        for (int vertex2 = vertex1; vertex2 < n; ++vertex2) {
            if (graph[vertex1][vertex2] == 1 && chromosome[vertex1] == chromosome[vertex2]) {
                penalty++;
            }
        }
    }
    return penalty;
}

vector<int> simulatedAnnealing(const vector<int>& chromosome, const vector<vector<int>>& graph, int n, int maxNumColors, double T_initial = 10000, double T_min = 1, double alpha = 0.95) {
    vector<int> currentSolution = chromosome;
    int currentFitness = calcFitness(graph, currentSolution, n);
    vector<int> bestSolution = currentSolution;
    int bestFitness = currentFitness;

    double T = T_initial;  // Temperatura inițială

    while (T > T_min) {
        // Generăm o soluție vecină prin schimbarea unei culori
        vector<int> neighbor = currentSolution;
        int position = rand() % n;
        neighbor[position] = rand() % maxNumColors + 1;

        int neighborFitness = calcFitness(graph, neighbor, n);
        int delta = neighborFitness - currentFitness;

        if (delta < 0) {
            // Dacă soluția vecină este mai bună, o acceptăm
            currentSolution = neighbor;
            currentFitness = neighborFitness;
            if (currentFitness < bestFitness) {
                bestSolution = currentSolution;
                bestFitness = currentFitness;
            }
        }
        else {
            // Dacă soluția vecină este mai rea, o acceptăm cu o probabilitate
            double probability = exp(-delta / T);
            if (rand() / double(RAND_MAX) < probability) {
                currentSolution = neighbor;
                currentFitness = neighborFitness;
            }
        }

        T *= alpha;  // Răcirea (scăderea temperaturii)
    }

    return bestSolution;
}

vector<vector<int>> tournamentSelection(vector<vector<int>>& population, const vector<vector<int>>& graph, int n) {
    vector<vector<int>> newPopulation;
    for (int k = 0; k < 2; ++k) {
        random_shuffle(population.begin(), population.end());
        for (int i = 0; i < popSize - 1; i += 2) {
            if (calcFitness(graph, population[i], n) < calcFitness(graph, population[i + 1], n)) {
                newPopulation.push_back(population[i]);
            }
            else {
                newPopulation.push_back(population[i + 1]);
            }
        }
    }
    return newPopulation;
}

pair<vector<int>, vector<int>> twoPointCrossover(const vector<int>& parent1, const vector<int>& parent2) {
    int firstPoint = rand() % (parent1.size() - 2) + 1;
    int secondPoint = rand() % (parent1.size() - firstPoint - 1) + firstPoint + 1;

    vector<int> child1(parent1.begin(), parent1.begin() + firstPoint);
    for (int i = firstPoint; i < secondPoint; ++i) {
        child1.push_back(parent2[i]);
    }
    for (int i = secondPoint; i < parent1.size(); ++i) {
        child1.push_back(parent1[i]);
    }

    vector<int> child2(parent2.begin(), parent2.begin() + firstPoint);
    for (int i = firstPoint; i < secondPoint; ++i) {
        child2.push_back(parent1[i]);
    }
    for (int i = secondPoint; i < parent2.size(); ++i) {
        child2.push_back(parent2[i]);
    }

    return { child1, child2 };
}

vector<int> mutation(vector<int>& chromosome, double chance, const vector<vector<int>>& graph, int maxNumColors, int n) {
    double possible = (double)rand() / RAND_MAX;
    if (chance <= possible) {
        for (int vertex1 = 0; vertex1 < n; ++vertex1) {
            for (int vertex2 = vertex1; vertex2 < n; ++vertex2) {
                if (graph[vertex1][vertex2] == 1 && chromosome[vertex1] == chromosome[vertex2]) {
                    chromosome[vertex1] = rand() % maxNumColors + 1;
                }
            }
        }
    }
    return chromosome;
}

int maxColoring(const map<int, int>& colorDict) {
    int maxColors = 0;
    for (map<int, int>::const_iterator it = colorDict.begin(); it != colorDict.end(); ++it) {
        maxColors = max(maxColors, it->second);
    }
    return maxColors;
}

int testColoring(const vector<vector<int>>& graph, int n) {
    map<int, int> colorGraph;

    // Parcurge fiecare vertex din graf
    for (int vertex = 0; vertex < graph.size(); ++vertex) {
        vector<bool> unusedColors(graph.size(), true);  // Marcăm toate culorile ca disponibile

        // Parcurge vecinii unui vertex
        for (int neighbor = 0; neighbor < graph.size(); ++neighbor) {
            if (graph[vertex][neighbor] == 1 && colorGraph.find(neighbor) != colorGraph.end()) {
                unusedColors[colorGraph[neighbor]] = false;  // Marcam culoarea ca fiind folosită
            }
        }

        // Atribuim prima culoare disponibilă
        for (int i = 0; i < unusedColors.size(); ++i) {
            if (unusedColors[i]) {
                colorGraph[vertex] = i;
                break;
            }
        }
    }
    return maxColoring(colorGraph);  // Calculăm numărul maxim de culori
}

double calculateMean(const vector<int>& fitnesses) {
    double sum = accumulate(fitnesses.begin(), fitnesses.end(), 0);
    return sum / fitnesses.size();
}

double calculateStandardDeviation(const vector<int>& fitnesses, double mean) {
    double variance = 0;
    for (int fitness : fitnesses) {
        variance += (fitness - mean) * (fitness - mean);
    }
    return sqrt(variance / fitnesses.size());
}

int main() {
    int n, m;
    cout << "Enter number of nodes (n): ";
    cin >> n;
    cout << endl;
    cout << "Enter number of edges (m): ";
    cin >> m;
    cout << endl;
    vector<int> fitnesses; // vector pentru a salva fitness-urile
    int maxim;

    vector<vector<int>> graph = createGraph(n, m);
    for (int run = 0; run < 30; run++) {
        cout << "Run: " << run << endl;
        int maxNumColors = getMaxColors(graph, n);

        int checkCount = 0;
        int failedColors = 0;
        maxim = INFINITY;
        cout << "Trying to color with " << maxNumColors << " colors\n";

        while (true) {
            vector<vector<int>> population = createPopulation(n, maxNumColors);

            int bestFitness = calcFitness(graph, population[0], n);
            vector<int> fittest = population[0];

            int generation = 0;
            int numGenerations = 1000;

            while (bestFitness != 0 && generation != numGenerations) {
                generation++;

                population = tournamentSelection(population, graph, n);

                if (population.size() % 2 != 0) {
                    population.pop_back();
                }

                vector<vector<int>> childrenPopulation;
                random_shuffle(population.begin(), population.end());
                for (size_t i = 0; i < population.size() - 1; i += 2) {
                    pair<vector<int>, vector<int>> children = twoPointCrossover(population[i], population[i + 1]);
                    childrenPopulation.push_back(children.first);
                    childrenPopulation.push_back(children.second);
                }

                for (size_t i = 0; i < childrenPopulation.size(); ++i) {
                    mutation(childrenPopulation[i], 0.007, graph, maxNumColors, n);
                    childrenPopulation[i] = simulatedAnnealing(childrenPopulation[i], graph, n, maxNumColors);
                }
                for (size_t i = population.size(); i < popSize; ++i) {
                    population.push_back(createChromosome(n, maxNumColors));
                }

                population = childrenPopulation;
                bestFitness = calcFitness(graph, population[0], n);
                fittest = population[0];

                for (size_t i = 0; i < population.size(); ++i) {
                    int fitness = calcFitness(graph, population[i], n);
                    if (fitness < bestFitness) {
                        bestFitness = fitness;
                        fittest = population[i];
                        if (maxim > bestFitness)
                            maxim = bestFitness;
                    }
                }

                if (bestFitness == 0) {
                    break;
                }
            }

            if (bestFitness == 0) {
                maxim = bestFitness;
                std::cout << maxNumColors << " colors succeeded! Trying " << maxNumColors - 1 << " colors\n";
                maxNumColors--;
                checkCount = 0;
            }
            else {
                if (checkCount != 2 && maxNumColors > 1) {
                    failedColors = maxNumColors;

                    if (checkCount == 0) {
                        std::cout << maxNumColors << " failed. For safety, checking for improvement with " << maxNumColors << " colors again\n";
                    }
                    if (checkCount == 1) {
                        std::cout << maxNumColors << " failed. For safety, checking for improvement with " << maxNumColors - 1 << " colors\n";
                        maxNumColors--;
                    }

                    checkCount++;
                    continue;
                }
                if (maxNumColors > 1) {
                    std::cout << "Graph is " << failedColors + 1 << " colorable\n";
                }
                else {
                    std::cout << "Graph is " << maxNumColors + 1 << " colorable\n";
                }
                std::cout << "Test Solution: " << testColoring(graph, n) << std::endl;
                break;
            }
        }
        fitnesses.push_back(maxim);
    }
    double mean = calculateMean(fitnesses);
    double stddev = calculateStandardDeviation(fitnesses, mean);
    int minFitness = *min_element(fitnesses.begin(), fitnesses.end());
    int maxFitness = *max_element(fitnesses.begin(), fitnesses.end());

    cout << "\nFitness statistics after 30 runs:\n";
    cout << "Mean: " << mean << endl;
    cout << "Standard Deviation: " << stddev << endl;
    cout << "Min Fitness: " << minFitness << endl;
    cout << "Max Fitness: " << maxFitness << endl;

    return 0;
}