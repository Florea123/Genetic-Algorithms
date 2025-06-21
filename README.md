# Genetic-Algorithms

This repository contains implementations of various optimization algorithms designed to solve mathematical and graph-based problems. The algorithms include Hill Climbing, Simulated Annealing, and Genetic Algorithms, applied to both mathematical functions and graph coloring problems.

## Overview

### Hill Climbing and Simulated Annealing
This module implements two optimization techniques:
- **Hill Climbing**: A local search algorithm with three improvement strategies:
  - Best Improvement
  - First Improvement
  - Worst Improvement
- **Simulated Annealing**: A probabilistic technique for approximating the global optimum of a given function.

Supported functions:
- De Jong
- Schwefel
- Rastrigin
- Michalewicz

### Genetic Algorithm for Mathematical Functions
This module applies Genetic Algorithms to optimize mathematical functions. It includes:
- **Genetic Operators**:
  - Mutation
  - Crossover
  - Selection
- **Supported Functions**:
  - De Jong
  - Schwefel
  - Rastrigin
  - Michalewicz

### Genetic Algorithm Applied to Graph Coloring
This module uses Genetic Algorithms to solve graph coloring problems, determining the chromatic number of a graph. Features include:
- **Graph Coloring**:
  - Chromatic number calculation
- **Genetic Operators**:
  - Tournament Selection
  - Two-Point Crossover
  - Mutation
- **Simulated Annealing**: Used for local optimization.

## Features
- **Optimization Algorithms**:
  - Hill Climbing
  - Simulated Annealing
  - Genetic Algorithm
- **Applications**:
  - Mathematical function optimization
  - Graph coloring problems
- **Statistical Analysis**:
  - Mean
  - Standard deviation
  - Minimum and maximum fitness values

## How to Use
1. Clone the repository:
   ```bash
   git clone https://github.com/your-repo/genetic-algorithms.git
   cd genetic-algorithms

2. Compile the .cpp files using a C++ compiler (e.g., g++):
```bash
    g++ -o HillClimbing Hill\ Climbing\ and\ Simulated\ Annealing/H1.cpp
    g++ -o GeneticAlgorithm Genetic\ Algorithm/H2.1.cpp
    g++ -o GraphColoring Genetic\ Algorithm\ Applied\ to\ Graphs/H3.cpp
```
3. Run the executables:
    - **Hill Climbing and Simulated Annealing**:
      ```bash
      ./HillClimbing
      ```
    - **Genetic Algorithm for Mathematical Functions**:
      ```bash
      ./GeneticAlgorithm
      ```
    - **Genetic Algorithm for Graph Coloring**:
      ```bash
      ./GraphColoring
      ```

## Folder Structure
- **Hill Climbing and Simulated Annealing**:
  - Contains `H1.cpp`, implementing Hill Climbing and Simulated Annealing algorithms.
- **Genetic Algorithm**:
  - Contains `H2.1.cpp`, implementing Genetic Algorithms for mathematical functions.
- **Genetic Algorithm Applied to Graphs**:
  - Contains `H3.cpp`, implementing Genetic Algorithms for graph coloring problems.


## Dependencies
- C++ Standard Library
- Random number generation utilities (`<random>`)
- Chrono utilities for timing (`<chrono>`)