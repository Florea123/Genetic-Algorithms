#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <ctime>
#include <thread>
#include <chrono>
#include <string>
#include <vector>
#include <utility>
using namespace std;
bool v[10001], v_vecin[10001], aux[30], aux_vecin[30];
double valori[101];
int kk = 0;
mt19937_64 g_randomGenerator;
void initRandomGenerator(mt19937_64& generator)
{
    mt19937_64 helper;
    helper.seed(time(NULL) + clock() * 1000 + hash<thread::id>{}(this_thread::get_id()) * 4 * 100000);
    helper.discard(31337);
    generator.seed(helper());
}///Functie preluata de la domnul profesor Eugen Croitoru
double get_value(const bool v[], const int poz_init, const int poz_end)
{
    int p = 0;
    double s = 0;
    for (int i = poz_end - 1; i >= poz_init; i--)
    {
        if (p == 0)
        {
            s = s + v[i];
            p = 2;
        }
        else
        {
            s = s + v[i] * p;
            p = p * 2;
        }
    }
    if (v[poz_end] == 1)
        s = s * (-1);
    return s / 100000;
}
void generate_vector(bool v[], const int d)
{
    initRandomGenerator(g_randomGenerator);
    for (int i = 1; i <= 20 * d; i++)
    {
        v[i] = g_randomGenerator() % 2;
    }
}
void generate_vector_Schwefel(bool v[], const int d)
{
    initRandomGenerator(g_randomGenerator);
    for (int i = 1; i <= 27 * d; i++)
    {
        v[i] = g_randomGenerator() % 2;
    }
}
void generate_vector_Michalewicz(bool v[], const int d)
{
    initRandomGenerator(g_randomGenerator);
    for (int i = 1; i <= 20 * d; i++)
    {
        v[i] = g_randomGenerator() % 2;
        if (i % 20 == 0)
        {
            v[i] = 0;
        }
    }
}
double De_Jong(const bool v[], const int d)
{
    int p = 1, i;
    double s = 0;
    double val;
    for (i = 1; i <= 20 * d; i++)
    {
        if (p == 21)
        {
            val = get_value(v, 20 * (i / 20 - 1) + 1, 20 * (i / 20));
            s = s + val * val;
            p = 1;
        }
        p++;
    }
    val = get_value(v, 20 * (i / 20 - 1) + 1, 20 * (i / 20));
    s = s + val * val;
    return s;
}
double Schwefel(const bool v[], const int d)
{
    int p = 1, i;
    double s = 0;
    double val;
    for (i = 1; i <= 27 * d; i++)
    {
        if (p == 28)
        {
            val = get_value(v, 27 * (i / 27 - 1) + 1, 27 * (i / 27));
            s = s - val * sin(sqrt(abs(val)));
            p = 1;
        }
        p++;
    }
    val = get_value(v, 27 * (i / 27 - 1) + 1, 27 * (i / 27));
    s = s - val * sin(sqrt(abs(val)));
    return s;
}
double Rastrigin(const bool v[], const int d)
{
    int p = 1, i;
    double s = 10 * d;
    double val;
    for (i = 1; i <= 20 * d; i++)
    {
        if (p == 21)
        {
            val = get_value(v, 20 * (i / 20 - 1) + 1, 20 * (i / 20));
            s = s + (val * val - 10 * cos(2 * 3.14 * val));
            p = 1;
        }
        p++;
    }
    val = get_value(v, 20 * (i / 20 - 1) + 1, 20 * (i / 20));
    s = s + (val * val - 10 * cos(2 * 3.14 * val));
    return s;
}
double Michalewicz(const bool v[], const int d)
{
    int p = 1, i;
    double s = 0;
    double val;
    for (i = 1; i <= 20 * d; i++)
    {
        if (p == 21)
        {
            val = get_value(v, 20 * (i / 20 - 1) + 1, 20 * (i / 20));
            s = s + sin(val) * pow((sin(i / 20 * val * val / 3.14)), (20));
            p = 1;
        }
        p++;
    }
    val = get_value(v, 20 * (i / 20 - 1) + 1, 20 * (i / 20));
    s = s + sin(val) * pow((sin(i / 20 * val * val / 3.14)), (20));
    return (-1) * s;
}
double get_function_value(const bool v[], const int f, const int d)
{
    if (f == 1)
    {
        return De_Jong(v, d);
    }
    else
        if (f == 2)
        {
            return Schwefel(v, d);
        }
        else
            if (f == 3)
            {
                return Rastrigin(v, d);
            }
            else
            {
                return Michalewicz(v, d);
            }
}
void get_elements(const bool v[], bool aux[], const int p)
{
    int j = 1;
    for (int i = (p - 1) * 20 + 1; i <= p * 20; i++, j++)
        aux[j] = v[i];
}
void get_elements_Schwefel(const bool v[], bool aux[], const int p)
{
    int j = 1;
    for (int i = (p - 1) * 27 + 1; i <= p * 27; i++, j++)
        aux[j] = v[i];
}
double get_function_value_Hill(const bool v[], const int f, const int p, const int d)
{
    double value;
    value = get_value(v, 1, 20);
    if (f == 1)
    {
        return value * value;
    }
    else
        if (f == 2)
        {
            return value * sin(sqrt(abs(value)));
        }
        else
            if (f == 3)
            {
                return value * value - 10 * cos(2 * 3.14 * value);
            }
            else
            {
                return (-1) * sin(value) * pow((sin(p * value * value / 3.14)), (20));
            }
}
double get_function_value_Hill_Schwefel(const bool v[], const int f, const int p, const int d)
{
    double value;
    value = get_value(v, 1, 27);
    if (f == 1)
    {
        return value * value;
    }
    else
        if (f == 2)
        {
            return -value * sin(sqrt(abs(value)));
        }
        else
            if (f == 3)
            {
                return value * value - 10 * cos(2 * 3.14 * value);
            }
            else
            {
                return (-1) * sin(value) * pow((sin(p * value * value / 3.14)), (20));
            }
}
void change(bool v[], bool v_vecin[], const int d)
{
    for (int i = 1; i <= d * 20; i++)
    {
        v_vecin[i] = v[i];
    }
}
void change_Schwefel(bool v[], bool v_vecin[], const int d)
{
    for (int i = 1; i <= d * 27; i++)
    {
        v_vecin[i] = v[i];
    }
}
void change_vector(bool v[], const bool aux[], const int p)
{
    int j = 1;
    for (int i = (p - 1) * 20 + 1; i <= p * 20; i++, j++)
        v[i] = aux[j];
}
void change_vector_Schwefel(bool v[], const bool aux[], const int p)
{
    int j = 1;
    for (int i = (p - 1) * 27 + 1; i <= p * 27; i++, j++)
        v[i] = aux[j];
}
void Improve(bool aux_vecin[], const int f, const int p, const int d)
{
    double best = get_function_value_Hill(aux_vecin, f, p, d);
    double value;
    bool vector[22];
    change(aux_vecin, vector, 1);
    for (int i = 1; i <= 20; i++)
    {
        vector[i] = abs(vector[i] - 1);
        value = get_function_value_Hill(vector, f, p, d);
        if (best > value)
        {
            change(vector, aux_vecin, 1);
            best = value;
        }
        vector[i] = abs(vector[i] - 1);
    }
}
void First_Improve(bool aux_vecin[], const int f, const int p, const int d)
{
    double best = get_function_value_Hill(aux_vecin, f, p, d);
    double value;
    bool vector[22];
    change(aux_vecin, vector, 1);
    for (int i = 1; i <= 20; i++)
    {
        vector[i] = abs(vector[i] - 1);
        value = get_function_value_Hill(vector, f, p, d);
        if (best > value)
        {
            change(vector, aux_vecin, 1);
            break;
        }
        vector[i] = abs(vector[i] - 1);
    }
}
void Worst_Improve(bool aux_vecin[], const int f, const int p, const int d)
{
    double best = get_function_value_Hill(aux_vecin, f, p, d);
    double value;
    bool vector[22];
    for (int i = 1; i <= 20; i++)
        vector[i] = g_randomGenerator() % 2;
    value = get_function_value_Hill(vector, f, p, d);
    if (best > value)
        change(vector, aux_vecin, 1);
}
void Improve_Schwefel(bool aux_vecin[], const int f, const int p, const int d)
{
    double best = get_function_value_Hill_Schwefel(aux_vecin, f, p, d);
    double value;
    bool vector[30];
    change_Schwefel(aux_vecin, vector, 1);
    for (int i = 1; i <= 27; i++)
    {
        vector[i] = abs(vector[i] - 1);
        value = get_function_value_Hill_Schwefel(vector, f, p, d);
        if (best > value)
        {
            change_Schwefel(vector, aux_vecin, 1);
            best = value;
        }
        vector[i] = abs(vector[i] - 1);
    }
}
void First_Improve_Schwefel(bool aux_vecin[], const int f, const int p, const int d)
{
    double best = get_function_value_Hill_Schwefel(aux_vecin, f, p, d);
    double value;
    bool vector[30];
    change_Schwefel(aux_vecin, vector, 1);
    for (int i = 1; i <= 27; i++)
    {
        vector[i] = abs(vector[i] - 1);
        value = get_function_value_Hill_Schwefel(vector, f, p, d);
        if (best > value)
        {
            change_Schwefel(vector, aux_vecin, 1);
            break;
        }
        vector[i] = abs(vector[i] - 1);
    }
}
void Worst_Improve_Schwefel(bool aux_vecin[], const int f, const int p, const int d)
{
    double best = get_function_value_Hill_Schwefel(aux_vecin, f, p, d);
    double value;
    bool vector[30];
    for (int i = 1; i <= 27; i++)
        vector[i] = g_randomGenerator() % 2;
    value = get_function_value_Hill_Schwefel(vector, f, p, d);
    if (best > value)
        change_Schwefel(vector, aux_vecin, 1);
}
double get_random_neighbor_value(bool v[], bool v_vecin[], const int f, const int d)
{
    initRandomGenerator(g_randomGenerator);
    double val;
    change(v, v_vecin, d);
    int i = g_randomGenerator() % (d * 20) + 1;
    v_vecin[i] = abs(v_vecin[i] - 1);
    val = get_function_value(v_vecin, f, d);
    return val;
}
double get_random_neighbor_value_Schwefel(bool v[], bool v_vecin[], const int f, const int d)
{
    initRandomGenerator(g_randomGenerator);
    double val;
    change_Schwefel(v, v_vecin, d);
    int i = g_randomGenerator() % (d * 27) + 1;
    v_vecin[i] = abs(v_vecin[i] - 1);
    val = get_function_value(v_vecin, f, d);
    return val;
}
double get_random()
{
    initRandomGenerator(g_randomGenerator);
    double r = g_randomGenerator() % 10;
    return r / 10;
}

int main()
{
    int t = 0, i, p, f, d, t_final, alg, l, media = 0, imp, minim = 999999, maxim = -999999;
    bool local;
    double T = 100.00;
    double best = 99999;
    double val_v;
    double val_vecin;
    double val_number;
    cout << "Introdu 1 pentru a folosi algoritmul Hill Climbing si 2 pentru a folosi algoritmul Simulated Annealing:";
    cin >> alg;
    if (alg == 1)
    {
        cout << '\n' << "Introdu 1 pentru best improvement, 2 pentru first improvement si 3 pentru worst improvement:";
        cin >> imp;
    }
    cout << '\n' << "Numarul de dimensiuni este:";
    cin >> d;
    cout << '\n' << "Introdu 1 pentru a genera functia De Jong, 2 pentru a genera functia Schwefel, 3 pentru a genera functia Rastrigin si 4 pentru a genera functia Michalewicz:";
    cin >> f;
    auto start = std::chrono::high_resolution_clock::now();
    t_final = d * 10;
    if (alg == 1)
    {
        if (f == 1 || f == 3 || f == 4)
        {
            for (int j = 1; j <= 30; j++) {
                best = 99999;
                t = 0;
                if (f == 4)
                    generate_vector_Michalewicz(v, d);
                else
                    generate_vector(v, d);
                do
                {
                    local = false;
                    p = g_randomGenerator() % d + 1;
                    get_elements(v, aux, p);
                    val_number = get_function_value_Hill(aux, f, p, d);
                    do
                    {
                        change(aux, aux_vecin, 1);
                        if (imp == 1)
                            Improve(aux_vecin, f, p, d);
                        else
                            if (imp == 2)
                                First_Improve(aux_vecin, f, p, d);
                            else
                                Worst_Improve(aux_vecin, f, p, d);
                        val_vecin = get_function_value_Hill(aux_vecin, f, p, d);
                        if (val_vecin < val_number)
                        {
                            val_number = val_vecin;
                            change(aux_vecin, aux, 1);
                        }
                        else
                        {
                            local = true;
                            change_vector(v, aux, p);
                        }
                    } while (local != true);
                    t++;
                    val_v = get_function_value(v, f, d);
                    if (val_v < best)
                        best = val_v;
                } while (t != t_final);
                if (best < minim)
                    minim = best;
                if (best > maxim)
                    maxim = best;
                media += best;
                valori[++kk] = best;
            }
            cout << "Minimul este: " << minim;
            cout << '\n' << "Maximul este: " << maxim;
            cout << '\n' << "Media este: " << media / 30;
            int ss = 0;
            for (int kkk = 1; kkk <= kk; kkk++)
            {
                ss = ss + (valori[kkk] - media / 30) * (valori[kkk] - media / 30);
            }
            ss = ss / 30;
            cout << '\n' << "Deviatia standard: " << sqrt(ss);
        }
        else if (f == 2)
        {
            for (int j = 1; j <= 30; j++) {
                best = 99999;
                t = 0;
                generate_vector_Schwefel(v, d);
                do
                {
                    local = false;
                    p = g_randomGenerator() % d + 1;
                    get_elements_Schwefel(v, aux, p);
                    val_number = get_function_value_Hill_Schwefel(aux, f, p, d);
                    do
                    {
                        change_Schwefel(aux, aux_vecin, 1);
                        if (imp == 1)
                            Improve_Schwefel(aux_vecin, f, p, d);
                        else
                            if (imp == 2)
                                First_Improve_Schwefel(aux_vecin, f, p, d);
                            else
                                Worst_Improve_Schwefel(aux_vecin, f, p, d);
                        val_vecin = get_function_value_Hill_Schwefel(aux_vecin, f, p, d);
                        if (val_vecin < val_number)
                        {
                            val_number = val_vecin;
                            change_Schwefel(aux_vecin, aux, 1);
                        }
                        else
                        {
                            local = true;
                            change_vector_Schwefel(v, aux, p);
                        }
                    } while (local != true);
                    t++;
                    val_v = get_function_value(v, f, d);
                    if (val_v < best)
                        best = val_v;
                } while (t != t_final);
                media += best;
                valori[++kk] = best;
                if (best < minim)
                    minim = best;
                if (best > maxim)
                    maxim = best;
            }
            cout << "Minimul este: " << minim;
            cout << '\n' << "Maximul este: " << maxim;
            cout << '\n' << "Media este: " << media / 30;
            int ss = 0;
            for (int kkk = 1; kkk <= kk; kkk++)
            {
                ss = ss + (valori[kkk] - media / 30) * (valori[kkk] - media / 30);
            }
            ss = ss / 30;
            cout << '\n' << "Deviatia standard: " << sqrt(ss);
        }
    }
    else if (alg == 2)
    {
        if (f == 1 || f == 3 || f == 4)
        {
            for (int j = 1; j <= 30; j++) {
                best = 99999;
                t = 0;
                if (f == 4)
                    generate_vector_Michalewicz(v, d);
                else
                    generate_vector(v, d);
                do
                {
                    local = false;
                    p = g_randomGenerator() % d + 1;
                    get_elements(v, aux, p);
                    val_number = get_function_value_Hill(aux, f, p, d);
                    do
                    {
                        change(aux, aux_vecin, 1);
                        Improve(aux_vecin, f, p, d);
                        val_vecin = get_function_value_Hill(aux_vecin, f, p, d);
                        if (val_vecin < val_number)
                        {
                            val_number = val_vecin;
                            change(aux_vecin, aux, 1);
                        }
                        else if (get_random() < exp((-1) * abs(val_vecin - val_number) / T))
                        {
                            val_number = val_vecin;
                            change_vector(v, aux, p);
                            local = true;
                        }
                    } while (local != true);
                    T = T * 0.9;
                    t++;
                    val_v = get_function_value(v, f, d);
                    if (val_v < best)
                        best = val_v;
                } while (t < t_final);
                media += best;
                valori[++kk] = best;
                if (best < minim)
                    minim = best;
                if (best > maxim)
                    maxim = best;
            }
            cout << "Minimul este: " << minim;
            cout << '\n' << "Maximul este: " << maxim;
            cout << '\n' << "Media este: " << media / 30;
            int ss = 0;
            for (int kkk = 1; kkk <= kk; kkk++)
            {
                ss = ss + (valori[kkk] - media / 30) * (valori[kkk] - media / 30);
            }
            ss = ss / 30;
            cout << '\n' << "Deviatia standard: " << sqrt(ss);
        }
        else
        {
            for (int j = 1; j <= 30; j++) {
                best = 99999;
                t = 0;
                generate_vector_Schwefel(v, d);
                do
                {
                    local = false;
                    p = g_randomGenerator() % d + 1;
                    get_elements_Schwefel(v, aux, p);
                    val_number = get_function_value_Hill_Schwefel(aux, f, p, d);
                    do
                    {
                        change_Schwefel(aux, aux_vecin, 1);
                        Improve_Schwefel(aux_vecin, f, p, d);
                        val_vecin = get_function_value_Hill_Schwefel(aux_vecin, f, p, d);
                        if (val_vecin < val_number)
                        {
                            val_number = val_vecin;
                            change_Schwefel(aux_vecin, aux, 1);
                        }
                        else if (get_random() < exp((-1) * abs(val_vecin - val_number) / T))
                        {
                            local = true;
                            change_vector_Schwefel(v, aux, p);
                        }
                    } while (local != true);
                    T = T * 0.9;
                    t++;
                    val_v = get_function_value(v, f, d);
                    if (val_v < best)
                        best = val_v;
                } while (t < t_final);
                media += best;
                valori[++kk] = best;
                if (best < minim)
                    minim = best;
                if (best > maxim)
                    maxim = best;
            }
            cout << "Minimul este: " << minim;
            cout << '\n' << "Maximul este: " << maxim;
            cout << '\n' << "Media este: " << media / 30;
            int ss = 0;
            for (int kkk = 1; kkk <= kk; kkk++)
            {
                ss = ss + (valori[kkk] - media / 30) * (valori[kkk] - media / 30);
            }
            ss = ss / 30;
            cout << '\n' << "Deviatia standard: " << sqrt(ss);
        }
    }

    cout << endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Timp de executie: " << duration.count() << " secunde." << std::endl;///Liniile 318,514,515,516 sunt preluate de pe chatgpt pentru a calcula timpul de cand a inceput programul si pana cand se termina
    return 0;
}