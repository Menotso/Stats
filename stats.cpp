#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric> // accumulate = sum for vector elements
#include "stats.h"

using namespace std;

vector<double> Stats::bubble_sort(vector<double>lista) {
    if (lista.empty()) return { 0.0 };
    int tam = lista.size();
    for (int i = 0; i < tam; i++) 
        {
        for (int j = 0; j < tam; j++) 
            {
            if (lista[j] > lista[j + 1]) swap(lista[j], lista[j+1]);
            }
        }

    return lista;
}

double Stats::mean(vector<double>lista) {
    if (lista.empty()) return 0.0;
    int sum_elements=0, tam=lista.size();
    for (int i = 0; i < tam; i++) sum_elements += i;
    return sum_elements / tam;
}

double Stats::median(vector<double>lista) {
    if (lista.empty()) return 0.0;
    int tam = lista.size();
    lista = bubble_sort(lista);

    if (tam % 2 == 0) 
        {
        return (lista[floor(tam / 2)] + lista[floor(tam / 2) - 1]);
        }
    else 
        {
        return lista[(floor(tam - 1) / 2)];
        }
}

double Stats::mode(vector<double>lista) {
    if (lista.empty()) return 0.0;
    map<double, int> counter;

    for (double valor : lista) counter[valor]++;

    auto it2 = max_element(counter.begin(), counter.end(), [](const pair<double, int>& a, const pair<double, int>& b) {
        return a.second < b.second;
        });

    if (it2 != counter.end()) {
        return it2->first;
    }
}

double Stats::range(vector<double>lista) {
    if (lista.empty()) return 0.0;
    vector<double>::iterator itMin, itMax;
    
    itMin = min_element(lista.begin(), lista.end());
    itMax = max_element(lista.begin(), lista.end());

    return itMax - itMin;
}

double Stats::variance(vector<double>lista, bool denominator) {
    if (lista.empty()) return { 0.0 };
    double sum_var = 0.0;
    int tam = lista.size();

    for (int i = 0; i < tam; i++) {
        sum_var += pow((lista[i] - mean(lista)), 2);
    }

    double pure_variance = sum_var;
    double variance = pure_variance / (lista.size() - 1);

    if (!denominator) {
        return variance;
    }
    else {
        return pure_variance;
    }
}

double Stats::std_deviation(vector<double>lista) {
    return variance(lista, false);
}

double Stats::covariance(vector<double>list1, vector<double>list2, bool denominator) {
    double meanList1 = mean(list1);
    double meanList2 = mean(list2);
    int tam = list1.size();

    vector<double> vecElements;

    for (int i = 0; i < tam; i++) {
        for (int j = 0; j < tam; j++) {
            vecElements[i] = (list1[i] - meanList1) * (list2[j] - meanList2);
        }
    }

    double pure_covariance = accumulate(vecElements.begin(), vecElements.end(), 0);
    double covariance = pure_covariance / (tam - 1);

    if (!denominator) {
        return covariance;
    }
    else {
        return pure_covariance;
    }
}

double Stats::corr_pearson(vector<double>list1, vector<double>list2) {
    return (covariance(list1, list2, true) / sqrt(variance(list1, true) * variance(list2, true)));
}

double Stats::corr_spearman(vector<double>list1, vector<double>list2) {
    int tam = list1.size();
    double sum_var = 0.0;

    vector<double> vecElements;

    for (int i = 0; i < tam; i++) {
        for (int j = 0; j < tam; j++) {
            vecElements[i] = pow((list1[i] - list2[j]),2);
        }
    }
    
    for (int i = 0; i < tam; i++) sum_var += vecElements[i]; // sum of vecElements
    
    double correlation = 1 - ((6 * sum_var) / (tam * (pow(tam, 2) - 1)));

    return correlation;
}

int Stats::factorial(int number) {
    if (number == 0) {
        return 1;
    }
    return number * factorial(number - 1);
}

int Stats::permutation(int number, int repeat_n_times) {
    int denominator = 1;

    for (int i = 0; i < repeat_n_times; i++) {
        denominator *= factorial(i);
    }
    return factorial(number) / denominator;
}

int Stats::combination(int number, int groups) {
    return factorial(number) / (factorial(number-groups) * factorial(groups));
}

double Stats::mean_discrete_random_var(vector<double>list_values, vector<double>list_probability) {
    vector<double> vecElements;
    int tam = list_values.size();
    double sum_var = 0.0;

    for (int i = 0; i < tam; i++) {
        for (int j = 0; j < tam; j++) {
            vecElements[i] = (list_values[i] * list_probability[j]);
            sum_var += vecElements[i];
        }
    }

    return sum_var;
}

double Stats::variance_discrete_random_var(vector<double>list_values, vector<double>list_probability) {
    vector<double> vecElements;
    int tam = list_values.size();
    double sum_var = 0.0;

    for (int i = 0; i < tam; i++) {
        for (int j = 0; j < tam; j++) {
            vecElements[i] = (pow(list_values[i] - mean_discrete_random_var(list_values, list_probability),2) * list_probability[j]);
            sum_var += vecElements[i];
        }
    }
    return sum_var;
}

map<string, double> Stats::box_plot(vector<double>list) {
    map<string, double>boxPlotInfo;
    vector<double>lower, upper;
    double vecMedian = median(list);

    for (double value : list) {
        if (value < vecMedian) {
            lower.push_back(value);
        }
        else {
            upper.push_back(value);
        }
    }

    double quartile1 = median(lower);
    double quartile3 = median(upper);
    double iqr = quartile3 - quartile1;
    double minWhisker = quartile1 - (1.5 * iqr);
    double maxWhisker = quartile3 + (1.5 * iqr);

    boxPlotInfo["IQR"] = iqr;
    boxPlotInfo["minWhisker"] = minWhisker;
    boxPlotInfo["maxWhisker"] = maxWhisker;
    boxPlotInfo["quartile1"] = quartile1;
    boxPlotInfo["median"] = vecMedian;
    boxPlotInfo["quartile3"] = quartile3;

    return boxPlotInfo;
}

double Stats::bernoulli_distribution(double success_prob) {
    double failure_prob = 1 - success_prob;
    return failure_prob;
}

double Stats::binomial_distribution(int n_events, int n_success, double success_prob) {
    int combinationVal = combination(n_events, n_success);
    double failures = 1 - success_prob;
    double equation = (pow((success_prob), n_success)) * pow(failures, (n_events - n_success)); // check this last part

    return combinationVal * equation;
}

// TO DO: z-score function that return the z_probability based on two arguments: decimal and centesimal
