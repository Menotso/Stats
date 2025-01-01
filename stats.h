#ifndef STATS_H
#define STATS_H

#include <vector>
#include <map>

// Statistic Functions

class Stats {
public:
	std::vector<double> bubble_sort(std::vector<double> lista);
	double mean(std::vector<double> lista);
	double median(std::vector<double> lista);
	double mode(std::vector<double> lista);
	double range(std::vector<double> lista);
	double variance(std::vector<double> lista, bool denominator);
	double std_deviation(std::vector<double> lista);
	double covariance(std::vector<double> list1, std::vector<double> list2, bool denominator);
	double corr_pearson(std::vector<double> list1, std::vector<double> list2);
	double corr_spearman(std::vector<double> list1, std::vector<double> list2);
	int factorial(int number);
	int permutation(int number, int repeat_n_times);
	int combination(int number, int groups);
	double mean_discrete_random_var(std::vector<double> list_values, std::vector<double> list_probability);
	double variance_discrete_random_var(std::vector<double> list_values, std::vector<double> list_probability);
	std::map<std::string, double> box_plot(std::vector<double> list);
	double bernoulli_distribution(double success_prob);
	double binomial_distribution(int n_events, int n_success, double success_prob);
};

#endif // STATS_H
