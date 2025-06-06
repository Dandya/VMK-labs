#include "labs.h"

#define BAN_DIMENSION_LIMITATION 20

static
std::pair<bool, int>
GetProhibition(const std::vector<int>& freq, const int number_current_f){
	std::vector<int> statistic(pow(2, number_current_f), 0);
	for(int i = 0; i < freq.size(); i++){
			statistic[freq[i]]++;
	}
	for(int i = 0; i < statistic.size(); i++){
			if(statistic[i]==0){
					return std::make_pair(true, i);
			}
	}
	return std::make_pair(false, 0);
}

bool
Labs::IsHighlyEquiprobable(const std::vector<int>& f_vec) {
	// 2^{number_current_f+1} = диапазон значений f
	int number_current_f = 0;
	std::vector<std::vector<int>> all_f;
	all_f.push_back(f_vec);

	// проверка, есть ли значение f с встречаемостью 0 и если да - это запрет
	std::pair<bool, int> check = GetProhibition(all_f[number_current_f],
					number_current_f+1);
	while(number_current_f <= BAN_DIMENSION_LIMITATION && check.first == false){
		std::vector<int> buf_new_f(2*all_f[number_current_f].size());
		for(int i = 0; i < all_f[number_current_f].size(); i++){
			// предыдущее значение f || значение f от последних n переменных(бит)
			buf_new_f[2*i]   = (all_f[number_current_f][i]<<1) + f_vec[(2*i) & ((int)pow(2, N_NUMBER)-1)];
			buf_new_f[2*i+1] = (all_f[number_current_f][i]<<1) + f_vec[(2*i+1) & ((int)pow(2, N_NUMBER)-1)];
		}
		all_f.push_back(buf_new_f);
		number_current_f++;
		check = GetProhibition(all_f[number_current_f], number_current_f+1);
	}
	check.first ? (std::cout << "запрет = " <<
			GetBits(check.second, number_current_f+1) << "\n") :
			std::cout << "запрет не найден\n";
	return !check.first;
}

void
Labs::Lab2() {
	std::cout << "----Лабораторная работа 2----\n";
	std::cout << "----1. Преобладание нулей над единицами для координатных функций----\n";
	for(int i = 0; i < N_NUMBER; i++) {
		std::cout << "f"+std::to_string(i+1)+": " <<
					(double)GetPredominanceZero(f_vec[i]) << "\n";
	}

	std::cout << "----2,3. Сильно равновероятность для координатных функций и их запреты----\n";
	for(int i = 0; i < N_NUMBER; i++) {
		std::cout << "f"+std::to_string(i+1) + ": \t";
		if (IsHighlyEquiprobable(f_vec[i]))
			std::cout << "\tСИЛЬНО равновероятная\n";
		else
			std::cout << "\tНЕ сильно равновероятная\n";
	}
}