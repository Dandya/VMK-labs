#include "labs.h"

std::vector<double>
Labs::CreateFurieKoef(const std::vector<int>& f_vec) {
	std::vector<double> furie(f_vec.size());
	for (int i = 0; i < f_vec.size(); ++i)
		furie[i] = static_cast<double>(f_vec[i]);
	int counter = 1;

	// используется быстрое преобразование Фурье
	for(int i = 0; i < N_NUMBER; i++){
		std::vector<double> buf(furie.size());
		for(int j = 0; j < f_vec.size(); j += 2 * counter){
			for(int k=j; k<j+counter; k++){
				buf[k] = (furie[k]+furie[k+counter])/2;
				buf[k+counter] = (furie[k]-furie[k+counter])/2;
			}
		}
		counter *= 2;
		furie = buf;
	}
	return furie;
}


std::pair<bool, int>
Labs::IsCorrelativeImmunity(const std::vector<double>& walsh) {
	std::pair<bool, int> result = std::make_pair(false, 0);

	// проверка по критерию
	for(int i = 1; i <= N_NUMBER; i++) { // порядок
		bool check_zero_walsh = true;
		for(int j = 0; j < walsh.size(); j++) { // вектор u
			if(i == GetWeight(GetBitsVec(j, N_NUMBER))){
				check_zero_walsh = walsh[j] == 0.0 ? check_zero_walsh : false;
			}
			if(check_zero_walsh == false){
				break;
			}
		}
		if(check_zero_walsh)
			result = std::make_pair(check_zero_walsh, i);
		else
			break;
	}
	return result;
}

std::pair<bool, int>
Labs::Elasticity (const std::vector<int>& f_vec, std::pair<bool, int> corr_immun,
		const std::vector<double>& walsh) {
	std::pair<bool, int> result = std::make_pair(false, 0);

	// проверка по критерию
	for(int i = 1; i <= corr_immun.second; i++) {
		bool check_zero_walsh = true;
		for(int j = 0; j < walsh.size(); j++) { // вектор u
			if(i > GetWeight(GetBitsVec(j, N_NUMBER))){
				check_zero_walsh = walsh[j] == 0.0 ? check_zero_walsh : false;
			}
			if(check_zero_walsh == false){
				break;
			}
		}
		if(check_zero_walsh)
			result = std::make_pair(check_zero_walsh, i);
		else
			break;
	}
	return result;
}

void
Labs::PrintSpectrum(const std::vector<std::vector<int>>& delta){
		std::cout<< "x1...x" << N_NUMBER << "\t";
		for(int i=0; i<N_NUMBER; i++){
				std::cout<< "W_(f" << i+1 <<")\t";
		}
		std::cout<< " *2^{-" << N_NUMBER - 1 <<"}\n";

		for(int i=0; i<delta[0].size(); i++){
				std::cout << GetBits(i, N_NUMBER) <<"\t";
				for(int j = 0; j < N_NUMBER; j++){
						std::cout << delta[j][i] <<"\t";
				}
				std::cout<< "\n";
		}
		std::cout<< "\n";
}

static
std::string
GetStringPoly(const std::vector<int>& bits, bool negative){
	bool string_empty = true;
	std::string poly = "";
	for(int i = 0; i < bits.size(); i++){
		if(bits[i] == 1){
			if(string_empty)
				string_empty = false;
			else
				poly += " + ";
			poly += "x" + std::to_string(i+1);
		}
	}

	if(negative){
		if(!string_empty)
			poly += " + ";
		poly += "1";
	}
	return poly;
}

std::vector<std::string>
Labs::BestLinearApproximations(const std::vector<int>& stat_coeff) {
	int max = *max_element(begin(stat_coeff), end(stat_coeff));
	int min = *min_element(begin(stat_coeff), end(stat_coeff));
	max = abs(min) > abs(max) ? abs(min) : abs(max);
	if(max == 0) return {""};

	std::vector<std::string> result;
	for(int i=0; i<stat_coeff.size(); i++){
		if(abs(stat_coeff[i]) == max) {
			auto poly = GetStringPoly(GetBitsVec(i, N_NUMBER), stat_coeff[i]<0);
			result.push_back(poly);
		}
	}
	return result;
}

static
bool
IsEqual(double a, double b, double pr) {
	return fabs(a-b) < pr;
}

static
bool
CheckBent(const std::vector<double>& walsh) {
	if(N_NUMBER%2 != 0) return false;
	bool result = true;
	double true_value = pow(2, -(static_cast<double>(N_NUMBER)/2));
	for(int i = 0; i < walsh.size(); i++) {
		result = true_value == fabs(walsh[i]) ? result : false;
		if(!result) return result;
	}
	return result;
}

void
Labs::Lab3() {
	std::vector<std::vector<double>> furie(N_NUMBER);
	std::vector<std::vector<double>> walsh_hadamard(N_NUMBER);
	std::vector<std::vector<int>> stat_coeff(N_NUMBER);
	int number_non_corr_immun_func = 0;

	std::cout << "----Лабораторная работа 3----\n";

	for(int i = 0; i < N_NUMBER; i++){
		furie[i] = CreateFurieKoef(f_vec[i]);
		walsh_hadamard[i].resize(furie[i].size());
		stat_coeff[i].resize(furie[i].size());

		for(int j = 0; j < furie[i].size(); j++){
			walsh_hadamard[i][j] = j == 0 ? 1 - 2*furie[i][j] : -2*furie[i][j];
			stat_coeff[i][j] = static_cast<int>(pow(2, N_NUMBER-1) * walsh_hadamard[i][j]);
		}
	}

	std::cout << "----1. Корреляционная иммунность и эластичность для координатных функций----\n";

	for(int i=0; i<N_NUMBER; i++){
		std::cout << "f"+std::to_string(i+1)+": \t";
		std::pair<bool, int> corr_immun = IsCorrelativeImmunity(walsh_hadamard[i]);
		std::pair<bool, int> elastic = Elasticity(f_vec[i], corr_immun, walsh_hadamard[i]);

		corr_immun.first ?
				std::cout << "функция корреляционно иммунна порядка " +
				std::to_string(corr_immun.second)
				: std::cout << "функция НЕ корреляционно иммунна порядка 1";
		elastic.first ?
				std::cout << "\n\tфункция эластична порядка " +
				std::to_string(corr_immun.second) << "\n"
				: std::cout << "\n\tфункция НЕ эластична\n";

		if(corr_immun.first == false){
			number_non_corr_immun_func = i;
		}
	}

	std::cout << "\n----2. Построение спектра булевых функций для координатных функций----\n";
	PrintSpectrum(stat_coeff);

	std::cout << "----3. Нахождение наилучшего линейного приближения для координатных булеых функций----\n";
	for(int i = 0; i < N_NUMBER; i++){
		std::vector<std::string> lin = BestLinearApproximations(stat_coeff[i]);
		std::cout << "\n\nf"+std::to_string(i+1)+":";
		for(int j = 0; j < lin.size(); j++)
				std::cout << "\tg_"+std::to_string(j+1)+"(x) = " << lin[j] << "\n";
	}

	std::cout << "----4. Является ли координатная функция бент-функцией----\n";
	for(int i=0; i<N_NUMBER; i++){
			std::cout << "f"+std::to_string(i+1)+":\tданная координатная функция ";
			if(!CheckBent(walsh_hadamard[i]))
				std::cout << "НЕ";
			std::cout << " является бент-функцией\n";
	}

	return;
}