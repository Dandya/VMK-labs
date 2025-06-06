#include "labs.h"

static
std::vector<int>
Jegalkin(std::vector<int> vec) {
	// Метод треугольника (метод Паскаля)
	const int vec_size = vec.size();
	int counter = vec.size() / 2;
	std::vector<int> buf(vec_size);
	for(int i = 0; i < N_NUMBER; i++) {
		for(int j = 0; j < vec_size; j += 2 * counter) {
			for(int k = j; k < j + counter; k++) {
				buf[k] = vec[k];
				buf[k+counter] = vec[k] ^ vec[k+counter];
			}
		}
		vec = buf;
		counter/=2;
	}
	return vec;
}

static
std::pair<std::string, std::vector<bool>>
Jegalkin_poly(std::vector<int> jeg){
	// Построение строки с многочленом Жегалкина
	// И определение фиктивных переменных
	std::vector<bool> fict_x(N_NUMBER, true);
	bool empty = true;
	std::string poly = "";
	for(int i = 0; i < jeg.size(); i++){
		if(jeg[i] == 1){
			if(i == 0){
				poly += "1";
				empty = false;
			}
			else{
				if(empty)
					empty = false;
				else
					poly += " + ";
				int bit = 1;
				bool mult_empty = true;
				for(int j=0; j<N_NUMBER; j++){
					int buf_check_j_bit = bit & i;
					if(buf_check_j_bit != 0){
						if(mult_empty)
							mult_empty = false;
						else
							poly += "*";
						poly += "x" + std::to_string(j+1);
						fict_x[j] = false;
					}
					bit <<= 1;
				}
			}
		}
	}
	return make_pair(poly, fict_x);
}

void
Labs::Lab1(){
	std::cout << "----Лабораторная работа 1----\n";
	int two_pow_n = static_cast<int>(std::pow(2, N_NUMBER));
	if (s_vec.empty()) {
		std::cout << "----1. Генерация случайной подстановки----\n";
		s_vec = CreateVecRand(two_pow_n);
	} else {
		std::cout << "----1. Генерация случайной подстановки (пропуск)----\n";
	}
	PrintVec(s_vec, "Подстановка");

	std::cout << "----2. Вектора значений координатных функций----\n";
	std::vector<std::vector<int>> f_vec_tmp(N_NUMBER, std::vector<int>(two_pow_n));
	for(int i = 0; i < N_NUMBER; i++){
		f_vec_tmp[N_NUMBER - 1 - i] = GetFuncF(i, s_vec);
	}
	f_vec = f_vec_tmp;
	for(int i = 0; i < N_NUMBER; i++)
		PrintVec(f_vec[i], "f"+std::to_string(i+1));

	std::cout << "----3. Вес координатных функций----\n";
	for(int i = 0; i < N_NUMBER; i++){
		std::cout << "||f"+std::to_string(i+1)+"|| = "<< GetWeight(f_vec[i])<<"\n";
	}

	std::cout << "----4. Многочлен Жегалкина для координатных функций----\n";
	std::vector<std::vector<bool>> f_fict_x(N_NUMBER);
	for(int i = 0; i < N_NUMBER; i++){
		auto f_vec_J = Jegalkin(f_vec[i]);
		auto buf_pair = Jegalkin_poly(f_vec_J);

		auto f_vec_J_poly = buf_pair.first;
		f_fict_x[i] = buf_pair.second;

		PrintVec(f_vec_J, "f"+std::to_string(i+1));
		std::cout << "f"+std::to_string(i+1)+": " << buf_pair.first << "\n";
	}

	std::cout << "----5. Фиктивные переменные для координатных функций----\n";
	for(int i = 0; i < N_NUMBER; i++){
		bool no_fict = true;
		std::cout << "f"+std::to_string(i+1)+": ";
		for(int j = 0; j < N_NUMBER; j++){
				if(f_fict_x[i][j] == true){
						no_fict ? std::cout << "" : std::cout << ", ";
						std::cout << "x"+std::to_string(j+1);
						no_fict = false;
				}
		}
		no_fict ? std::cout <<"нет фиктивных переменных\n" : std::cout <<"\n";
	}
}

