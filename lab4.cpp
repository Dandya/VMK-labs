#include "labs.h"
#include "SNDT.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
// #include <random>
#include <algorithm>
#include <omp.h>

#define MAX_GAMMA_SIZE pow(2, 32)

std::vector<int> func_f_non_corr_immun;

// Веса вектора бит
static
int
GetWeightVecBit(bool* vec, unsigned long long int vec_size) {
	int counter = 0;
	for (unsigned long long i = 0; i < vec_size; i++) {
		if (vec[i])
			counter++;
	}
	return counter;
}

// Подсчёт коэффициентов статической структуры для функции
std::vector<int>
Labs::CreateStatStructCoeff(const std::vector<int>& f_vec) {
	std::vector<double> Furie = CreateFurieKoef(f_vec);
	size_t sz = Furie.size();
	std::vector<double> Walsh_Hadamard(sz);
	std::vector<int> Stat_Struct_Coeff(sz);

	for (size_t j = 0; j < sz; j++) {
		Walsh_Hadamard[j] = (j == 0 ? 1 - 2 * Furie[j] : -2 * Furie[j]);
		Stat_Struct_Coeff[j] = (int)(pow(2, N_NUMBER - 1) * Walsh_Hadamard[j]);
	}
	return Stat_Struct_Coeff;
}

// Создание первоначальных регистров
std::vector<std::vector<int>>
Labs::CreateReg(const std::vector<int>& register_lengths) {
	std::vector<std::vector<int>> Registers(register_lengths.size());
	int buf;
	srand((unsigned)time(nullptr));
	for (size_t i = 0; i < register_lengths.size(); i++) {
			int len = register_lengths[i];
			int max_val = (int)pow(2, len) - 1;
			do {
					buf = rand() % max_val;
			} while (buf < 1);
			Registers[i] = GetBitsVec(buf, len);
	}
	return Registers;
}

// Создание истинного ключа из регистров
static
std::vector<int>
CreateTrueKey(const std::vector<std::vector<int>>& Registers,
		const std::vector<int>& register_lengths) {
	std::vector<int> key;
	for (size_t i = 0; i < register_lengths.size(); i++) {
		for (int j = register_lengths[i] - 1; j >= 0; j--) {
			key.push_back(Registers[i][j]);
		}
	}
	return key;
}

// Генерация нового значения внутри регистра
static
int
FuncReg(const std::vector<int>& vec, int reg_size) {
	int start = vec.size() - reg_size;
	int result = 0;
	switch (reg_size) {
		case 5:  // Для регистра 1: x^5 + x^2 + 1 (примитивный)
			result = vec[start] ^ vec[start+3];
			break;
		case 7:  // Для регистра 2: x^7 + x^6 + 1 (примитивный)
			result = vec[start] ^ vec[start+1];
			break;
		case 9:  // Для регистра 3: x^9 + x^4 + 1 (примитивный)
			result = vec[start] ^ vec[start+5];
			break;
		case 11: // Для регистра 4: x^11 + x^2 + 1 (примитивный)
			result = vec[start] ^ vec[start+9];
			break;
		case 13: // Для регистра 5: x^13 + x^4 + x^3 + x + 1 (неприводимый)
			result = vec[start] ^ vec[start+9] ^ vec[start+10] ^ vec[start+12];
			break;
		case 19: // Для регистра 6: x^19 + x^5 + x^2 + x + 1 (неприводимый)
			result = vec[start] ^ vec[start+14] ^ vec[start+17] ^ vec[start+18];
			break;
		default: // Стандартное правило для других длин
			result = vec[start] ^ vec[start+1];
	}
	return result;
}

// Получение индекса числа из бита-вектора: перевод на целочисленный номер
static
int
VecToInt(const std::vector<int>& x) {
	int number = 0;
	for (size_t i = 0; i < x.size(); i++) {
		number = (number << 1) | (x[i] ? 1 : 0);
	}
	return number;
}

// Генерация нового элемента гаммы
static
int
GenGammaF(const std::vector<int>& x) {
	int number = VecToInt(x);
	return func_f_non_corr_immun[number];
}

// Функция подсчёта НЕ совпадений в векторах-значениях xi и f в таблице истинности
static std::vector<double> TruthTableF(int N) {
	std::vector<double> xn(N, 0.0);
	int size = (int)pow(2, N);
	// Параллелим по i, но аккуратно собираем суммы xn
	std::vector<std::vector<int>> local_counts;
	int num_threads = omp_get_max_threads();
	local_counts.assign(num_threads, std::vector<int>(N, 0));

	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		std::vector<int> i_vec(N);
		for (int i = tid; i < size; i += num_threads) {
			// Преобразуем i в бит-вектор длины N
			int tmp = i;
			for (int bit = N - 1; bit >= 0; bit--) {
				i_vec[bit] = tmp & 1;
				tmp >>= 1;
			}
			int f = GenGammaF(i_vec);
			for (int j = 0; j < N; j++) {
				if (i_vec[j] == f) {
					local_counts[tid][j]++;
				}
			}
		}
	}

	// Суммируем из каждого потока
	for (int t = 0; t < num_threads; t++) {
		for (int j = 0; j < N; j++) {
			xn[j] += local_counts[t][j];
		}
	}

	for (int i = 0; i < N; i++) {
		xn[i] = 1.0 - xn[i] / pow(2.0, N);
	}
	return xn;
}

// Генерируем не корреляционно-иммунную функцию длины 2^N_NUMBER
std::vector<int>
Labs::CreateSpecificVecRand() {
	printf("Подбираем функцию f, удовлетворяющую условиям...\n");
	bool check_corr_immun;
	bool check_stat_struct_coeff;
	bool check_p;
	std::vector<std::vector<int>> result(N_NUMBER);
	do {
		check_corr_immun = false;
		check_stat_struct_coeff = false;
		check_p = false;
		std::vector<int> s_vec = CreateVecRand((int)pow(2, N_NUMBER));
		for (int i = 0; i < N_NUMBER && !check_corr_immun; i++) {
			result[i] = GetFuncF(i, s_vec);
			std::vector<double> Furie = CreateFurieKoef(result[i]);
			size_t sz = Furie.size();
			std::vector<double> Walsh_Hadamard(sz);
			std::vector<int> Stat_Struct_Coeff(sz);

			for (size_t j = 0; j < sz; j++) {
				Walsh_Hadamard[j] = (j == 0 ? 1 - 2 * Furie[j] : -2 * Furie[j]);
				Stat_Struct_Coeff[j] = (int)(pow(2, N_NUMBER - 1) * Walsh_Hadamard[j]);
			}
			std::pair<bool, int> corr_immun = IsCorrelativeImmunity(Walsh_Hadamard);
			check_corr_immun = corr_immun.first;
			if (!check_corr_immun) {
				check_stat_struct_coeff = true;
				for (int k = 0; k < N_NUMBER; k++) {
					int idx = (0b1 << (N_NUMBER - 1 - k));
					if (Stat_Struct_Coeff[idx] <= 0) {
						check_stat_struct_coeff = false;
						break;
					}
				}
				if (check_stat_struct_coeff) {
					func_f_non_corr_immun = result[i];
					std::vector<double> p = TruthTableF(N_NUMBER);
					check_p = true;
					for (int k = 0; k < N_NUMBER; k++) {
						if (p[k] == 0.5) {
							check_p = false;
							break;
						}
					}
					if (check_p) {
						return result[i];
					}
				}
			}
		}
	} while (check_corr_immun || !check_stat_struct_coeff || !check_p);
	return {}; // Должно никогда сюда не дойти
}

// Функция, создающая гамму длины gamma_size на основе регистров
static
std::vector<int>
CreateGamma(int gamma_size, int reg_size, const std::vector<std::vector<int>>& Registers_const, const std::vector<int>& register_lengths) {
	std::vector<std::vector<int>> Registers = Registers_const; // локальная копия
	std::vector<int> gamma(gamma_size);
	std::vector<int> current_vec_x(reg_size);

	for (int i = 0; i < gamma_size; i++) {
		for (int j = 0; j < reg_size; j++) {
			current_vec_x[j] = Registers[j][i];
			if ((int)Registers[j].size() < gamma_size) {
				std::vector<int> buf_reg(register_lengths[j]);
				copy_n(Registers[j].end() - register_lengths[j], register_lengths[j], buf_reg.begin());
				Registers[j].push_back(FuncReg(buf_reg, register_lengths[j]));
			}
		}
		gamma[i] = GenGammaF(current_vec_x);
	}
	return gamma;
}

// Поиск кандидатов в регистры с учетом вероятностей p, q
std::vector<std::vector<int>>
Labs::SearchPartKey(double p, double q, int candidate_len, int number_candidate, const std::vector<std::vector<int>>& Registers, const std::vector<int>& register_lengths) {
	Standard_Normal_Distribution_Table SNDT = Standard_Normal_Distribution_Table();
	double alpha = 0.0, beta = 0.0;
	double T = 0.0, C = 0.0;
	double Qantil_n_alpha = 0.0, Qantil_n_beta = 0.0;

	std::cout << "\n\tПодбираем кандидатов в регистр " << number_candidate << "\n";

	// Подбор порогов alpha, beta
	do {
			alpha += 0.0000001;
			beta += 0.0000001;
			Qantil_n_alpha = SNDT.FindX(1 - alpha);
			Qantil_n_beta = SNDT.FindX(1 - beta);
			T = pow((Qantil_n_alpha * sqrt(p * (1 - p)) + Qantil_n_beta * sqrt(q * (1 - q))), 2) / pow(q - p, 2);
	} while (Qantil_n_alpha == ERROR_CODE || Qantil_n_beta == ERROR_CODE);

	std::cout << "\tT = " << T << "\n";
	double min_C = p < q ? T * p : T * q;
	double max_C = p > q ? T * p : T * q;
	C = Qantil_n_alpha * sqrt(T * p * (1 - p)) + T * p;
	if (C < min_C || C > max_C) {
			// C = (-Qantil_n_beta) * sqrt(T * p * (1 - q)) + T * q; // при необходимости
	}
	std::cout << "\tC = " << C << "\n";

	std::vector<std::vector<int>> result; // найденные кандидаты
	std::vector<int> mismatches_vec; // число несовпадений для каждого найденного кандидата

	int max_candidate = (int)pow(2, candidate_len);
	int reg_size = (int)register_lengths.size();

	// Генерация гаммы один раз (последняя длина регистров — передаётся в функцию)
	std::vector<int> gamma = CreateGamma((int)T, reg_size, Registers, register_lengths);

	int positive_count = 0;
	int min_mismatches = (int)T; // инициализируем максимально возможным
	// Параллельный цикл по всем потенциальным базовым значениям регистра
	#pragma omp parallel
	{
			std::vector<int> buf_res_local;            // локальный буфер для буферизации кандидата
			std::vector<int> buf_candidate_local(candidate_len);
			for (int i = 1; i < max_candidate; i++) {
					// Разделение по потокам: каждый поток обрабатывает свой диапазон i
					int thread_id = omp_get_thread_num();
					int num_threads = omp_get_num_threads();
					if (i % num_threads != thread_id) continue;

					// Инициализируем начальный вектор длины candidate_len
					buf_res_local = GetBitsVec(i, candidate_len);

					// Достраиваем до длины T
					for (int j = candidate_len; j < (int)T; j++) {
							copy_n(buf_res_local.end() - candidate_len, candidate_len, buf_candidate_local.begin());
							buf_res_local.push_back(FuncReg(buf_candidate_local, candidate_len));
					}

					// Считаем несовпадения с гаммой
					int mismatches = 0;
					for (int j = 0; j < (int)T; j++) {
							if (gamma[j] != buf_res_local[j]) {
									mismatches++;
									if (mismatches >= C) break;
							}
					}

					if (mismatches < C) {
							// Кандидат найден, сохраняем его (в критической секции)
							#pragma omp critical
							{
									// Сохраняем только базовую часть (первые candidate_len бит)
									std::vector<int> base_candidate(candidate_len);
									for (int k = 0; k < candidate_len; k++) {
											base_candidate[k] = buf_res_local[k];
									}
									result.push_back(base_candidate);
									mismatches_vec.push_back(mismatches);
									positive_count++;
							}
					}
			}
	}

	if (positive_count > 0) {
			// Определяем индекс кандидата с минимальным количеством несовпадений
			int min_idx = 0;
			for (int i = 1; i < positive_count; i++) {
					if (mismatches_vec[i] < mismatches_vec[min_idx]) {
							min_idx = i;
					}
			}
			// Выводим всех кандидатов, помечая минимальный
			for (int i = 0; i < positive_count; i++) {
					std::string label = "\n\t\tcandidate " + std::to_string(i + 1);
					if (i == min_idx) {
							label += " (min mismatches)";
					}
					PrintVec(result[i], label);
			}
			std::cout << "\n\tУспешно!" << "\n";
	}
	else {
			std::cout << "\n\tНе удалось найти кандидатов! Перезапустите программу, пожалуйста!" << "\n";
	}

	return result;
}

// Функция отсеивания ложных кандидатов путём непосредственной проверки
std::vector<std::vector<int>>
SearchTrueKey(int count, std::vector<std::vector<int>> keys, const std::vector<std::vector<std::vector<int>>>& supposed_parts_key, std::vector<int> vec_numb, int reg_size, const std::vector<int>& register_lengths) {
	if (count == reg_size) {
		std::vector<std::vector<int>> current_key(reg_size);
		for (int i = 0; i < reg_size; i++) {
				current_key[i] = supposed_parts_key[i][vec_numb[i]];
		}
		keys.push_back(CreateTrueKey(current_key, register_lengths));
		return keys;
	}

	for (size_t i = 0; i < supposed_parts_key[count].size(); i++) {
		vec_numb[count] = (int)i;
		keys = SearchTrueKey(count + 1, keys, supposed_parts_key, vec_numb, reg_size, register_lengths);
	}
	return keys;
}

// main-функция лабораторной работы 4
void
Labs::Lab4() {
	// Генерируем не корреляционно-иммунную функцию
	func_f_non_corr_immun = CreateSpecificVecRand();

	std::vector<int> Stat_Struct_Coeff = CreateStatStructCoeff(func_f_non_corr_immun);
	PrintVec(func_f_non_corr_immun, "Функция генерации гаммы");

	std::vector<int> register_lengths = {5, 7, 9, 11, 13, 19};
	int reg_size = (int)register_lengths.size();
	if (reg_size != N_NUMBER) {
		std::cout << "Ошибка в размерности. Корректный запуск программы невозможен" << "\n";
		return;
	}

	printf("\n");
	const std::vector<std::vector<int>> Registers = CreateReg(register_lengths);
	for (int i = 0; i < reg_size; i++) {
		PrintVec(Registers[i], "Регистр " + std::to_string(i + 1));
	}

	std::vector<int> key = CreateTrueKey(Registers, register_lengths);
	PrintVec(key, "\nИстинный ключ");
	printf("\n");

	// p — частота НЕ совпадений в таблице истинности, q — частота совпадений y и f
	std::vector<double> P_f_no_equals_xn = TruthTableF(reg_size);
	std::vector<double> Q_f_no_equals_y(reg_size);

	for (int i = 0; i < reg_size; i++) {
		if (Stat_Struct_Coeff[(0b1 << (reg_size - 1 - i))] > 0) {
			Q_f_no_equals_y[i] = 0.5 + (double)abs(Stat_Struct_Coeff[(0b1 << (reg_size - 1 - i))]) / pow(2, reg_size);
		}
		else {
			Q_f_no_equals_y[i] = 1.0 - P_f_no_equals_xn[i];
		}
		printf("q%d = %lf\tp%d = %lf\n", i + 1, Q_f_no_equals_y[i], i + 1, P_f_no_equals_xn[i]);
	}

	std::vector<std::vector<std::vector<int>>> supposed_parts_key(reg_size);
	for (int i = 0; i < reg_size; i++) {
		supposed_parts_key[i] = SearchPartKey(1.0 - Q_f_no_equals_y[i], 0.5, register_lengths[i], i + 1, Registers, register_lengths);
		if (supposed_parts_key[i].empty()) {
			std::cout << "\nВеликий бог рандома не даровал Вам удачу в этот раз... Чтож, попробуйте снова!" << "\n";
			return;
		}
	}

	std::vector<std::vector<int>> keys;
	keys = SearchTrueKey(0, keys, supposed_parts_key, std::vector<int>(reg_size), reg_size, register_lengths);

	for (size_t i = 0; i < keys.size(); i++) {
			PrintVec(keys[i], "Выбранный ключ " + std::to_string(i + 1));
			std::vector<int> check_true_key(key.size());
			for (size_t j = 0; j < key.size(); j++) {
					check_true_key[j] = key[j] ^ keys[i][j];
			}
			std::cout << "Алгоритм ошибся на " << GetWeight(check_true_key) << " бит" << "\n";
	}
	return;
}
