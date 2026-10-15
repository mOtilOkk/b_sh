#pragma once
//#pragma warning(disable : 4996)
#include "Header.h"


class S_box
{
private:
	uint8_t n;
	std::vector<uint8_t> Pi;
public:
	S_box(uint8_t n_, std::vector<uint8_t> Pi_)
	{
		n = n_;
		Pi.resize(1 << n);
		for (uint8_t i = 0; i < (1 << n); i++)
			Pi[i] = Pi_[i];
	}

	S_box(uint8_t n_)
	{
		n = n_;
		Pi.resize(1 << n);
		for (uint8_t i = 0; i < (1 << n); i++)
			Pi[i] = i;
	}

	S_box(const std::string& filename)
	{
		std::ifstream file(filename);
		file >> std::hex;
		unsigned tmp;
		file >> tmp;
		n = tmp;
		Pi.resize(1 << n);
		for (unsigned i = 0; i < (1u << n); i++)
		{
			file >> std::hex >> tmp;
			Pi[i] = tmp;
		}
	}

	uint8_t Pi_x(uint8_t x) { return Pi[x]; }

	S_box rev_Pi()
	{
		std::vector<uint8_t> rev_Pi(1 << n);
		for (uint8_t i = 0; i < (1 << n); i++)
			rev_Pi[Pi[i]] = i;
		S_box rev(n, rev_Pi);
		return rev;
	}


	void print()
	{
		std::cout << '(';
		for (uint8_t i = 0; i < (1 << n) - 1; i++)
			std::cout << std::hex << static_cast<unsigned int>(Pi[i]) << ", ";
		std::cout << static_cast<unsigned int>(Pi[(1 << n) - 1]) << std::dec << ")\n";
	}

	bool is_permutation()
	{
		std::vector<unsigned short> p(1 << n);
		for (unsigned i = 0; i < (1 << n); ++i)
			p[i] = 1 << n;
		for (unsigned i = 0; i < (1 << n); ++i)
			p[i] = Pi[i];
		for (unsigned i = 0; i < (1 << n); ++i)
			if (p[i] == (1 << n))
				return false;
		return true;
	}

	std::vector<bool> coord_f(uint8_t i) // возвращается i-ая координатная функция
	{
		std::vector<bool> f_i(1 << n);
		for (unsigned j = 0; j < (1 << n); j++)
			f_i[j] = (Pi[j] >> i) & 1;
		return f_i;
	}


	std::vector<unsigned short> differential_uniformity()
	{
		std::vector<unsigned short> res = { 0, 0, 0 }; // {delta, a, b}
		unsigned short d = 0;
		for (unsigned short a = 1; a < (1 << n); ++a)
			for (unsigned short b = 0; b < (1 << n); ++b)
			{
				d = 0;
				for (unsigned short x = 0; x < (1u << n); ++x)
					if ((Pi[x] ^ Pi[x ^ a]) == b)
						d++;

				if (d > res[0])
				{
					res[0] = d;
					res[1] = a;
					res[2] = b;
				}
			}
		return res;
	}

	std::vector<unsigned short> algebraic_degree()
	{
		std::vector<unsigned short> res = {0 , n}; // {d, d_min}
		std::vector<std::vector<bool>> coord_func(n);

		for (unsigned i = 0; i < n; ++i)
			coord_func[i] = coord_f(i);

		std::vector<bool> sum_j_f(1 << n);
		unsigned short d_j;
		for (unsigned j = 1; j < (1 << n); ++j) // Считаем степень каждой возможной суммы координатных функций (кроме тождественного 0)
		{
			for (unsigned s = 0; s < (1 << n); ++s) // Зануляем вектор значений суммы функций
				sum_j_f[s] = false;

			for (unsigned k = 0; k < n; ++k) // Считаем сумму
				if ((j >> k) & 1)
					for (unsigned s = 0; s < (1 << n); ++s)
						sum_j_f[s] = sum_j_f[s] ^ coord_func[k][s];
			d_j = alg_deg(sum_j_f, n);
			if (d_j > res[0])
				res[0] = d_j;
			if(d_j < res[1])
				res[1] = d_j;
		}
		return res;
	}

	std::vector<unsigned short> linear_characteristic()
	{
		std::vector<unsigned short> res = { 0, 0, 0, 0 };   // {Lambda, a, b, NL(S)}

		for (unsigned short b = 1; b < (1 << n); ++b)
		{
			
			std::vector<bool> g(1 << n); // g = b · S(x) 
			for (unsigned x = 0; x < (1 << n); ++x)
			{
				unsigned short t = b & Pi[x];
				for (unsigned i = 0; i < 4; ++i)
					t ^= t >> (1 << (3 - i));
				g[x] = t&1;
			}

			std::vector<short> G = bpf(n, g);

			for (unsigned short a = 0; a < (1 << n); ++a)
			{
				short W = (a == 0) ? (1 << n) - 2 * G[a] : -2 * G[a];

				unsigned short absW = (unsigned short)(W < 0 ? -W : W);
				if (absW > res[0])
				{
					res[0] = absW;
					res[1] = a;
					res[2] = b;
				}
			}
		}
		res[3] = (1 << (n - 1)) - res[0] / 2;
		return res;
	}
};