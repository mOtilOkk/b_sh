#pragma once
#include <iostream>
#include <fstream>
#include <vector>

//unsigned pow(unsigned x, unsigned y)
//{
//	unsigned a = x;
//	unsigned res = 1;
//	for (uint8_t i = 0; i < 32; i++)
//	{
//		if ((y >> i) & 1)
//			res *= a;
//		a *= a;
//	}
//	return res;
//}

inline std::vector<short> bpf(uint8_t n, std::vector<bool>& f) // Быстрое преобразование Фурье. n - размерность БФ
{
	std::vector<short> Fur1(1 << n);
	std::vector<short> Fur2(1 << n);
	for (unsigned i = 0; i < (1 << n); i++)
		Fur1[i] = f[i];
	for (unsigned i = 0; i < n; i++)
	{
		for (unsigned j = 0; j < 1 << (n - 1 - i); j++)
		{
			for (unsigned k = 0; k < 1 << i; k++)
			{
				Fur2[j * (1 << (i + 1)) + k] = Fur1[j * (1 << (i + 1)) + k] + Fur1[j * (1 << (i + 1)) + (1 << i) + k];
				Fur2[j * (1 << (i + 1)) + (1 << i) + k] = Fur1[j * (1 << (i + 1)) + k] - Fur1[j * (1 << (i + 1)) + (1 << i) + k];
			}
		}
		for (unsigned j = 0; j < (1 << n); j++)
		{
			Fur1[j] = Fur2[j];
		}
	}
	return Fur2;
}


inline unsigned wt(unsigned u) // Вес Хэмминга
{
	unsigned c = 0;
	while (u) { c++; u &= u - 1; }
	return c;
}


inline std::vector<bool> mobius(const std::vector<bool>& f, uint8_t n) // Преобразование Мёбиуса
{
	std::vector<bool> c = f;
	unsigned N = 1u << n;

	for (unsigned i = 0; i < n; ++i)
		for (unsigned u = 0; u < N; ++u)
			if (u & (1u << i))
				c[u] = c[u] ^ c[u ^ (1u << i)];

	return c;
}


inline unsigned short alg_deg(const std::vector<bool>& f, uint8_t n) // Алгебраическая степень булевой функции
{
	std::vector<bool> c = mobius(f, n);

	unsigned short deg = 0;
	for (unsigned u = 0; u < (1u << n); ++u)
		if (c[u] && wt(u) > deg)
			deg = wt(u);

	return deg;
}

