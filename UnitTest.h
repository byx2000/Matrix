#include <iostream>
#include <fstream>
#include "BigInt.h"
#include "Fraction.h"
#include "Matrix.h"

using namespace std;

//矩阵加法测试
void TestAdd()
{
	ifstream fin("add.txt");

	int cnt, r, c;
	fin >> cnt >> r >> c;

	Matrix<double> a(r, c), b(r, c), res(r, c);

	cout << "begin add test..."  << endl;
	for (int i = 0; i < cnt; ++i)
	{
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				fin >> a.m[i][j];
			}
		}
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				fin >> b.m[i][j];
			}
		}
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				fin >> res.m[i][j];
			}
		}

		if ((a + b) != res)
		{
			cout << "Add test failed!" << endl;
			cout << "Failed case: " << i + 1 << endl;
			fin.close();
			return;
		}
	}

	cout << "Add test passed!" << endl;
	fin.close();
	return;
}

//矩阵减法测试
void TestSub()
{
	ifstream fin("sub.txt");

	int cnt, r, c;
	fin >> cnt >> r >> c;

	Matrix<double> a(r, c), b(r, c), res(r, c);

	cout << "begin sub test..."  << endl;
	for (int i = 0; i < cnt; ++i)
	{
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				fin >> a.m[i][j];
			}
		}
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				fin >> b.m[i][j];
			}
		}
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				fin >> res.m[i][j];
			}
		}

		if ((a - b) != res)
		{
			cout << "Sub test failed!" << endl;
			cout << "Failed case: " << i + 1 << endl;
			fin.close();
			return;
		}
	}

	cout << "Sub test passed!" << endl;
	fin.close();
	return;
}

//矩阵乘法测试
void TestMul()
{
	ifstream fin("mul.txt");

	int cnt, r1, r2, r3;
	fin >> cnt >> r1 >> r2 >> r3;

	Matrix<double> a(r1, r2), b(r2, r3), res(r1, r3);

	cout << "begin mul test..."  << endl;
	for (int i = 0; i < cnt; ++i)
	{
		for (int i = 0; i < r1; ++i)
		{
			for (int j = 0; j < r2; ++j)
			{
				fin >> a.m[i][j];
			}
		}
		for (int i = 0; i < r2; ++i)
		{
			for (int j = 0; j < r3; ++j)
			{
				fin >> b.m[i][j];
			}
		}
		for (int i = 0; i < r1; ++i)
		{
			for (int j = 0; j < r3; ++j)
			{
				fin >> res.m[i][j];
			}
		}

		if ((a * b) != res)
		{
			cout << "Mul test failed!" << endl;
			cout << "Failed case: " << i + 1 << endl;
			fin.close();
			return;
		}
	}

	cout << "Mul test passed!" << endl;
	fin.close();
	return;
}

//矩阵数乘测试
void TestSca()
{
	ifstream fin("sca.txt");

	int cnt, r, c;
	fin >> cnt >> r >> c;

	Matrix<double> a(r, c), res(r, c);

	cout << "begin sca test..."  << endl;
	for (int i = 0; i < cnt; ++i)
	{
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				fin >> a.m[i][j];
			}
		}
		double k;
		fin >> k;
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				fin >> res.m[i][j];
			}
		}

		if ((k * a) != res)
		{
			cout << "Sca test failed!" << endl;
			cout << "Failed case: " << i + 1 << endl;
			fin.close();
			return;
		}
	}

	cout << "Sca test passed!" << endl;
	fin.close();
	return;
}

//行列式测试
void TestDet()
{
	ifstream fin("det.txt");

	int cnt, r;
	fin >> cnt >> r;

	Matrix<Fraction> a(r, r);
	Fraction res;
	string s;

	cout << "begin det test..."  << endl;
	for (int i = 0; i < cnt; ++i)
	{
		cout << i << endl;
		//if (i < 24)
		//	continue;
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < r; ++j)
			{
				fin >> s;
				a.m[i][j] = Fraction(s);
			}
		}
		fin >> s;

		Fraction t = a.det();
		if (t.toString() != s)
		{
			cout << "Det test failed!" << endl;
			cout << "Failed case: " << i + 1 << endl;
			cout << a << endl;
			cout << s << endl;
			cout << t << endl;
			fin.close();
			return;
		}
	}

	cout << "Det test passed!" << endl;
	fin.close();
	return;
}

//初等变换测试
void TestRow()
{
	ifstream fin("row.txt");

	int cnt, r, c;
	fin >> cnt >> r >> c;

	Matrix<Fraction> a(r, c);
	Matrix<Fraction> res(r, c);
	string s, s1;

	cout << "begin row test..."  << endl;
	for (int i = 0; i < cnt; ++i)
	{
		cout << i << endl;
		//if (i < 24)
		//	continue;
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				fin >> s;
				//cout << s << endl;
				//cout << s1 << endl;
				a.m[i][j] = Fraction(s);
			}
		}
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				fin >> s >> s1;
				res.m[i][j] = Fraction(s, s1);
			}
		}

		if (a.leastLine() != res)
		{
			cout << "Row test failed!" << endl;
			cout << "Failed case: " << i + 1 << endl;
			cout << a << endl;
			cout << res << endl;
			cout << a.leastLine() << endl;
			//cout << a << endl;
			//cout << s << endl;
			//cout << t << endl;
			fin.close();
			return;
		}
	}

	cout << "Row test passed!" << endl;
	fin.close();
	return;
}

//逆矩阵测试
void TestInv()
{
	ifstream fin("inv.txt");

	int cnt, r;
	fin >> cnt >> r;

	Matrix<Fraction> a(r, r);
	Matrix<Fraction> res(r, r);
	string s, s1;

	cout << "begin inv test..."  << endl;
	for (int i = 0; i < cnt; ++i)
	{
		cout << i << endl;
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < r; ++j)
			{
				fin >> s;
				a.m[i][j] = Fraction(s);
			}
		}
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < r; ++j)
			{
				fin >> s >> s1;
				res.m[i][j] = Fraction(s, s1);
			}
		}

		Matrix<Fraction> X = a.inverse();
		//a.inverse(X);
		if (X != res)
		{
			cout << "Inv test failed!" << endl;
			cout << "Failed case: " << i + 1 << endl;
			cout << a << endl;
			cout << res << endl;
			cout << X << endl;
			//cout << X.m[2][1].a << " " << X.m[2][1].b << " " << X.m[2][1].sign << endl;
			//cout << res.m[2][1].a << " " << res.m[2][1].b << " " << res.m[2][1].sign << endl;
			fin.close();
			return;
		}
	}

	cout << "Inv test passed!" << endl;
	fin.close();
	return;
}

void TestInv2()
{
	ifstream fin("inv2.txt");

	int cnt, r;
	fin >> cnt >> r;

	Matrix<double> a(r, r);
	//Matrix<Fraction> res(r, r);
	double x;

	for (int i = 0; i < cnt; ++i)
	{
		cout << i << endl;
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < r; ++j)
			{
				fin >> a.m[i][j];
				//cout << a.m[i][j] << endl;
			}
		}
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < r; ++j)
			{
				fin >> x;
				//res.m[i][j] = Fraction(s, s1);
			}
		}

		//cout << a << endl;
		Matrix<double> X = a.inverse();
		cout << X << endl;
		//cout << a.det() << endl;
	}

	cout << "Inv test passed!" << endl;
	fin.close();
	return;
}

//求幂测试
void TestPow()
{
    ifstream fin("pow.txt");

	int cnt, r;
	fin >> cnt >> r;
	//cout << cnt << r << endl;
	Matrix<BigInt> a(r, r), b(r, r);
	string s;
	long long n;

	cout << "begin pow test..."  << endl;
	for (int i = 0; i < cnt; ++i)
    {
        cout << i << endl;
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < r; ++j)
			{
			    fin >> s;
				a.m[i][j] = BigInt(s);
				//cout << a.m[i][j] << endl;
			}
		}
		fin >> n;
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < r; ++j)
			{
			    fin >> s;
				b.m[i][j] = BigInt(s);
			}
		}

		//cout << a << endl;
		//cout << b << endl;

		if ((a ^ n) != b)
        {
            cout << "Pow test failed!" << endl;
			cout << "Failed case: " << i + 1 << endl;
			fin.close();
			return;
        }
    }

    cout << "Pow test passed!" << endl;
	fin.close();
}
