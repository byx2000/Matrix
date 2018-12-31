#pragma once

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>

using namespace std;

//double�ж����
bool DoubleEqual(double x1, double x2)
{
    return fabs(x1 - x2) < 1e-8;
}

//float�ж����
bool FloatEqual(float f1, float f2)
{
    return fabs(f1 - f2) < 1e-3;
}

template<typename T>
struct Matrix
{
	int r, c;
	vector<vector<T>> m;

	//���캯��
	Matrix()
	{
		init(1, 1, T(0));
	}

	//r��c�������
	Matrix(int r, int c)
	{
		init(r, c, T(0));
	}

	//r��c��Ԫ��ȫ����num�ľ���
	Matrix(int r, int c, T num)
	{
		init(r, c, num);
	}

	//���ݸ�����ά���鹹�����
	Matrix(const T *a, int r, int c)
	{
		init(r, c, 0);
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				m[i][j] = *(a + i * c + j);
			}
		}
	}

	/*Matrix(const T **a, int r, int c)
	{
        const T *b = (T*)a;
        init(r, c, 0);
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				m[i][j] = *(b + i * c + j);
			}
		}
	}*/

	//��ʼ��r��c�о��� Ԫ��ȫ����num
	void init(int r, int c, const T &num)
	{
		assert(r >= 1 && c >= 1);

		this->r = r;
		this->c = c;
		m.resize(r);
		for (int i = 0; i < r; ++i)
		{
			m[i].resize(c);
			for (int j = 0; j < c; ++j)
			{
				m[i][j] = num;
			}
		}
	}

	//ת��Ϊ�ַ���
	string toString() const
	{
		stringstream ss;
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				ss << m[i][j] << " ";
			}
			ss << endl;
		}
		return ss.str();
	}

	//���þ����i��j��Ԫ��Ϊval
	void setElem(int i, int j, T val)
	{
        assert(i >= 0 && i < r);
        assert(j >= 0 && j < c);
        m[i][j] = val;
	}

	//����==
	bool operator==(const Matrix<T> &mat) const
	{
		if (r != mat.r || c != mat.c)
		{
			return false;
		}

		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				if (m[i][j] != mat.m[i][j])
				{
					//cout << i << " " << j << endl;
					return false;
				}
			}
		}

		return true;
	}

	//����!=
	bool operator!=(const Matrix<T> &mat) const
	{
		return !((*this) == mat);
	}

	//����+
	Matrix<T> operator+(const Matrix<T> &mat) const
	{
		assert(r == mat.r && c == mat.c);

		Matrix<T> res(r, c);
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				res.m[i][j] = m[i][j] + mat.m[i][j];
			}
		}

		return res;
	}

	//����+=
	Matrix<T> operator+=(const Matrix<T> &mat)
	{
        return (*this) = (*this) + mat;
	}

	//����-
	Matrix<T> operator-(const Matrix<T> &mat) const
	{
		assert(r == mat.r && c == mat.c);

		Matrix<T> res(r, c);
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				res.m[i][j] = m[i][j] - mat.m[i][j];
			}
		}

		return res;
	}

	//����-=
	Matrix<T> operator-=(const Matrix<T> &mat)
	{
        return (*this) = (*this) - mat;
	}

	//����* (����˷�)
	Matrix<T> operator*(const Matrix<T> &mat) const
	{
		assert(c == mat.r);

		Matrix<T> res(r, mat.c);
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < mat.c; ++j)
			{
				for (int k = 0; k < c; ++k)
				{
					res.m[i][j] += (m[i][k] * mat.m[k][j]);
				}
			}
		}

		return res;
	}

	//����*=
	Matrix<T> operator*=(const Matrix<T> &mat)
	{
        return (*this) = (*this) * mat;
	}

	//����
	static Matrix<T> power(const Matrix<T> &A, long long n)
	{
	    if (n == 0)
        {
            return A.unit();
        }
        else if (n == 1)
        {
            return A;
        }

        if (n % 2 == 0)
        {
            return power(A, n / 2) * power(A, n / 2);
        }
        else
        {
            return A * power(A, (n - 1) / 2) * power(A, (n - 1) / 2);
        }
	}

	//����^ ����
    Matrix<T> operator^(long long n)
    {
        assert(r == c);
        assert(n >= 0);
        return power(*this, n);
    }

    //����^=
	Matrix<T> operator^=(long long n) const
	{
        return (*this) = (*this) ^ n;
	}

	//����ʽ
	Matrix complement(int row, int col) const
	{
		assert(r == c);
		assert(r >= 2 && c >= 2);
		assert(row >= 0 && row < r && col >= 0 && col < c);

		Matrix<T> res(r - 1, c - 1);

		int ti = 0, tj = 0;
		for (int i = 0; i < r; ++i)
		{
			if (i != row)
			{
				for (int j = 0; j < c; ++j)
				{
					if (j != col)
					{
						res.m[ti][tj] = m[i][j];
						tj++;
					}
				}
				ti++;
				tj = 0;
			}
		}

		return res;
	}

	//��������ʽչ����������ʽ
	T expand(const Matrix<T> &mat) const
	{
		if (mat.r == 1)
		{
			return mat.m[0][0];
		}

		T res = 0;
		for (int i = 0; i < mat.c; ++i)
		{
			if (i % 2 == 0)
			{
				res += mat.m[0][i] * expand(mat.complement(0, i));
			}
			else
			{
				res -= mat.m[0][i] * expand(mat.complement(0, i));
			}
		}

		return res;
	}

	//ת��
	Matrix<T> transpose()
	{
		Matrix<T> res(c, r);

		for (int i = 0; i < c; ++i)
		{
			for (int j = 0; j < r; ++j)
			{
				res.m[i][j] = m[j][i];
			}
		}

		return res;
	}

	//���ȱ任��������ʽ
	T det() const
	{
		assert(r == c);

		T res = T(1);

		Matrix<T> t = *this;
		for (int i = 0; i < c - 1; ++i) //����i�������Ԫ�ر�Ϊ0
		{
			if (t.m[i][i] == T(0))
			{
				int j;
				for (j = i + 1; j < r; ++j)
				{
					if (t.m[j][i] != T(0))
					{
						swap(t.m[i], t.m[j]);
						res *= T(-1);
						break;
					}
				}
				if (j == r)
				{
					return res = T(0);
				}
			}

			for (int j = i + 1; j < r; ++j)
			{
				if (t.m[j][i] != T(0))
				{
					T s = t.m[j][i] * T(-1) / t.m[i][i];
					for (int k = 0; k < c; ++k)
					{
						t.m[j][k] += (s * t.m[i][k]);
					}
				}
				//cout << t << endl;
			}
		}

		for (int i = 0; i < c; ++i)
		{
			res *= (t.m[i][i]);
		}

		//cout << t << endl;

		return res;
	}

	//���ȱ任 ������Ϊ�������
	Matrix<T> leastLine()
	{
		Matrix<T> res = *this;
		//cout << res.r << " " << res.c << endl;
		//cout << res.m[1][0] << endl;

		for (int i = 0; i < r; ++i)
		{
			int master = i; //��Ԫ����
			while (master < c && res.m[i][master] == T(0))
			{
				int j;
				for (j = i + 1; j < r; ++j)
				{
					if (res.m[j][master] != T(0))
					{
						swap(res.m[i], res.m[j]);
						break;
					}
				}
				if (j == r)
				{
					master++;
				}
			}

			//cout << res << endl;

			if (master == c)
			{
				break;
			}

			for (int j = i + 1; j < r; ++j)
			{
				//cout << i << " " << j << endl;
				if (res.m[j][master] != T(0))
				{
					T s = res.m[j][master] * T(-1) / res.m[i][master];
					for (int k = master; k < c; ++k)
					{
						res.m[j][k] += (res.m[i][k] * s);
					}
				}
				//cout << res << endl;
			}

			for (int j = i - 1; j >= 0; --j)
			{
				if (res.m[j][master] != T(0))
				{
					T s = res.m[j][master] * T(-1) / res.m[i][master];
					for (int k = master; k < c; ++k)
					{
						res.m[j][k] += (res.m[i][k] * s);
					}
				}
				//cout << res << endl;
			}

			T head = res.m[i][master];
			for (int j = master; j < c; ++j)
			{
				res.m[i][j] /= head;
			}
			//cout << res << endl;
		}

		return res;
	}

	//��ȡͬ�͵�λ����
	Matrix<T> unit() const
	{
		assert(r == c);

		Matrix<T> res(r, r);
		for (int i = 0; i < r; ++i)
		{
			res.m[i][i] = T(1);
		}

		return res;
	}

	//��ȡͬ�������
	Matrix<T> zero() const
	{
	    assert(r == c);

        return Matrix(r, r, T(0));
	}

	//��������ͬ�����ľ���ϲ�Ϊһ�������
	static Matrix<T> merge(const Matrix<T> &A, const Matrix<T> &B)
	{
		assert(A.r == B.r);

		Matrix<T> res(A.r, A.c + B.c);
		for (int i = 0; i < A.r; ++i)
		{
			for (int j = 0; j < A.c; ++j)
			{
				res.m[i][j] = A.m[i][j];
			}
			for (int j = 0; j < B.c; ++j)
			{
				res.m[i][j + A.c] = B.m[i][j];
			}
		}

		return res;
	}

	//�жϾ����ǰ���Ƿ������λ���� (�����ж����󷽳��Ƿ��н�)
	static bool includeUnit(const Matrix<T> &A, int r)
    {
        assert(A.r >= r && A.c >= r);

        for (int i = 0; i < r; ++i)
        {
            for (int j = 0; j < r; ++j)
            {
                if (i == j)
                {
                    if (A.m[i][j] != T(1))
                    {
                        return false;
                    }
                }
                else
                {
                    if (A.m[i][j] != T(0))
                    {
                        return false;
                    }
                }
            }
        }

        return true;
    }

	//����� ��������������ʾ������
	Matrix<T> inverse()
	{
        assert(r == c);
        Matrix<T> res(r, r);
        if (solve(*this, unit(), res))
        {
            return res;
        }
        else
        {
            return zero();
        }
	}

	//����󷽳� AX=B
	static bool solve(const Matrix<T> &A, const Matrix<T> &B, Matrix<T> &X)
	{
		assert(A.r == A.c);
		assert(A.r == B.r);
		assert(A.c == X.r && X.c == B.c);

		Matrix<T> t = merge(A, B);
		//cout << t << endl;

		//cout << t << endl;

		t = t.leastLine();

		//cout << t << endl;

		//�ж�ǰ�벿���Ƿ�Ϊ��λ����
		if (includeUnit(t, A.r))
        {
            for (int i = 0; i < X.r; ++i)
            {
                for (int j = 0; j < X.c; ++j)
                {
                    X.m[i][j] = t.m[i][j + A.c];
                }
            }
            return true;
        }
        else
        {
            return false;
        }
	}
};

//�ж��Ƿ������λ����
/*template<typename T>
bool Matrix<T>::includeUnit(const Matrix<T> &A, int r)
{
    assert(A.r >= r && A.c >= r);

    for (int i = 0; i < r; ++i)
    {
        for (int j = 0; j < r; ++j)
        {
            if (i == j)
            {
                if (A.m[i][j] != T(1))
                {
                    return false;
                }
            }
            else
            {
                if (A.m[i][j] != T(0))
                {
                    return false;
                }
            }
        }
    }

    return true;
}*/

//�ж��Ƿ������λ���� ��double���⴦��
template<>
bool Matrix<double>::includeUnit(const Matrix<double> &A, int r)
{
    assert(A.r >= r && A.c >= r);

    for (int i = 0; i < r; ++i)
    {
        for (int j = 0; j < r; ++j)
        {
            if (i == j)
            {
                //if (A.m[i][j] != T(1))
                //if (fabs(A.m[i][j] - 1.0) > 1e-10)
                if (!DoubleEqual(A.m[i][j], 1.0))
                {
                    return false;
                }
            }
            else
            {
                //if (A.m[i][j] != T(0))
                //if (fabs(A.m[i][j] - 0.0) > 1e-10)
                if (!DoubleEqual(A.m[i][j], 0.0))
                {
                    return false;
                }
            }
        }
    }

    return true;
}

//�ж��Ƿ������λ���� ��float���⴦��
template<>
bool Matrix<float>::includeUnit(const Matrix<float> &A, int r)
{
    assert(A.r >= r && A.c >= r);

    //cout << A.toString() << endl;

    for (int i = 0; i < r; ++i)
    {
        for (int j = 0; j < r; ++j)
        {
            if (i == j)
            {
                //if (A.m[i][j] != T(1))
                //if (fabs(A.m[i][j] - 1.0) > 1e-10)
                if (!FloatEqual(A.m[i][j], 1.0))
                {
                    //cout << "f1" << endl;
                    return false;
                }
            }
            else
            {
                //if (A.m[i][j] != T(0))
                //if (fabs(A.m[i][j] - 0.0) > 1e-10)
                if (!FloatEqual(A.m[i][j], 0.0))
                {
                    //cout << A.m[i][j] - 0.0 << endl;
                    //cout << "f2" << endl;
                    return false;
                }
            }
        }
    }

    return true;
}

//����== ��double���⴦��
template<>
bool Matrix<double>::operator==(const Matrix<double> &mat) const
{
    if (r != mat.r || c != mat.c)
    {
        return false;
    }

    for (int i = 0; i < r; ++i)
    {
        for (int j = 0; j < c; ++j)
        {
            //if (m[i][j] != mat.m[i][j])
            //if (fabs(m[i][j] - mat.m[i][j]) > 1e-10)
            if (!DoubleEqual(m[i][j], mat.m[i][j]))
            {
                return false;
            }
        }
    }

    return true;
}

//����== ��float���⴦��
template<>
bool Matrix<float>::operator==(const Matrix<float> &mat) const
{
    if (r != mat.r || c != mat.c)
		{
			return false;
		}

		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				//if (m[i][j] != mat.m[i][j])
				//if (fabs(m[i][j] - mat.m[i][j]) > 1e-10)
				if (!FloatEqual(m[i][j], mat.m[i][j]))
				{
					return false;
				}
			}
		}

		return true;
}

//����ʽ ��double���⴦��
template<>
double Matrix<double>::det() const
{
    assert(r == c);

    double res = 1.0;

    Matrix<double> t = *this;
    for (int i = 0; i < c - 1; ++i) //����i�������Ԫ�ر�Ϊ0
    {
        //if (t.m[i][i] == T(0))
        //if (fabs(t.m[i][i] - 0.0) < 1e-10)
        if (DoubleEqual(t.m[i][i], 0.0))
        {
            int j;
            for (j = i + 1; j < r; ++j)
            {
                //if (t.m[j][i] != T(0))
                //if (fabs(t.m[j][i] - 0.0) > 1e-10)
                if (!DoubleEqual(t.m[j][i], 0.0))
                {
                    swap(t.m[i], t.m[j]);
                    res *= -1.0;
                    break;
                }
            }
            if (j == r)
            {
                return res = 0.0;
            }
        }

        for (int j = i + 1; j < r; ++j)
        {
            //if (t.m[j][i] != T(0))
            //if (fabs(t.m[j][i] - 0.0) > 1e-10)
            if (!DoubleEqual(t.m[j][i], 0.0))
            {
                double s = t.m[j][i] * (-1.0) / t.m[i][i];
                for (int k = 0; k < c; ++k)
                {
                    t.m[j][k] += (s * t.m[i][k]);
                }
            }
            //cout << t << endl;
        }
    }

    for (int i = 0; i < c; ++i)
    {
        res *= (t.m[i][i]);
    }

    //cout << t << endl;

    return res;
}

//����ʽ ��float���⴦��
template<>
float Matrix<float>::det() const
{
    assert(r == c);

    float res = 1.0;

    Matrix<float> t = *this;
    for (int i = 0; i < c - 1; ++i) //����i�������Ԫ�ر�Ϊ0
    {
        //if (t.m[i][i] == T(0))
        //if (fabs(t.m[i][i] - 0.0) < 1e-10)
        if (FloatEqual(t.m[i][i], 0.0))
        {
            int j;
            for (j = i + 1; j < r; ++j)
            {
                //if (t.m[j][i] != T(0))
                //if (fabs(t.m[j][i] - 0.0) > 1e-10)
                if (!FloatEqual(t.m[j][i], 0.0))
                {
                    swap(t.m[i], t.m[j]);
                    res *= -1.0;
                    break;
                }
            }
            if (j == r)
            {
                return res = 0.0;
            }
        }

        for (int j = i + 1; j < r; ++j)
        {
            //if (t.m[j][i] != T(0))
            //if (fabs(t.m[j][i] - 0.0) > 1e-10)
            if (!FloatEqual(t.m[j][i], 0.0))
            {
                double s = t.m[j][i] * (-1.0) / t.m[i][i];
                for (int k = 0; k < c; ++k)
                {
                    t.m[j][k] += (s * t.m[i][k]);
                }
            }
            //cout << t << endl;
        }
    }

    for (int i = 0; i < c; ++i)
    {
        res *= (t.m[i][i]);
    }

    //cout << t << endl;

    return res;
}

//���ȱ任 ������Ϊ������� ��double���⴦��
template<>
Matrix<double> Matrix<double>::leastLine()
{
    Matrix<double> res = *this;
    //cout << res.r << " " << res.c << endl;
    //cout << res.m[1][0] << endl;

    for (int i = 0; i < r; ++i)
    {
        int master = i; //��Ԫ����
        //while (master < c && res.m[i][master] == T(0))
        //while (master <c && fabs(res.m[i][master] - 0.0) < 1e-10)
        while (master < c && DoubleEqual(res.m[i][master], 0.0))
        {
            int j;
            for (j = i + 1; j < r; ++j)
            {
                //if (res.m[j][master] != T(0))
                //if (fabs(res.m[j][master] - 0.0) > 1e-10)
                if (!DoubleEqual(res.m[j][master], 0.0))
                {
                    swap(res.m[i], res.m[j]);
                    break;
                }
            }
            if (j == r)
            {
                master++;
            }
        }

        //cout << res << endl;

        if (master == c)
        {
            break;
        }

        for (int j = i + 1; j < r; ++j)
        {
            //cout << i << " " << j << endl;
            //if (res.m[j][master] != T(0))
            //if (fabs(res.m[j][master] - 0.0) > 1e-10)
            if (!DoubleEqual(res.m[j][master], 0.0))
            {
                double s = res.m[j][master] * (-1.0) / res.m[i][master];
                for (int k = master; k < c; ++k)
                {
                    res.m[j][k] += (res.m[i][k] * s);
                }
            }
            //cout << res << endl;
        }

        for (int j = i - 1; j >= 0; --j)
        {
            //if (res.m[j][master] != T(0))
            //if (fabs(res.m[j][master] - 0.0) > 1e-10)
            if (!DoubleEqual(res.m[j][master], 0.0))
            {
                double s = res.m[j][master] * (-1.0) / res.m[i][master];
                for (int k = master; k < c; ++k)
                {
                    res.m[j][k] += (res.m[i][k] * s);
                }
            }
            //cout << res << endl;
        }

        double head = res.m[i][master];
        for (int j = master; j < c; ++j)
        {
            res.m[i][j] /= head;
        }
        //cout << res << endl;
    }

    return res;
}

//���ȱ任 ������Ϊ������� ��float���⴦��
template<>
Matrix<float> Matrix<float>::leastLine()
{
    Matrix<float> res = *this;
    //cout << res.r << " " << res.c << endl;
    //cout << res.m[1][0] << endl;

    for (int i = 0; i < r; ++i)
    {
        int master = i; //��Ԫ����
        //while (master < c && res.m[i][master] == T(0))
        //while (master <c && fabs(res.m[i][master] - 0.0) < 1e-10)
        while (master < c && FloatEqual(res.m[i][master], 0.0))
        {
            int j;
            for (j = i + 1; j < r; ++j)
            {
                //if (res.m[j][master] != T(0))
                //if (fabs(res.m[j][master] - 0.0) > 1e-10)
                if (!FloatEqual(res.m[j][master], 0.0))
                {
                    swap(res.m[i], res.m[j]);
                    break;
                }
            }
            if (j == r)
            {
                master++;
            }
        }

        //cout << res << endl;

        if (master == c)
        {
            break;
        }

        for (int j = i + 1; j < r; ++j)
        {
            //cout << i << " " << j << endl;
            //if (res.m[j][master] != T(0))
            //if (fabs(res.m[j][master] - 0.0) > 1e-10)
            if (!FloatEqual(res.m[j][master], 0.0))
            {
                float s = res.m[j][master] * (-1.0) / res.m[i][master];
                for (int k = master; k < c; ++k)
                {
                    res.m[j][k] += (res.m[i][k] * s);
                }
            }
            //cout << res << endl;
        }

        for (int j = i - 1; j >= 0; --j)
        {
            //if (res.m[j][master] != T(0))
            //if (fabs(res.m[j][master] - 0.0) > 1e-10)
            if (!FloatEqual(res.m[j][master], 0.0))
            {
                float s = res.m[j][master] * (-1.0) / res.m[i][master];
                for (int k = master; k < c; ++k)
                {
                    res.m[j][k] += (res.m[i][k] * s);
                }
            }
            //cout << res << endl;
        }

        float head = res.m[i][master];
        for (int j = master; j < c; ++j)
        {
            res.m[i][j] /= head;
        }
        //cout << res << endl;
    }

    return res;
}

//�������
template<typename T>
ostream& operator<<(ostream &out, const Matrix<T> &mat)
{
	out << mat.toString();
	return out;
}

//����* (�����˾���)
template<typename T>
Matrix<T> operator*(const T &k, const Matrix<T> &mat)
{
	Matrix<T> res(mat.r, mat.c);
	for (int i = 0; i < mat.r; ++i)
	{
		for (int j = 0; j < mat.c; ++j)
		{
			res.m[i][j] = k * mat.m[i][j];
		}
	}

	return res;
}
