#include <iostream>
#include "Matrix.h"
#include "UnitTest.h"

using namespace std;

int main()
{
    /*��Ԫ����*/

    //ȥ������ע�Ϳ����е�Ԫ����
    //��Ԫ������Ҫ�ϳ�ʱ�䣬�����ĵȴ�
    /*TestAdd();
    TestSub();
    TestMul();
    TestSca();
    TestDet();
    TestRow();
    TestInv();
    TestPow();*/

    /*�÷�ʾ��*/

    double a[3][3] =
    {
        1, 2, 3,
        9, 8, 7,
        3, 2, 5
    };

    Matrix<double> mat((double*)a, 3, 3); //���ö�ά���鹹�����
    mat.setElem(1, 1, 0); //�������1�е�1������Ϊ0 (����������������0��ʼ)
    cout << "mat = " << endl << mat << endl; //�������
    cout << "mat + mat = " << endl << mat + mat << endl; //�ӷ�
    cout << "mat - mat = " << endl << mat - mat << endl; //����
    cout << "mat * mat = " << endl << mat * mat << endl; //�˷�
    cout << "-6 * mat = " << endl << (-6.0) * mat << endl; //����
    cout << "mat ^ 3 = " << endl << (mat ^ 3) << endl; //����
    cout << "det = " << mat.det() << endl << endl; //����ʽ
    cout << "inv = " << endl << mat.inverse() << endl; //�����
    cout << "row elem = " << endl << mat.leastLine() << endl; //�������
    cout << "transpose = " << endl << mat.transpose() << endl; //ת��
    cout << "complement = " << endl << mat.complement(0, 1) << endl; //����ʽ

    return 0;
}
