#include <iostream>
#include "Matrix.h"
#include "UnitTest.h"

using namespace std;

int main()
{
    /*单元测试*/

    //去掉以下注释可运行单元测试
    //单元测试需要较长时间，请耐心等待
    /*TestAdd();
    TestSub();
    TestMul();
    TestSca();
    TestDet();
    TestRow();
    TestInv();
    TestPow();*/

    /*用法示例*/

    double a[3][3] =
    {
        1, 2, 3,
        9, 8, 7,
        3, 2, 5
    };

    Matrix<double> mat((double*)a, 3, 3); //利用二维数组构造矩阵
    mat.setElem(1, 1, 0); //将矩阵第1行第1列设置为0 (矩阵行数和列数从0开始)
    cout << "mat = " << endl << mat << endl; //输出矩阵
    cout << "mat + mat = " << endl << mat + mat << endl; //加法
    cout << "mat - mat = " << endl << mat - mat << endl; //减法
    cout << "mat * mat = " << endl << mat * mat << endl; //乘法
    cout << "-6 * mat = " << endl << (-6.0) * mat << endl; //数乘
    cout << "mat ^ 3 = " << endl << (mat ^ 3) << endl; //求幂
    cout << "det = " << mat.det() << endl << endl; //行列式
    cout << "inv = " << endl << mat.inverse() << endl; //逆矩阵
    cout << "row elem = " << endl << mat.leastLine() << endl; //行最简型
    cout << "transpose = " << endl << mat.transpose() << endl; //转置
    cout << "complement = " << endl << mat.complement(0, 1) << endl; //余子式

    return 0;
}
