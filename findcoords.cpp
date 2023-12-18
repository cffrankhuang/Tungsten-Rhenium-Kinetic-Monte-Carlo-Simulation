#include<iostream>
using namespace std;



int main()
{
	int ltcp, x, y, z;
	int cor_x, cor_y, cor_z;

	ltcp=89594;

	x=(int)(ltcp/64)/64;
	y=(int)(ltcp/64)%64;
	z=(int)(ltcp%64);

	cor_x=-15787-x;
	cor_y=-23753-y;
	cor_z=-23942-z;

	cout << x << " " << y << " " << z << endl;
	cout << cor_x << " " << cor_y << " " << cor_z << endl;

	ltcp=41;

	x=(int)(ltcp/64)/64;
	y=(int)(ltcp/64)%64;
	z=(int)(ltcp%64);


	cout << x << " " << y << " " << z << endl;
	cout << -1088-x << " " << 2112-y << " " << -1815-z << endl;


	return 0;
}