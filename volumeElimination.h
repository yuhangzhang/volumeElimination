#include <map>
#include <vector>
#include "../poly/poly.h"
#include "../QPBO-v1.31.src/QPBO.h"

using namespace std;

#ifndef _VOLUMEELIMINATION_H_
#define _VOLUMEELIMINATION_H_


class volumeElimination
{
public:
	typedef vector<vector<vector<double>>> vector3d;
	typedef vector<vector<vector<bool>>> vector3b;

	volumeElimination(int row, int col, int height);
	void addDataterm(int i, int j, int k, double v);
	double getDataterm(int i, int j, int k);
	void addEdgeterm(int i, int j, int k, int i2, int j2, int k2, double v);
	double getEdgeterm(int i, int j, int k, int i2, int j2, int k2);
	void minimize();
	void setLabel(int i, int j, int k, bool v);
	bool getLabel(int i, int j, int k);
	vector3b getLabel();
	double energy();

private:
	vector3d _dataterm;
	poly<double> _edgeterm;

	vector3b _label;


};


#endif