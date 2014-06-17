#include "volumeElimination.h"


volumeElimination::volumeElimination(int row, int col, int height)
{
	_dataterm.resize(row);
	_label.resize(row);

	for(int i=0;i<row;i++)
	{
		_dataterm[i].resize(col);
		_label[i].resize(col);

		for(int j=0;j<col;j++)
		{
			_dataterm[i][j].resize(height,0);
			_label[i][j].resize(height,0);
		}
	}

	return;
}

void volumeElimination::addDataterm(int i, int j, int k, double v)
{
	_dataterm[i][j][k] = v;

	return;
}

double volumeElimination::getDataterm(int i, int j, int k)
{
	return _dataterm[i][j][k];
}

void volumeElimination::addEdgeterm(int i, int j, int k, int i2, int j2, int k2, double v)
{
	vector<int> index;

	index.push_back(i);
	index.push_back(j);
	index.push_back(k);
	index.push_back(i2);
	index.push_back(j2);
	index.push_back(k2);

	_edgeterm.addTerm(index,v);

	return;
}

double volumeElimination::getEdgeterm(int i, int j, int k, int i2, int j2, int k2)
{
	vector<int> index;

	index.push_back(i);
	index.push_back(j);
	index.push_back(k);
	index.push_back(i2);
	index.push_back(j2);
	index.push_back(k2);

	return _edgeterm.getTerm(index);
}

void volumeElimination::minimize()
{
	int numnode = _dataterm.size()*_dataterm[0].size()*_dataterm[0][0].size();
	QPBO<double>* solver = new QPBO<double>(numnode,numnode*3);

	solver->AddNode(numnode);

	for(int i=0;i<_dataterm.size();i++)
	{
		for(int j=0;j<_dataterm[i].size();j++)
		{
			for(int k=0;k<_dataterm[i][j].size();k++)
			{
				solver->AddUnaryTerm(i*_dataterm[i].size()+j*_dataterm[i][j].size()+k,0,_dataterm[i][j][k]);
			}
		}
	}

	for(poly<double>::TERMS::iterator it = _edgeterm.firstTerm();it!=_edgeterm.lastTerm();it++)
	{
		vector<int> index = it->first;

		solver->AddPairwiseTerm(
			index[0]*_dataterm[0].size()+index[1]*_dataterm[0][0].size()+index[2],
			index[3]*_dataterm[0].size()+index[4]*_dataterm[0][0].size()+index[5],
			0,0,0,it->second
			);
		
	}

	solver->Solve();
	solver->ComputeWeakPersistencies();

	for(int i=0;i<_label.size();i++)
	{
		for(int j=0;j<_label[i].size();j++)
		{
			for(int k=0;k<_label[i][j].size();k++)
			{
				_label[i][j][k] = solver->GetLabel(i*_label[i].size()+j*_label[i][j].size()+k);
				
			}
		}
	}

	solver->Reset();

	return;
}

void volumeElimination::setLabel(int i, int j, int k, bool v)
{
	 _label[i][j][k] = v;
	 return;
}

bool volumeElimination::getLabel(int i, int j, int k)
{
	return _label[i][j][k];
}

volumeElimination::vector3b volumeElimination::getLabel()
{
	return _label;
}

double volumeElimination::energy()
{
	double total = 0;

	for(int i=0;i<_dataterm.size();i++)
	{
		for(int j=0;j<_dataterm[i].size();j++)
		{
			for(int k=0;k<_dataterm[i][j].size();k++)
			{
				total+= _dataterm[i][j][k]*double(_label[i][j][k]);
			}
		}
	}

	for(poly<double>::TERMS::iterator it=_edgeterm.firstTerm();it!=_edgeterm.lastTerm();it++)
	{
		vector<int> index = it->first;
		total+= it->second*_label[index[0]][index[1]][index[2]]*_label[index[3]][index[4]][index[5]];
	}

	return total;
}