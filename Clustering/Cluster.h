#ifndef CLUSTER_H
#define CLUSTER_H

#include "SpaceObject.h"
#include <iostream>
#include <string>
using namespace std;

class Cluster : public SpaceObject
{
public:
	Cluster();
	Cluster(string);
	Cluster(string, double, double, double, double);
	virtual void setRA(double newRA){ra=newRA;}
	virtual void setDec(double newDec){dec=newDec;}
	virtual double getRA(){return ra;}
	virtual double getDec(){return dec;}
	friend istream& operator >>(istream& ins, Cluster& source);
	friend ostream& operator <<(ostream& outs, const Cluster& source);
	virtual double getTheta();
	virtual double getPhi();
private:
	double ra, dec;
};

#endif
