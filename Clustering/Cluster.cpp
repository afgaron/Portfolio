#include "SpaceObject.h"
#include "Cluster.h"
#include <iostream>
#include <string>
using namespace std;

const double PI=3.1415926535;

Cluster::Cluster()
{
	name="no_name";
	catalog="no_catalog";
	ra=0;
	dec=0;
	z=0;
	zErr=0;
	calcDistance();
}

Cluster::Cluster(string newName)
{
	name=newName;
	setCatalog();
	ra=0;
	dec=0;
	z=0;
	zErr=0;
	calcDistance();
}

Cluster::Cluster(string newName, double newRA, double newDec, double newZ, double newZErr)
{
	name=newName;
	setCatalog();
	ra=newRA;
	dec=newDec;
	z=newZ;
	zErr=newZErr;
	calcDistance();
}

istream& operator >>(istream& ins, Cluster& group)
{
	ins >> group.name >> group.ra >> group.dec >> group.z >> group.zErr;
	group.setCatalog();
	group.calcDistance();
	return ins;
}

ostream& operator <<(ostream& outs, const Cluster& group)
{
	outs << group.name << '\t' << group.ra << '\t' << group.dec << '\t' << group.z << '\t' << group.zErr << '\t' << group.ADD << '\t' << group.CRD;
	return outs;
}

double Cluster::getTheta()
{
	return (90.-dec) * PI / 180;
}

double Cluster::getPhi()
{
	return ra * PI / 180;
}
