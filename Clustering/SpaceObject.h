#ifndef SPACEOBJECT_H
#define SPACEOBJECT_H

#include <string>
using namespace std;

class SpaceObject
{
public:
	SpaceObject();
	SpaceObject(string);
	SpaceObject(string, double, double);
	void setName(string newName){name=newName;}
	void setCatalog();
	void setZ(double newZ){z=newZ;}
	void setZErr(double newZErr){zErr=newZErr;}
	void setADD(double newADD){ADD=newADD;}
	void setCRD(double newCRD){CRD=newCRD;}
	virtual void setRA(double newRA);
	virtual void setDec(double newDec);
	string getName(){return name;}
	string getCatalog(){return catalog;}
	double getZ(){return z;}
	double getZErr(){return zErr;}
	double getADD(){return ADD;}
	double getCRD(){return CRD;}
	virtual double getTheta();
	virtual double getPhi();
	virtual double getRA();
	virtual double getDec();
	void calcDistance();
protected:
	string name, catalog;
	double z, zErr, ADD, CRD;
};

#endif
