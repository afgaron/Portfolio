#include "SpaceObject.h"
#include <iostream>
#include <string>
#include <cmath>
using namespace std;

const double H0=71.;
const double h=H0/100;
const double omegaM=0.270;
const double omegaV=1.-omegaM-0.4165/(H0*H0);
const double omegaR=4.165e-5/(h*h);
const double omegaK=1.-omegaM-omegaV-omegaR;
const double c=299792.458;
const double DH=c/H0;

double calcADD(double z); //angular diameter distance
double calcTCD(double z); //transverse comoving distance
double calcCRD(double z); //comoving radial distance

SpaceObject::SpaceObject()
{
	name="no_name";
	catalog="no_catalog";
	z=0;
	zErr=0;
	calcDistance();
}

SpaceObject::SpaceObject(string newName)
{
	name=newName;
	setCatalog();
	z=0;
	zErr=0;
	calcDistance();
}

SpaceObject::SpaceObject(string newName, double newZ, double newZErr)
{
	name=newName;
	setCatalog();
	z=newZ;
	zErr=newZErr;
	calcDistance();
}

void SpaceObject::setCatalog()
{
	int pos = name.find_first_of("1234567890");
	if(name[pos-1]=='J')
		catalog = name.substr(0, pos-1);
	else
		catalog = name.substr(0, pos);
}

void SpaceObject::setRA(double newRA)
{
	cout << "RA value for " << name << " is undefined.\n";
}

void SpaceObject::setDec(double newDec)
{
	cout << "Dec value for " << name << " is undefined.\n";
}

double SpaceObject::getTheta()
{
	cout << "Theta value for " << name << " is undefined. Using 0 instead.\n";
	return 0;
}

double SpaceObject::getPhi()
{
	cout << "Phi value for " << name << " is undefined. Using 0 instead.\n";
	return 0;
}

double SpaceObject::getRA()
{
	cout << "RA value for " << name << " is undefined. Using 0 instead.\n";
	return 0;
}

double SpaceObject::getDec()
{
	cout << "Dec value for " << name << " is undefined. Using 0 instead.\n";
	return 0;
}

void SpaceObject::calcDistance()
{
	ADD = calcADD(z);
	CRD = calcCRD(z);
}

double calcADD(double z)
{
	double az, DM, DA;
	az = 1./(1.+z);
	DM = calcTCD(z);
	DA = az*DM;
	return DA;
}

double calcTCD(double z)
{
	double DCMR=calcCRD(z), ratio, x;
	x = sqrt(abs(omegaK)) * DCMR / DH;
		if(x > 0.1)
		{
			if(omegaK > 0)
				ratio = 0.5 * (exp(x) - exp(-x)) / x;
			else
				ratio = sin(x) / x;
		}
		else
		{
			double y = x*x;
			if(omegaK < 0)
				y = -y;
			ratio = 1. + y/6. + y*y/120.;
		}
	double DM = ratio*DCMR;
	return DM;
}

double calcCRD(double z)
{
	double az, DCMR=0.;
	az = 1./(1.+z);
	for(int i=0; i<1000; i++)
	{
		double a = az + (1-az) * (i+0.5) / 1000.;
		double adot = sqrt(omegaK + (omegaM/a) + (omegaR/(a*a)) + (omegaV*a*a));
		DCMR += 1./(a*adot);
	}
	DCMR = DH * (1. - az) * DCMR / 1000.;
	return DCMR;
}
