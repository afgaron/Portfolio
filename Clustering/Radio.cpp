#include "SpaceObject.h"
#include "Radio.h"
#include <iostream>
#include <string>
#include <cmath>
using namespace std;

const double PI=3.1415926535;

int center(double ra1, double dec1, double ra2, double dec2, double ra3, double dec3);
double distance(double ra1, double dec1, double ra2, double dec2);

Radio::Radio()
{
	name="no_name";
	catalog="no_catalog";
	bent=0;
	ra1=0;
	dec1=0;
	ra2=0;
	dec2=0;
	ra3=0;
	dec3=0;
	photo_z=0;
	spec_z=0;
	z=0;
	zErr=0;
	calcDistance();
}

Radio::Radio(string newName)
{
	name=newName;
	setCatalog();
	bent=0;
	ra1=0;
	dec1=0;
	ra2=0;
	dec2=0;
	ra3=0;
	dec3=0;
	photo_z=0;
	spec_z=0;
	z=0;
	zErr=0;
	calcDistance();
}

Radio::Radio(string newName, double newBent, double newRA1, double newDec1, double newRA2, double newDec2, double newRA3, double newDec3, double newPhotoZ, double newSpecZ, double newBestZ, double newZErr)
{
	name=newName;
	setCatalog();
	bent=newBent;
	ra1=newRA1;
	dec1=newDec1;
	ra2=newRA2;
	dec2=newDec2;
	ra3=newRA3;
	dec3=newDec3;
	photo_z=newPhotoZ;
	spec_z=newSpecZ;
	z=newBestZ;
	zErr=newZErr;
	calcDistance();
}

void Radio::setRA(double newRA)
{
	ra1=newRA;
	ra2=newRA;
	ra3=newRA;
}

void Radio::setDec(double newDec)
{
	dec1=newDec;
	dec2=newDec;
	dec3=newDec;
}

istream& operator >>(istream& ins, Radio& source)
{
	ins >> source.name >> source.bent >> source.ra1 >> source.dec1 >> source.ra2 >> source.dec2 >> source.ra3 >> source.dec3 >> source.photo_z >> source.spec_z >> source.z >> source.zErr;
	source.setCatalog();
	source.calcDistance();
	return ins;
}

ostream& operator <<(ostream& outs, const Radio& source)
{
	outs << source.name << '\t' << source.bent << '\t';
	int centerPt = center(source.ra1, source.dec1, source.ra2, source.dec2, source.ra3, source.dec3);
	switch(centerPt)
	{
		case 1:
			outs << source.ra1 << '\t' << source.dec1 << '\t' << source.ra2 << '\t' << source.dec2 << '\t' << source.ra3 << '\t' << source.dec3;
			break;
		case 2:
			outs << source.ra2 << '\t' << source.dec2 << '\t' << source.ra1 << '\t' << source.dec1 << '\t' << source.ra3 << '\t' << source.dec3;
			break;
		case 3:
			outs << source.ra3 << '\t' << source.dec3 << '\t' << source.ra2 << '\t' << source.dec2 << '\t' << source.ra1 << '\t' << source.dec1;
	}
	outs << '\t' << source.photo_z << '\t' << source.spec_z << '\t' << source.z << '\t' << source.zErr << '\t' << source.ADD << '\t' << source.CRD;
	return outs;
}

double Radio::getRA()
{
	int centerPt = center(ra1, dec1, ra2, dec2, ra3, dec3);
	switch(centerPt)
	{
		case 1: return ra1;
		case 2: return ra2;
		case 3: return ra3;
	}
}

double Radio::getDec()
{
	int centerPt = center(ra1, dec1, ra2, dec2, ra3, dec3);
	switch(centerPt)
	{
		case 1: return dec1;
		case 2: return dec2;
		case 3: return dec3;
	}
}

double Radio::getTheta()
{
	int centerPt = center(ra1, dec1, ra2, dec2, ra3, dec3);
	switch(centerPt)
	{
		case 1: return (90.-dec1) * PI / 180;
		case 2: return (90.-dec2) * PI / 180;
		case 3: return (90.-dec3) * PI / 180;
	}
}

double Radio::getPhi()
{
	int centerPt = center(ra1, dec1, ra2, dec2, ra3, dec3);
	switch(centerPt)
	{
		case 1: return ra1 * PI / 180;
		case 2: return ra2 * PI / 180;
		case 3: return ra3 * PI / 180;
	}
}

int center(double ra1, double dec1, double ra2, double dec2, double ra3, double dec3)
{
	double d1, d2, d3;
	d1 = distance(ra2, dec2, ra3, dec3);
	d2 = distance(ra1, dec1, ra3, dec3);
	d3 = distance(ra1, dec1, ra2, dec2);

	if(d1 > d2 && d1 > d3)
	{
		return 1;
	}
	else if(d2 > d1 && d2 > d3)
	{
		return 2;
	}
	else
	{
		return 3;
	}
}

double distance(double ra1, double dec1, double ra2, double dec2)
{
	return sqrt((ra1-ra2) * (ra1-ra2) + (dec1-dec2) * (dec1-dec2));
}
