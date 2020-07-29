#ifndef RADIO_H
#define RADIO_H

#include "SpaceObject.h"
#include <iostream>
#include <string>
using namespace std;

class Radio : public SpaceObject
{
public:
	Radio();
	Radio(string);
	Radio(string, double, double, double, double, double, double, double, double, double, double, double);
	void setBent(double newBent){bent=newBent;}
	virtual void setRA(double newRA);
	virtual void setDec(double newDec);
	void setRA1(double newRA1){ra1=newRA1;}
	void setDec1(double newDec1){dec1=newDec1;}
	void setRA2(double newRA2){ra2=newRA2;}
	void setDec2(double newDec2){dec2=newDec2;}
	void setRA3(double newRA3){ra3=newRA3;}
	void setDec3(double newDec3){dec3=newDec3;}
	void setPhotoZ(double newPhotoZ){photo_z=newPhotoZ;}
	void setSpecZ(double newSpecZ){spec_z=newSpecZ;}
	double getBent(){return bent;}
	virtual double getRA();
	virtual double getDec();
	double getRA1(){return ra1;}
	double getDec1(){return dec1;}
	double getRA2(){return ra2;}
	double getDec2(){return dec2;}
	double getRA3(){return ra3;}
	double getDec3(){return dec3;}
	double getPhotoZ(){return photo_z;}
	double getSpecZ(){return spec_z;}
	friend istream& operator >>(istream& ins, Radio& source);
	friend ostream& operator <<(ostream& outs, const Radio& source);
	virtual double getTheta();
	virtual double getPhi();
private:
	double bent, ra1, dec1, ra2, dec2, ra3, dec3, photo_z, spec_z;
};

#endif
