#ifndef MATCHEDPAIR_H
#define MATCHEDPAIR_H

#include "SpaceObject.h"
#include "Cluster.h"
#include "Radio.h"
#include <iostream>
using namespace std;

class MatchedPair
{
public:
	MatchedPair(){}
	MatchedPair(const Cluster&, const Radio&);
	MatchedPair(const Radio&, const Cluster&);
	void setCluster(const Cluster&);
	void setRadio(const Radio&);
	Cluster getCluster(){return cMatch;}
	Radio getRadio(){return rMatch;}
	friend istream& operator >>(istream&, MatchedPair&);
	friend ostream& operator <<(ostream&, const MatchedPair&);
private:
	Cluster cMatch;
	Radio rMatch;
};

#endif
