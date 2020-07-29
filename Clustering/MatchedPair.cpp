#include "SpaceObject.h"
#include "Cluster.h"
#include "Radio.h"
#include "MatchedPair.h"
#include <iostream>
using namespace std;

MatchedPair::MatchedPair(const Cluster& cNew, const Radio& rNew)
{
	cMatch = cNew;
	rMatch = rNew;
}

MatchedPair::MatchedPair(const Radio& rNew, const Cluster& cNew)
{
	cMatch = cNew;
	rMatch = rNew;
}

void MatchedPair::setCluster(const Cluster& cNew)
{
	cMatch = cNew;
}

void MatchedPair::setRadio(const Radio& rNew)
{
	rMatch = rNew;
}

istream& operator >>(istream& ins, MatchedPair& pair)
{
	ins >> pair.cMatch >> pair.rMatch;
	return ins;
}

ostream& operator <<(ostream& outs, const MatchedPair& pair)
{
	outs << pair.cMatch << '\t' << pair.rMatch;
	return outs;
}
