#include <iostream>
#include <fstream>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <list>
#include <vector>
#include <algorithm>
#include <ctime>
#include "SpaceObject.h"
#include "Cluster.h"
#include "Radio.h"
#include "MatchedPair.h"
using namespace std;

const double PI=3.1415926535;

double separation(SpaceObject&, SpaceObject&); //returns physical distance
double sphereSep(SpaceObject&, SpaceObject&); //returns apparent distance in plane of sky
double losDistance(SpaceObject&, SpaceObject&); //returns distance along line-of-sight
double zDiff(SpaceObject&, SpaceObject&); //returns difference in redshift
bool zDiffMax(SpaceObject&, SpaceObject&); //returns true if one redshift is within the other's error
double degrees(SpaceObject&, SpaceObject&); //returns arc separation
void makePair(list<MatchedPair>&, list<Cluster>&, list<Radio>&, int crit, int max); //pairs clusters with radio sources
void unDuplicate(list<MatchedPair>&); //removes duplicate radio matches; only closest pair remains
void printList(list<MatchedPair>&, ofstream&); //outputs pairs to file

int main()
{
	time_t begin=time(0);

	ifstream cluster_ins, radio_ins;
	ofstream mp_outs, control_outs;
	cluster_ins.open("CLUSTER_INPUT.txt");
	if(cluster_ins.fail())
	{
		cout << "Cluster input file opening failed.\n";
		exit(1);
	}
	radio_ins.open("RADIO_INPUT.txt");
	if(radio_ins.fail())
	{
		cout << "Radio input file opening failed.\n";
		exit(1);
	}
	mp_outs.open("MATCHED_PAIRS.xls");
	if(mp_outs.fail())
	{
		cout << "Output file opening failed.\n";
		exit(1);
	}

	double maxSep;
	while(true)
	{
		cout << "Max separation (in Mpc): ";
		cin >> maxSep;
		if(cin.good())
			break;
		cout << "Invalid choice.\n";
		cin.clear();
		while(cin.get() != '\n'){}
	}

	cout << "Criteria for matching:\n    (1) Max physical separation of " << maxSep << " Mpc\n    (2) Max separation of " << maxSep << " Mpc in plane of the sky and z within 0.01\n    (3) Max of " << maxSep << " Mpc in plane and z within given error\n    (4) Max of " << maxSep << " Mpc in plane and no z condition\n";
	int criterion;
	while(true)
	{
		cout << "Choose criterion 1-4: ";
		cin >> criterion;
		if(criterion==1 || criterion==2 || criterion==3 || criterion==4)
			break;
		cout << "Invalid choice.\n";
		cin.clear();
		while(cin.get() != '\n'){}
	}

	char ctrl;
	while(true)
	{
		cout << "Control group (y/n): ";
		cin >> ctrl;
		ctrl = tolower(ctrl);
		if(ctrl=='y' || ctrl=='n')
			break;
		cout << "Invalid choice.\n";
		cin.clear();
		while(cin.get() != '\n'){}
	}

	char unDup;
	while(true)
	{
		cout << "Remove duplicates (y/n): ";
		cin >> unDup;
		unDup = tolower(unDup);
		if(unDup=='y' || unDup=='n')
			break;
		cout << "Invalid choice.\n";
		cin.clear();
		while(cin.get() != '\n'){}
	}

	string header;
	getline(cluster_ins, header);
	getline(radio_ins, header);

	list<Cluster> cluster_list;
	Cluster group;
	cluster_ins >> group;
	while(!cluster_ins.eof())
	{
		if(group.getZ() <= 0.6)
			cluster_list.push_back(group);
		cluster_ins >> group;
	}
	cluster_ins.close();

	list<Radio> radio_list;
	Radio source;
	radio_ins >> source;
	while(!radio_ins.eof())
	{
		if(source.getZ() <= 0.6)
			radio_list.push_back(source);
		radio_ins >> source;
	}
	radio_ins.close();

	header = "Cluster name\tra\tdec\tRedshift\tzErr\tADD (Mpc)\tCRD (Mpc)\tRadio name\tBent score\tra 1\t'dec 1\tra 2\t'dec 2\tra 3\t'dec 3\tPhoto z\tSpec z\tBest z\tzErr\tADD (Mpc)\tCRD (Mpc)\tDistance (Mpc)\tRadial (Mpc)\tTangential (Mpc)\tZ diff\tArc sep (°)\n";

	mp_outs << header;
	list<MatchedPair> pairs;
	makePair(pairs, cluster_list, radio_list, criterion, maxSep);
	if(unDup=='y')
		unDuplicate(pairs);
	printList(pairs, mp_outs);
	cout << pairs.size() << " matches under criterion " << criterion << endl;
	mp_outs.close();

	if(ctrl=='y')
	{
		control_outs.open("CONTROL.xls");
		if(control_outs.fail())
		{
			cout << "Control file opening failed.\n";
			exit(1);
		}

		control_outs << header;

		srand(time(0));
		vector<Radio> randomVec(radio_list.begin(), radio_list.end());
		random_shuffle(randomVec.begin(), randomVec.end());
		list<Radio> random1(randomVec.begin(), randomVec.end());
		random_shuffle(randomVec.begin(), randomVec.end());
		list<Radio> random2(randomVec.begin(), randomVec.end());

		list<Radio>::iterator rl = radio_list.begin();
		for(list<Radio>::iterator ran = random1.begin(); ran!=random1.end(); ran++)
		{
			double newRA=ran->getRA();
			rl->setRA(newRA);
			rl++;
		}
		rl = radio_list.begin();
		for(list<Radio>::iterator ran = random2.begin(); ran!=random2.end(); ran++)
		{
			double newZ=ran->getZ();
			rl->setZ(newZ);
			rl++;
		}

		list<MatchedPair> control;
		makePair(control, cluster_list, radio_list, criterion, maxSep);
		if(unDup=='y')
			unDuplicate(control);
		printList(control, control_outs);
		cout << control.size() << " matches under criterion " << criterion << endl;
		control_outs.close();
	}

	time_t end=time(0);
	cout << "Time elapsed: " << end-begin << " sec" << endl;

	return 0;
}

double separation(SpaceObject& arg1, SpaceObject& arg2)
{
	double arcSep, lawCos;
	arcSep = cos(arg1.getTheta())*cos(arg2.getTheta()) + sin(arg1.getTheta())*sin(arg2.getTheta())*cos( arg1.getPhi()-arg2.getPhi() );
	lawCos = sqrt( arg1.getCRD()*arg1.getCRD() + arg2.getCRD()*arg2.getCRD() - 2*arg1.getCRD()*arg2.getCRD()*arcSep );
	return lawCos;
}

double sphereSep(SpaceObject& arg1, SpaceObject& arg2)
{
	double arcSep, cosApprox;
	arcSep = cos(arg1.getTheta())*cos(arg2.getTheta()) + sin(arg1.getTheta())*sin(arg2.getTheta())*cos( arg1.getPhi()-arg2.getPhi() );
	cosApprox = (arg1.getCRD() + arg2.getCRD()) * sqrt((1. - arcSep)/4.);
	return cosApprox;
}

double losDistance(SpaceObject& arg1, SpaceObject& arg2)
{
	return abs(arg1.getCRD()-arg2.getCRD());
}

double zDiff(SpaceObject& arg1, SpaceObject& arg2)
{
	return abs(arg1.getZ()-arg2.getZ());
}

bool zDiffMax(SpaceObject& arg1, SpaceObject& arg2)
{
	bool z1 = arg1.getZErr()>=zDiff(arg1, arg2);
	bool z2 = arg2.getZErr()>=zDiff(arg1, arg2);
	if(arg1.getZErr()==0 && arg2.getZErr()==0)
		return true;
	else
		return (z1 || z2);
}

double degrees(SpaceObject& arg1, SpaceObject& arg2)
{
	double arcSep;
	arcSep = cos(arg1.getTheta())*cos(arg2.getTheta()) + sin(arg1.getTheta())*sin(arg2.getTheta())*cos(arg1.getPhi()-arg2.getPhi());
	return acos(arcSep) * 180/PI;
}

void makePair(list<MatchedPair>& pairList, list<Cluster>& cList, list<Radio>& rList, int crit, int max)
{
	for(list<Cluster>::iterator i = cList.begin(); i!=cList.end(); i++)
	{
		for(list<Radio>::iterator j = rList.begin(); j!=rList.end(); j++)
		{
			switch(crit)
			{
				case 1:
					if(separation(*i, *j)<=max)
					{
						MatchedPair pair(*i, *j);
						pairList.push_back(pair);
					}
					break;
				case 2:
					if(sphereSep(*i, *j)<=max && zDiff(*i, *j)<=0.01)
					{
						MatchedPair pair(*i, *j);
						pairList.push_back(pair);
					}
					break;
				case 3:
					if(sphereSep(*i, *j)<=max && zDiffMax(*i, *j))
					{
						MatchedPair pair(*i, *j);
						pairList.push_back(pair);
					}
					break;
				case 4:
					if(sphereSep(*i, *j)<=max)
					{
						MatchedPair pair(*i, *j);
						pairList.push_back(pair);
					}
			}
		}
	}
}

void unDuplicate(list<MatchedPair>& pList)
{
	int i_dup=0;
	for(list<MatchedPair>::iterator i = pList.begin(); i!=pList.end();)
	{
		for(list<MatchedPair>::iterator j = pList.begin(); j!=pList.end();)
		{
			if(i!=j && i->getRadio().getName() == j->getRadio().getName())
			{
				Cluster c_i=i->getCluster(), c_j=j->getCluster();
				Radio r_i=i->getRadio(), r_j=j->getRadio();
				if(separation(c_i, r_i) < separation(c_j, r_j))
					pList.erase(j++);
				else
				{
					pList.erase(i++);
					i_dup++;
				}
			}
			else
				j++;
		}
		if(i_dup==0)
			i++;
		else
			i_dup=0;
	}
}

void printList(list<MatchedPair>& pList, ofstream& outs)
{
	for(list<MatchedPair>::iterator p = pList.begin(); p!=pList.end(); p++)
	{
		Cluster cTemp=p->getCluster();
		Radio rTemp=p->getRadio();
		outs << *p << '\t' << separation(cTemp, rTemp) << '\t' << losDistance(cTemp, rTemp) << '\t' << sphereSep(cTemp, rTemp) << '\t' << zDiff(cTemp, rTemp) << '\t' << degrees(cTemp, rTemp) << endl;
	}
}
