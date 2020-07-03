/*
 * mapperStruct.cpp
 *
 *  Created on: May 9, 2020
 *      Author: pluto
 */

#include "mapperStruct.h"

void generate_seeds(const char *read, uint32_t tau, sFMindex nfmindx, vector<ms_seed> & vseed)
{
	uint32_t readlen = strlen(read);
	uint32_t cnt = tau + 1;
	struct ms_seed seedtmp;
	uint32_t start = 0;
	while(cnt)
	{
		uint32_t lentmp = (double)readlen / cnt + 0.5;
		seedtmp.seed_segment = new char[lentmp+1]();
		strncpy(seedtmp.seed_segment,read+start,lentmp);
		seedtmp.start_pos = start;
		seedtmp.end_pos = start + lentmp - 1;
		seedtmp.sa_ary = calc_SArangeSeq(nfmindx,seedtmp.seed_segment);
		vseed.push_back(seedtmp);
		start += lentmp;
		readlen -= lentmp;
		--cnt;
	}

}

void free_seed(vector<ms_seed> & vseed)
{
	vector<ms_seed>::iterator ite = vseed.begin();
	for(; ite != vseed.end(); ++ite)
	{
		delete [] ite->seed_segment;
		free(ite->sa_ary);
	}
}

void generate_PHNArray(struct PH_Node **p2_node, uint32_t tau)
{
	uint32_t level = ceil(log(tau+1) / log(2));
	uint32_t nodenum = pow(2,level+1);
	struct PH_Node * PHNArray = new PH_Node[nodenum];
	for(uint32_t i = 0; i < nodenum; ++i)
	{
		PHNArray[i].start_seed_id = -1;
		PHNArray[i].end_seed_id = -1;
		PHNArray[i].tau = -1;
	}
	PHNArray[1].start_seed_id = 0;		//index start by 1
	PHNArray[1].end_seed_id = tau;
	PHNArray[1].tau = tau;
	uint32_t cnt = 0;
	for(uint32_t i = 2; i < nodenum; ++i)
	{
		if(PHNArray[i/2].tau)
		{
			if(i % 2 == 0)
			{
				PHNArray[i].tau = ceil((PHNArray[i/2].tau + 1) / 2.0) - 1;
				PHNArray[i].start_seed_id = PHNArray[i/2].start_seed_id;
				PHNArray[i].end_seed_id = PHNArray[i].start_seed_id + PHNArray[i].tau;
			}
			else
			{
				PHNArray[i].tau = PHNArray[i/2].tau - PHNArray[i-1].tau - 1;
				PHNArray[i].start_seed_id = PHNArray[i-1].end_seed_id + 1;
				PHNArray[i].end_seed_id = PHNArray[i].start_seed_id + PHNArray[i].tau;
			}
			if(PHNArray[i].tau == 0)
			{
				++cnt;
				if(cnt == tau + 1)
				{
					break;
				}
			}
		}

	}
	*p2_node =  PHNArray;
	int p = 1;
	for(uint32_t i = 1; i < nodenum; ++i)
	{
		if(i == pow(2,p))
		{
			cout << endl;
			p++;
		}
		cout <<(*p2_node)[i].tau << " " << (*p2_node)[i].start_seed_id << "-" << (*p2_node)[i].end_seed_id << "\t";
	}

}

uint32_t PHNode_index(struct PH_Node *p_nodeA, struct PH_Node NodeT)
{
	uint32_t tau = p_nodeA[1].end_seed_id;
	uint32_t level = ceil(log(tau+1) / log(2));
	uint32_t nodenum = pow(2,level+1);
	for(uint32_t i = 1; i < nodenum; ++i)
	{
		if(NodeT.start_seed_id == p_nodeA[i].start_seed_id && NodeT.end_seed_id == p_nodeA[i].end_seed_id && NodeT.tau == p_nodeA[i].tau)
		{
			return i;
		}
	}
	return 0;  //not found
}

void generate_candidatePath(const vector<ms_seed> & vseed, sFMindex nfmindx, vector<candiPath> & vcadipth, uint32_t freqthld)
{
	uint32_t cnt = 0;
	struct candiPath candPathtmp;
	vector<ms_seed>::const_iterator ite = vseed.begin();
	for(; ite != vseed.end(); ++ite)
	{
		if(ite->sa_ary[1] - ite->sa_ary[0] <= freqthld)
		{
			candPathtmp.start_seed_id = cnt;
			candPathtmp.end_seed_id = cnt;
			candPathtmp.p_candidatepath = new char[strlen(ite->seed_segment) + 1]();
			strcpy(candPathtmp.p_candidatepath,ite->seed_segment);
			candPathtmp.sa_ary = new uint32_t[2];
			candPathtmp.sa_ary[0] = ite->sa_ary[0];
			candPathtmp.sa_ary[1] = ite->sa_ary[1];
			vcadipth.push_back(candPathtmp);
		}
		else
		{

		}
		++cnt;
	}
}

void generate_candidate(const vector<ms_seed> & vseed, sFMindex nfmindx, vector<ms_candidate> & vcand)
{
	uint32_t cnt = 0;
	struct ms_candidate candtmp;
	vector<ms_seed>::const_iterator ite = vseed.begin();
	for(; ite != vseed.end(); ++ite)
	{
		candtmp.seed_id = cnt++;		//seed_id从0开始
		for(uint32_t i = ite->sa_ary[0]; i <= ite->sa_ary[1]; ++i)
		{
			candtmp.ref_pos = calc_SA(nfmindx, i);
			vcand.push_back(candtmp);
		}
	}
}

void printvec_seed(const vector<ms_seed> & vseed)
{
	vector<ms_seed>::const_iterator ite = vseed.begin();
	for(; ite != vseed.end(); ++ite)
	{
		cout << ite->seed_segment << ":" << ite->start_pos << "~" << ite->end_pos << " sa:"<< ite->sa_ary[1] - ite->sa_ary[0] << endl;
	}
}

void printvec_candidate(const vector<ms_candidate> & vcand)
{
	vector<ms_candidate>::const_iterator ite = vcand.begin();
	for(; ite != vcand.end(); ++ite)
	{
		cout << ite->seed_id << "~" << ite->ref_pos << endl;
	}
}

void printvec_result(const vector<ms_result> & vrslt)
{
	vector<ms_result>::const_iterator ite = vrslt.begin();
	for(; ite != vrslt.end(); ++ite)
	{
		cout << ite->ref_start << "~" << ite->ref_end << " ed:"<< ite->ed << endl;
	}
}

void mapper()
{
}
