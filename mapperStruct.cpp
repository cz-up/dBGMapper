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
		seedtmp.sa_arry = calc_SArangeSeq(nfmindx,seedtmp.seed_segment);
		vseed.push_back(seedtmp);
		start += lentmp;
		readlen -= lentmp;
		--cnt;
	}

}

void printvec_seed(const vector<ms_seed> & vseed)
{
	vector<ms_seed>::const_iterator ite = vseed.begin();
	for(; ite != vseed.end(); ++ite)
	{
		cout << ite->seed_segment << ":" << ite->start_pos << "~" << ite->end_pos << " sa:"<< ite->sa_arry[1] - ite->sa_arry[0] << endl;
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
		for(uint32_t i = ite->sa_arry[0]; i <= ite->sa_arry[1]; ++i)
		{
			candtmp.ref_pos = calc_SA(nfmindx, i);
			vcand.push_back(candtmp);
		}
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
