/*
 * verification.cpp
 *
 *  Created on: May 25, 2020
 *      Author: pluto
 */

#include "verification.h"


void verification(const vector<ms_seed> & vseed, const vector<ms_candidate> & vcand, vector<ms_result> & vrslt,
				const char *read, const char *ref, uint32_t tau)
{
	vector<ms_seed>::const_iterator ites = --vseed.end();
	vector<ms_candidate>::const_iterator itec = vcand.begin();
	uint32_t read_len = ites->end_pos + 1;
	struct ms_result rsltmp;
	char *readtmp = const_cast<char *>(read);
	char *reftmp = const_cast<char *>(ref);
	char *readcur;
	uint32_t seed_len;
	uint32_t linel,liner;
	for(;itec != vcand.end(); ++itec)
	{
		uint32_t tautmp = tau;
		readcur = readtmp + vseed[itec->seed_id].start_pos;
		seed_len = strlen(vseed[itec->seed_id].seed_segment);
		uint32_t remain_len = read_len - seed_len;
		if(!itec->seed_id)	//first seed of read and expand to right
		{
			readcur += seed_len;
			liner = min_lastlineEd(readcur, reftmp+itec->ref_pos+seed_len, remain_len, remain_len+tautmp, tautmp);
			if(liner != -1)
			{
				rsltmp.ref_start = itec->ref_pos;
				rsltmp.ref_end = rsltmp.ref_start + seed_len - 1 + liner;
				rsltmp.ed = tau - tautmp;
				vrslt.push_back(rsltmp);
			}
		}
		else if(itec->seed_id == vseed.size()-1)	//last seed of read and expand to left
		{
			readcur = readtmp;
			char *read_r = new char[remain_len+1]();
			strncpy(read_r,readcur,remain_len);
			reverseq(read_r);
			char *reftmp_r = new char[remain_len+tautmp+1]();
			strncpy(reftmp_r,reftmp+itec->ref_pos-remain_len-tautmp,remain_len+tautmp);
			reverseq(reftmp_r);
			linel = min_lastlineEd(read_r, reftmp_r, remain_len, remain_len+tautmp, tautmp);
			if(linel != -1)
			{
				rsltmp.ref_start = itec->ref_pos - linel;   //bug
				rsltmp.ref_end = itec->ref_pos + seed_len - 1;
				rsltmp.ed = tau - tautmp;
				vrslt.push_back(rsltmp);
			}
			delete [] read_r;
			delete [] reftmp_r;
		}
		else	//middle seed of read and expand to bidirection
		{
			remain_len -= vseed[itec->seed_id].start_pos;
			readcur += seed_len;
			liner = min_lastlineEd(readcur, reftmp+itec->ref_pos+seed_len, remain_len, remain_len+tautmp, tautmp);
			if(liner != -1)
			{
				readcur = readtmp;
//				rextlen -= liner;
				remain_len = vseed[itec->seed_id].start_pos;
				char *read_r = new char[remain_len+1]();
				strncpy(read_r,readcur,remain_len);
				reverseq(read_r);
				char *reftmp_r = new char[remain_len+tautmp+1]();
				strncpy(reftmp_r,reftmp+itec->ref_pos-remain_len-tautmp,remain_len+tautmp);
				reverseq(reftmp_r);
				linel = min_lastlineEd(read_r, reftmp_r, remain_len, remain_len+tautmp, tautmp);
				if(linel != -1)
				{
					rsltmp.ref_start = itec->ref_pos - linel; //bug
					rsltmp.ref_end = itec->ref_pos + seed_len-1 + liner;
					rsltmp.ed = tau - tautmp;
					vrslt.push_back(rsltmp);
				}
				delete [] read_r;
				delete [] reftmp_r;
			}
		}
	}

}
