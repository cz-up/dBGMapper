/*
 * mapperStruct.h
 *
 *  Created on: May 9, 2020
 *      Author: pluto
 */

#ifndef MAPPERSTRUCT_H_
#define MAPPERSTRUCT_H_

#include "basic.h"
#include "FMindex_ExactMatch.h"

struct ms_seed{
	char *seed_segment;
	uint32_t start_pos;
	uint32_t end_pos;
	uint32_t *sa_ary;
};

struct candiPath{
	char *p_candidatepath;
	uint32_t start_seed_id;
	uint32_t end_seed_id;
	uint32_t *sa_ary;

};

struct PH_Node{
	uint32_t start_seed_id;
	uint32_t end_seed_id;
	uint32_t tau;
};

struct ms_candidate
{
	uint32_t seed_id;
	uint32_t ref_pos;
};

struct ms_result
{
	uint32_t ref_start;
	uint32_t ref_end;
	uint32_t ed;
};

void generate_PHNArray(struct PH_Node **p2_node, uint32_t tau);
uint32_t PHNode_index(struct PH_Node *p_nodeA, struct PH_Node NodeT);

void generate_seeds(const char *read, uint32_t tau, sFMindex nfmindx, vector<ms_seed> & vseed);
void free_seed(vector<ms_seed> & vseed);
void generate_candidate(const vector<ms_seed> & vseed, sFMindex nfmindx, vector<ms_candidate> & vcand);
void verification(const vector<ms_seed> & vseed, const vector<ms_candidate> & vcand, vector<ms_result> & vrslt,
				const char *read, const char *ref, uint32_t tau);
void printvec_seed(const vector<ms_seed> & vseed);
void printvec_candidate(const vector<ms_candidate> & vcand);
void printvec_result(const vector<ms_result> & vrslt);

#endif /* MAPPERSTRUCT_H_ */
