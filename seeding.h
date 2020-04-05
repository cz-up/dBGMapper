/*
 * seeding.h
 *
 *  Created on: Sep 18, 2019
 *      Author: bio
 */

#ifndef SEEDING_H_
#define SEEDING_H_

#include "basic.h"
#include "ExactMatch.h"

struct seed_single
{
	uint32_t start;
	uint32_t end;
	uint32_t len;
};

struct seed
{
	char * p_read;
	char * p_read_reverse;
	struct seed_single * p_seed;
	struct seed_single * p_seed_reverse;
	uint32_t seed_total_num;
};

struct TPTnode
{
	char c;
	uint8_t level;
	uint8_t offset;
	struct TPTnode * p_parent;
	struct TPTnode * p_child[4];
	uint32_t * edarry;
	uint32_t * saarry;
};

struct seedext
{
	char dir;
	uint32_t num;
	char *seed;
	char **seqext;
	uint32_t *** p3_extedarr;
	uint32_t ** p2_extchcnt;
};

void ext_treenode(struct bit256KmerPara bit_para, struct TPTnode *pnode, struct para_dBGindex sdBGidx, sFMindex FMidx,\
		char dir, uint32_t extlen, char *seq, char *alignseq, uint8_t tau);
void init_seedext(struct TPTnode node, struct seedext *p_seedext, char *seed);
void calc_seedextpara(struct seedext *p_seedext, char *calcseq, uint32_t tau, sFMindex n_index, sFMindex r_index);
void free_seedext(struct seedext *p_seedext);
void print_extree(struct TPTnode node,char *seq, struct seedext *p_seedext);
void destory_extree(struct TPTnode *pnode);

#endif /* SEEDING_H_ */
