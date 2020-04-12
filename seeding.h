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
#include "load_DBG_full.h"

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

struct ref_extpara
{
	struct bit256KmerPara bit_para;
	struct para_dBGindex sdBGidx;
	struct sFMindex FMidx;
};

struct ext_setting
{
	char dir;
	char *alignseq;
	uint8_t tau;
	bool extdone;
};

struct TPTnode
{
	bool extdone;
	char c;
	uint8_t level;
	uint8_t offset;			//假设unipath最大长度不超过255   统计如果大于255  改为uint16_t
	struct TPTnode * p_parent;
	struct TPTnode * p_child[4];
	uint32_t * edarry;
	uint32_t * saarry;
};

void init_rootnode(struct TPTnode *pnode, char *seq, struct ext_setting ext_para, sFMindex *pFMinx, sFMindex nFMidx, sFMindex rFMidx);
void ext_treenode(struct ref_extpara ref_para, struct TPTnode *pnode,struct ext_setting ext_set, char *seq, uint32_t extlen);
void print_extree(struct TPTnode node,char *seq);
void destory_extree(struct TPTnode *pnode);

#endif /* SEEDING_H_ */
