/*
 * load_DBG_full.h
 *
 *  Created on: Oct 9, 2019
 *      Author: bio
 */

#ifndef LOAD_DBG_FULL_H_
#define LOAD_DBG_FULL_H_

#include <stdint.h>

struct para_merge
{
	uint64_t *p_start_kmer;
	uint64_t *p_des_kmer;
	uint32_t *p_start_id;
	uint32_t *p_des_id;
	uint32_t blocksize;
	uint32_t thread_num;
	uint32_t thread_id;
	uint64_t ukN;
	uint32_t unitperkmer;
	uint32_t N_block;

};

struct para_dBGindex
{
	uint64_t bkN;
	uint64_t *p_branchedkmer;
	uint8_t *p_branchedkmerad;
	uint32_t uN;
	char **upath_arr;
	uint32_t upathid_start0;
	uint64_t ukN;
	uint64_t *p_unbranchedkmer;
	uint32_t *p_unbranchedkmerid;
	uint64_t **p2_bkmer;
	uint64_t **p2_ukmer;
};

void Save_ukmer(struct dBG * p_dBG, char * p_dbg_path);
void Divid_umers(struct dBG * p_dBG,char * p_dbg_path);
void SortA_umers(struct dBG * p_dBG,char * p_dbg_path,uint32_t thread_num);
void SortFile_umers(struct dBG * p_dBG,char * p_file_path,uint32_t thread_num);
void Gen_navigatSeq(struct dBG * p_dBG, char * p_dbg_path);
void gen_dBG_index(struct bit256KmerPara bit_para, struct para_dBGindex * p_sdBGidx, char * p_dbg_path, uint32_t thread_num);
void free_dBGindex(struct para_dBGindex * p_sdBGidx);


#endif /* LOAD_DBG_FULL_H_ */
