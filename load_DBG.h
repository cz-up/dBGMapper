/*
 * load_DBG.h
 *
 *  Created on: Sep 16, 2019
 *      Author: bio
 */

#ifndef LOAD_DBG_H_
#define LOAD_DBG_H_
#include "basic.h"

struct kmer_detail
{
	uint8_t is_branched;
	uint8_t ad;
	uint32_t unipath_id;
	uint16_t unipath_offset;
};

struct dBG
{
	struct NodeBit** p_kmer_root;
	struct kmer_detail *p_kmer_detail;
	char ** p_unipath;
	uint32_t kN;
	uint32_t uN;
	uint32_t L;

};
struct unipath
{
	uint32_t start;
	uint32_t len;
	uint8_t ref_id;
};

uint32_t get_Dbg_file_name(char* p_dbg_path, char *** p_dbg_file);
uint64_t get_total_kmer_number(char ** p_dbg_file, uint64_t *kN, uint64_t *ukN);
uint32_t get_total_unipath_number(char ** p_dbg_file);
void get_unipath_struct(uint64_t a,struct unipath *b);
void get_unipath_kmer_ad(char* p_cur,uint32_t kmer_len,uint8_t * ad);
void save_kmer_details(struct kmer_detail * a,\
		uint8_t is_branched,\
		uint8_t ad,\
		uint32_t unipath_id,\
		uint16_t unipath_offset);
void loadDbg(struct dBG * p_dBG, char * p_dbg_path);
void insukmer2BplusTree(struct dBG * p_dBG, char * p_dbg_path, uint32_t adsel, uint32_t offsetsel);
#endif /* LOAD_DBG_H_ */
