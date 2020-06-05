/*
 * Test_dBG_para.h
 *
 *  Created on: Dec 25, 2019
 *      Author: pluto
 */

#ifndef ANALYSE_DBG_H_
#define ANALYSE_DBG_H_

#include "basic.h"
#include "Binary_Search.h"
#include "load_DBG_full.h"
#include "BplusTreeBit.h"

void test_2MethodRate(struct dBG * p_dBG, char * p_file_path);
void test_bkmer_index(struct dBG * p_dBG, char * p_dbg_path);
void test_kmernum(struct dBG * p_dBG, char * p_dbg_path);
void test_lendis(struct dBG * p_dBG, char * p_dbg_path, int unidivnum);
void Gen_unipathLenInf(struct dBG * p_dBG, char * p_dbg_path);
void merge_superunipath(struct dBG * p_dBG, char * p_dbg_path,uint8_t ** p_uni_adinf,uint32_t **p_pre_arrlen);
void Test_dBG_Attribute(struct dBG * p_dBG, char * p_dbg_path);

#endif /* ANALYSE_DBG_H_ */
