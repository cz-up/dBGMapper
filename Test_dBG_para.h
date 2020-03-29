/*
 * Test_dBG_para.h
 *
 *  Created on: Dec 25, 2019
 *      Author: pluto
 */

#ifndef TEST_DBG_PARA_H_
#define TEST_DBG_PARA_H_

#include "basic.h"

void test_2MethodRate(struct dBG * p_dBG, char * p_file_path);
void test_bkmer_index(struct dBG * p_dBG, char * p_dbg_path);
void test_kmernum(struct dBG * p_dBG, char * p_dbg_path);
void test_lendis(struct dBG * p_dBG, char * p_dbg_path, int unidivnum);
void Gen_unipathLenInf(struct dBG * p_dBG, char * p_dbg_path);
void merge_superunipath(struct dBG * p_dBG, char * p_dbg_path,uint8_t ** p_uni_adinf,uint32_t **p_pre_arrlen);
void Test_dBG_Attribute(struct dBG * p_dBG, char * p_dbg_path);

#endif /* TEST_DBG_PARA_H_ */
