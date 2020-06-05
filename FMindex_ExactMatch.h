/*
 * ExactMatch.h
 *
 *  Created on: 2019年7月13日
 *      Author: jupiter
 */

#ifndef FMINDEX_EXACTMATCH_H_
#define FMINDEX_EXACTMATCH_H_

#include "basic.h"
#include "read.h"

char * get_queryfrag(uint32_t k,char *query);
uint32_t calc_C(sFMindex mem,char c);
uint32_t calc_OCC(sFMindex mem,char c,uint32_t pos);
uint32_t LF_Mapping(sFMindex mem,char c,uint32_t pos);
uint32_t* calc_SArangeSeq(sFMindex mem,char *read);
uint32_t* calc_SArangeChar(sFMindex mem,uint32_t *pre, char ch);
uint32_t calc_SA(sFMindex mem,uint32_t pos);

#endif /* FMINDEX_EXACTMATCH_H_ */
