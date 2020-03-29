/*
 * ExactMatch.h
 *
 *  Created on: 2019年7月13日
 *      Author: jupiter
 */

#ifndef EXACTMATCH_H_
#define EXACTMATCH_H_

#include "basic.h"
#include "read.h"
//#include "build.h"

char * get_queryfrag(uint32_t k,char *query);
uint32_t calc_C(sFMindex mem,char c);
uint32_t calc_OCC(sFMindex mem,struct build_para para,char c,uint32_t pos);
uint32_t LF_Mapping(sFMindex mem,struct build_para para,char c,uint32_t pos);
uint32_t* calc_SArangeSeq(sFMindex mem,struct build_para para,char *read);
uint32_t* calc_SArangeChar(sFMindex mem,struct build_para para,uint32_t *pre, char ch);
uint32_t calc_SA(sFMindex mem,struct build_para para,uint32_t pos);

#endif /* EXACTMATCH_H_ */
