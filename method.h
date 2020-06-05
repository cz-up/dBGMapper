/*
 * method.h
 *
 *  Created on: Feb 16, 2019
 *      Author: bio
 */

#ifndef METHOD_H_
#define METHOD_H_
#include "basic.h"
using namespace std;

uint32_t min_lastlineEd(char *a,char * b,int32_t x,int32_t y, uint32_t &tau);
uint32_t **originalEd(char *a,char * b,int32_t x,int32_t y);
void free_edmatrix(uint32_t ***m, int32_t x);
uint32_t rangeEd(char *a,char * b,int32_t x,int32_t y,int32_t ed);
uint32_t **levelEd(char *a,char * b,int32_t x,int32_t y,int32_t ed);
int recur_ed(string s,string r,int tor);

			


#endif /* METHOD_H_ */
