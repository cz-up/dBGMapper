/*
 * candidatePath.h
 *
 *  Created on: Jun 19, 2020
 *      Author: pluto
 */

#ifndef CANDIDATEPATH_H_
#define CANDIDATEPATH_H_
#include "basic.h"

struct candiPath{
	char *p_candidatepath;
	uint32_t start_seed_id;
	uint32_t end_seed_id;
	uint32_t *sa_array;

};



#endif /* CANDIDATEPATH_H_ */
