/*
 * SortedArray.h
 *
 *  Created on: Dec 5, 2019
 *      Author: pluto
 */

#ifndef SORTEDARRAY_H_
#define SORTEDARRAY_H_
#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <limits.h>

class SortedArray {
public:
	SortedArray();
	SortedArray(uint32_t n);
	void SetMemory(uint32_t n);
	virtual ~SortedArray();
	uint64_t * array;
	uint32_t * arrvalue;
	uint32_t arrlen;
	uint32_t index;
private:
};

#endif /* SORTEDARRAY_H_ */
