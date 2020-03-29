/*
 * SortedArray.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: pluto
 */

#include "SortedArray.h"

SortedArray::SortedArray(){
	this->array = NULL;
	this->arrvalue = NULL;
	this->arrlen = 0;
	this->index = 0;
}

SortedArray::SortedArray(uint32_t n):arrlen(n), index(0){
	// TODO Auto-generated constructor stub
	this->array = new uint64_t[n];
	this->arrvalue = new uint32_t[n];
}

void SortedArray::SetMemory(uint32_t n){
	this->array = new uint64_t[n];
	this->arrvalue = new uint32_t[n];
}

SortedArray::~SortedArray() {
	// TODO Auto-generated destructor stub
	if(array){
		delete [] array;
		array = NULL;
		delete [] arrvalue;\
		arrvalue = NULL;
	}
}

void SortVector(std::vector<SortedArray*> ver, uint64_t *instkey, uint32_t *instvalue){
    int32_t len = 0;
    uint64_t min_val = 0;
    int32_t min_idx = -1;
    uint32_t pos = 0;
    //多个有序数组进行归并操作
    while((len = ver.size()) > 0){
        min_val = ULONG_MAX;
        min_idx = -1;
        //选取当前轮次的最小值
        for(uint32_t i = 0; i < len; ++i){
            if (ver[i]->index >= ver[i]->arrlen){
                ver.erase(ver.begin() + i);
                break;
            }

            uint64_t tmp_val = ver[i]->array[ver[i]->index];
            if (tmp_val <= min_val){
                min_val = tmp_val;
                min_idx = i;
            }
        }

        //打印
        if (min_idx != -1){
        	instkey[pos++] = *(ver[min_idx]->array);
        	instvalue[pos++] = *(ver[min_idx]->arrvalue);
            ++(ver[min_idx]->index);
        }
    }
}
