/*
 * seeding.cpp
 *
 *  Created on: Sep 18, 2019
 *      Author: bio
 */

//不同的seeding的方法
//seeding的操作是将原read序列拆成包含tau+1个seeds的划分
//目标是拆分出的seed的目标
//1）不能精确匹配的seed数量越多越好
//2）精确匹配的seed的长度越长越好

#include "seeding.h"
#include "Hash.h"
#include "BPlusTree_full.h"
#include "basic.h"
#include "method.h"

void reverse(char *a, uint32_t l, char ** ra)
{
	//求一个字符串的反串。
	char * p_tmp;
	p_tmp=(char*)malloc(sizeof(char)*l);
	for(uint32_t i=0;i<l;i++)
	{
		p_tmp[i]=a[l-1-i];
	}
	*ra=p_tmp;
}
void generate_seed_array(uint32_t * p_array_seed_len, struct seed * p_seed, uint32_t tau, uint32_t read_len)
{
	uint32_t start;
	start=0;
	p_seed->p_seed=(struct seed_single*)malloc(sizeof(struct seed_single)*(tau+1));
	p_seed->p_seed_reverse=(struct seed_single*)malloc(sizeof(struct seed_single)*(tau+1));
	for(uint32_t i=0;i<tau+1;i++)
	{
		p_seed->p_seed[i].len=p_array_seed_len[i];
		p_seed->p_seed[i].start=start;
		p_seed->p_seed[i].end=start+p_array_seed_len[i]-1;
		start=start+p_array_seed_len[i];
	}

	for(uint32_t i=0;i<tau+1;i++)
	{
		p_seed->p_seed_reverse[i].len=p_seed->p_seed[i].len;
		p_seed->p_seed_reverse[i].start=read_len-1-p_seed->p_seed[i].end;
		p_seed->p_seed_reverse[i].end=read_len-1-p_seed->p_seed[i].start;
	}
}
void generate_read_kmer_present_vec(struct NodeBit** p_kmer_root,uint32_t ** kmer_present_label,char * p_read, uint32_t read_len,uint32_t kmer_len)
{
	struct bit256KmerPara bit_para;
	get_para(&bit_para,kmer_len);

	uint32_t * kmer_present_label_tmp;
	kmer_present_label_tmp=(uint32_t *)malloc(sizeof(uint32_t)*(read_len-kmer_len+1));
	uint64_t * p_kmer_tmp;
	p_kmer_tmp=(uint64_t *)malloc(sizeof(uint64_t)*4);

	cal_hash_value_directly_256bit(p_read,p_kmer_tmp,bit_para);
	kmer_present_label_tmp[0]=getHashFTableValue(p_kmer_root,p_kmer_tmp,bit_para);
	for(uint32_t i=1;i<read_len-kmer_len+1;i++)
	{
		cal_hash_value_indirectly_256bit(p_read+i,p_kmer_tmp,p_kmer_tmp,bit_para);
		kmer_present_label_tmp[i]=getHashFTableValue(p_kmer_root,p_kmer_tmp,bit_para);
	}

	*kmer_present_label=kmer_present_label_tmp;
	free(p_kmer_tmp);
}
void ave_seeding(char * p_read, uint32_t read_len, struct seed * p_seed, uint32_t tau,uint32_t kmer_len)
{
	//按照长度平均拆分
	char * p_read_reverse;
	reverse(p_read, read_len, &p_read_reverse);
	p_seed->p_read=p_read;
	p_seed->p_read_reverse=p_read_reverse;
	p_seed->seed_total_num=tau+1;

	//如果不能产生长度都大于kmer_len的kmer，那么就报告错误，并返回。
	if(read_len/(tau+1)<kmer_len)
	{
		cout << "error: from ave_seeding() method!" << endl;
		cout << "description: the kmer length is too large, ave_seeding() can't generate seeds with length larger than kmer_len." <<endl;
		return;
	}

	//产生seed的长度
	uint32_t * p_array_seed_len;
	p_array_seed_len=(uint32_t *)malloc(sizeof(uint32_t)*(tau+1));

	uint32_t div=read_len/(tau+1);
	uint32_t mod=read_len%(tau+1);

	for(uint32_t i=0;i<tau+1;i++)
	{
		p_array_seed_len[i]=div;
	}
	for(uint32_t i=0;i<mod;i++)
	{
		p_array_seed_len[i]++;
	}

	generate_seed_array(p_array_seed_len, p_seed, tau, read_len);

	free(p_array_seed_len);
}

void non_ave_seeding(char * p_read, uint32_t read_len, struct seed * p_seed, uint32_t tau,uint32_t kmer_len, uint32_t * kmer_present_label)
{
	//选择一组没有出现的kmer，如果将其作为分割边界，能够将整个read分割为tau+1个seeds，那么称这组没有出现的kmer是符合条件的。
	//选择一组没有出现的kmer，其中包含的kmer的数量最多。然后基于这组没有出现的kmer对read进行seeding。
}

void gen_candidateSet()
{

}


//  2020.3.3
void gen_truepathtree(uint64_t *kmer, uint32_t unitperkmer, struct para_dBGindex sdBGidx, uint32_t extlen)
{
	uint64_t **bkmer_ptr;
	bkmer_ptr = Tgenerate_array<uint64_t>(sdBGidx.bkN);
	uint64_t **ukmer_ptr;
	ukmer_ptr = Tgenerate_array<uint64_t>(sdBGidx.ukN);
	uint64_t bkmerfindret,ukmerfindret;
	bkmerfindret = Tfind_arrindexN<uint64_t>(bkmer_ptr, sdBGidx.p_branchedkmer, kmer,unitperkmer);
	ukmerfindret = Tfind_arrindexN<uint64_t>(ukmer_ptr, sdBGidx.p_unbranchedkmer, kmer,unitperkmer);
	uint64_t *pbkmer,*pukmer;
	pbkmer = sdBGidx.p_branchedkmer;
	pukmer = sdBGidx.p_unbranchedkmer;
	for(uint64_t i = 0; i < sdBGidx.bkN; i++)
	{

	}
//	if(bkmerfindret != -1)
//	{
//		struct TPTnode root;
//		root.c = 'R';
//		root.dir = 'R';
//		if(sdBGidx.p_branchedkmerad[bkmerfindret] != 0)
//		{
//
//		}
//	}
//	else
//	{
//
//	}

	Tfree_genarray<uint64_t>(&ukmer_ptr);
	Tfree_genarray<uint64_t>(&bkmer_ptr);
}

void calc_edarray(struct TPTnode *pnode, char *seq, uint8_t tau)
{
	pnode->edarry = (uint32_t *)malloc(sizeof(uint32_t)*(2*tau+1));
	for(uint32_t i = 0; i < 2*tau+1; i++)
	{
		pnode->edarry[i] = tau+1;
	}
	uint32_t flag,flagr0;
	flagr0 = 1;
	int upper = pnode->level - tau;
	int seqlen = strlen(seq);
	int start = 0;
	if(upper > 1)
	{
		start = upper - 1;
	}
	uint32_t left, top ,leftop;
	for(int32_t i = 0; i < 2*tau+1; i++)
	{
		flag = 1;
		if(upper + i <= 0)
		{
//			flagr0 = 0;
			continue;
		}
		if(upper + i == seqlen + 1)
		{
			break;
		}
		if(pnode->c == seq[start++])
		{
			flag = 0;
		}
		left = pnode->p_parent->edarry[i+1] + 1;
		top = pnode->edarry[i-1] + 1;
		leftop = pnode->p_parent->edarry[i] + flag;
		if(!i)
		{
			pnode->edarry[i] = min(leftop,left);
			continue;
		}
		if(i == 2*tau)
		{
			pnode->edarry[i] = min(leftop,top);
			break;
		}
		pnode->edarry[i] = min(leftop,top);
		pnode->edarry[i] = min(left,pnode->edarry[i]);
	}
//	for(int i = 0; i < 2*tau+1; i++)
//	{
//		cout << pnode->edarry[i] << " ";
//	}
//	cout << endl;
}

void init_rootnode(struct TPTnode *pnode,struct seed_extpara ext_para, char *seq)
{
	pnode->extdone = false;
	pnode->c = 'R';
	pnode->level = 0;
	pnode->offset = 0;
	pnode->edarry = (uint32_t*)malloc(sizeof(uint32_t)*(2*ext_para.tau+1));
	uint32_t i = 0;
	pnode->p_parent = NULL;
	for(i = 0 ; i < 4; i++)
	{
		pnode->p_child[i] = NULL;
	}
	i = 0;
	for(; i < ext_para.tau+1; i++)
	{
		pnode->edarry[i] = 0;
	}
	int tmp = 0;
	for(; i < 2*ext_para.tau+1; i++)
	{
		pnode->edarry[i] = ++tmp;
	}
	char *tmpseq = (char *)malloc(sizeof(char)*(strlen(seq)+1));
	memset(tmpseq,0,strlen(seq)+1);
	strcpy(tmpseq,seq);
	if(ext_para.dir == 'O')
	{
		reverseq(tmpseq);
	}
	pnode->saarry = calc_SArangeSeq(ext_para.FMidx,tmpseq);
	free(tmpseq);
}

bool init_childnode(struct seed_extpara ext_para, struct TPTnode *pnodec, struct TPTnode *pnodep)  //返回真表示继续扩展
{
	pnodec->p_parent = pnodep;
	pnodec->level = pnodep->level + 1;
	calc_edarray(pnodec, ext_para.alignseq, ext_para.tau);
	pnodec->saarry = calc_SArangeChar(ext_para.FMidx,pnodep->saarry,pnodec->c);
	if(pnodec->saarry[1] >= pnodec->saarry[0])
	{
		for(uint32_t i = 0; i < 2*ext_para.tau+1; i++)
		{
			if(pnodec->edarry[i] <= ext_para.tau)
			{
				for(uint32_t ii = 0 ; ii < 4; ii++)
				{
					pnodec->p_child[ii] = NULL;
				}
				return true;
			}
		}
	}
	return false;
}

void ext_treenode(struct seed_extpara ref_para, struct TPTnode *pnode,\
		char *seq, uint32_t extlen, bool *ponunipath)
{
	if(extlen != 0)
	{
		bool extflag = false;
		if(pnode->offset > 0 && pnode->offset < strlen(seq) - ref_para.bit_para.kmer1Len/2) //offset > 1?
		{
			if(ref_para.dir == 'I')
			{
				pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
				pnode->p_child[0]->offset = pnode->offset - 1;
				pnode->p_child[0]->c = seq[pnode->p_child[0]->offset];  //这里可能有问题
			}
			else
			{
				pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
				pnode->p_child[0]->offset = pnode->offset + 1;
				pnode->p_child[0]->c = seq[pnode->offset+ref_para.bit_para.kmer1Len/2];
			}
			extflag = init_childnode(ref_para, pnode->p_child[0], pnode);
			if(extflag == true)
			{
				*ponunipath = true;
				ext_treenode(ref_para, pnode->p_child[0], seq, extlen-1, ponunipath);
			}
			else
			{
				free(pnode->p_child[0]);
				pnode->p_child[0] = NULL;
			}
		}
		else
		{
			uint64_t * const hashvalue_tmp = (uint64_t *)malloc(sizeof(uint64_t) * ref_para.bit_para.kmer64Len);
			uint32_t seqlen = ref_para.bit_para.kmer1Len/2;
			char *seqe = (char *)malloc(sizeof(char) * (seqlen+1));
			memset(seqe,0,seqlen+1);
			if(pnode->c == 'R')
			{
				strcpy(seqe,seq);
			}
			else
			{
				if(true == *ponunipath)
				{
					strncpy(seqe,seq+pnode->offset,seqlen);
				}
				else
				{
					if(ref_para.dir == 'I')
					{
						seqe[0] = pnode->c;
						strncpy(seqe+1,seq,seqlen-1);
					}
					else
					{
						strncpy(seqe,seq+1,seqlen-1);
						seqe[seqlen-1] = pnode->c;
					}
				}
			}
			cal_hash_value_directly_256bit(seqe,hashvalue_tmp,ref_para.bit_para);
			uint64_t bkmerfindret,ukmerfindret;
			bkmerfindret = Tfind_arrindexN<uint64_t>(ref_para.sdBGidx.p2_bkmer, ref_para.sdBGidx.p_branchedkmer, \
					hashvalue_tmp,ref_para.bit_para.kmer64Len);
			if(bkmerfindret != ULLONG_MAX)
			{
				uint8_t adinfo = ref_para.sdBGidx.p_branchedkmerad[bkmerfindret];
				if(ref_para.dir == 'I')
				{
					adinfo >>= 4;
				}
				for(uint32_t i = 0; i < 4; i++)
				{
					if(adinfo & 0x01)
					{
						pnode->p_child[i] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
						pnode->p_child[i]->offset = 0;
						switch(i)
						{
							case 0 :
								pnode->p_child[i]->c = 'A';break;
							case 1 :
								pnode->p_child[i]->c = 'C';break;
							case 2 :
								pnode->p_child[i]->c = 'G';break;
							case 3 :
								pnode->p_child[i]->c = 'T';break;
							default:
								pnode->p_child[i]->c = 'E';
						}
						extflag = init_childnode(ref_para, pnode->p_child[i], pnode);
						if(extflag == true)
						{
							*ponunipath = false;
							ext_treenode(ref_para, pnode->p_child[i], seqe, extlen-1, ponunipath);
						}
						else
						{
							free(pnode->p_child[i]);
							pnode->p_child[i] = NULL;
						}
					}
					adinfo >>= 1;
				}
			}
			ukmerfindret = Tfind_arrindexN<uint64_t>(ref_para.sdBGidx.p2_ukmer, ref_para.sdBGidx.p_unbranchedkmer, \
					hashvalue_tmp,ref_para.bit_para.kmer64Len);
			if(ukmerfindret != ULLONG_MAX)//如果是unipath上的kmer  判断在unipath上的offset 根据是向出度方向还是入度方向扩展来找
			{
				uint32_t ukmerid = ref_para.sdBGidx.p_unbranchedkmerid[ukmerfindret];
				char *pos = NULL;
				pos = strstr(ref_para.sdBGidx.upath_arr[ukmerid],seqe);
				pnode->offset = pos - ref_para.sdBGidx.upath_arr[ukmerid];
				char ch;
				if(ref_para.dir == 'I')
				{
					ch = *(pos-1);
					pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
					pnode->p_child[0]->offset = pnode->offset - 1;
					pnode->p_child[0]->c = ch;
				}
				else
				{
					ch = *(pos + strlen(seqe));
					pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
					pnode->p_child[0]->offset = pnode->offset + 1;
					pnode->p_child[0]->c = ch;
				}
				extflag = init_childnode(ref_para, pnode->p_child[0], pnode);
				if(extflag == true)
				{
					*ponunipath = true;
					ext_treenode(ref_para, pnode->p_child[0], ref_para.sdBGidx.upath_arr[ukmerid], extlen-1, ponunipath);

				}
				else
				{
					free(pnode->p_child[0]);
					pnode->p_child[0] = NULL;
				}
			}
			free(seqe);
			free(hashvalue_tmp);
		}
	}
}

void print_extree(struct TPTnode node,char *seq)
{
	if(node.p_child[0] == NULL && node.p_child[1] == NULL && node.p_child[2] == NULL && node.p_child[3] == NULL)
	{
		seq[strlen(seq)] = node.c;
		printf("%s\n",seq);
		seq[strlen(seq)-1] = '\0';
	}
	else
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			if(node.p_child[i] != NULL)
			{
				seq[strlen(seq)] = node.c;
				print_extree(*node.p_child[i],seq);
				seq[strlen(seq)-1] = '\0';
			}
		}
	}
}

void print_specificlen(struct TPTnode node,struct seed_extpara ext_set, uint32_t extlen, char *seq)
{
	if(node.p_child[0] == NULL && node.p_child[1] == NULL && node.p_child[2] == NULL && node.p_child[3] == NULL)
	{
		seq[strlen(seq)] = node.c;
		if(node.level >= extlen - ext_set.tau && node.level <= extlen + ext_set.tau)
		{
			printf("%s\n",seq);
		}
		seq[strlen(seq)-1] = '\0';
	}
	else
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			if(node.p_child[i] != NULL)
			{
				seq[strlen(seq)] = node.c;
				print_specificlen(*node.p_child[i], ext_set, extlen, seq);
				seq[strlen(seq)-1] = '\0';
			}
		}
	}
}

void destory_extree(struct TPTnode *pnode)
{
	if(pnode)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			if(pnode->p_child[i] != NULL)
			{
				destory_extree(pnode->p_child[i]);
				free(pnode->p_child[i]);
				pnode->p_child[i] = NULL;
			}
		}
		if(pnode->saarry)
		{
			free(pnode->saarry);
			pnode->saarry = NULL;
		}
		if(pnode->edarry)
		{
			free(pnode->edarry);
			pnode->edarry = NULL;
		}
	}
}

