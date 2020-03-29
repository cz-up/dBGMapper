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
#include "load_DBG_full.h"
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

void ext_treenode(struct bit256KmerPara bit_para, struct TPTnode *pnode, struct para_dBGindex sdBGidx,\
		uint32_t extlen, char *seq, uint32_t unitperkmer)
{
	if(extlen != 0)
	{
		if(pnode->offset != 0 && pnode->offset < strlen(seq) - bit_para.kmer1Len/2) //offset > 1?
		{
			if(pnode->dir == 'I')
			{
				pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
				pnode->p_child[0]->dir = pnode->dir;
				pnode->p_child[0]->p_parent = pnode;
				pnode->p_child[0]->offset = pnode->offset - 1;
				pnode->p_child[0]->c = seq[pnode->p_child[0]->offset];  //这里可能有问题
				for(uint32_t ii = 0 ; ii < 4; ii++)
				{
					pnode->p_child[0]->p_child[ii] = NULL;
				}
			}
			else
			{
				pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
				pnode->p_child[0]->dir = pnode->dir;
				pnode->p_child[0]->p_parent = pnode;
				pnode->p_child[0]->offset = pnode->offset + 1;
				pnode->p_child[0]->c = seq[pnode->offset+bit_para.kmer1Len/2];
				for(uint32_t ii = 0 ; ii < 4; ii++)
				{
					pnode->p_child[0]->p_child[ii] = NULL;
				}
			}
			ext_treenode(bit_para, pnode->p_child[0], sdBGidx, extlen-1, seq, unitperkmer);
		}
		else
		{
			uint64_t * const hashvalue_tmp = (uint64_t *)malloc(sizeof(uint64_t) * unitperkmer);
			uint32_t seqlen = strlen(seq);
			char *seqe = (char *)malloc(sizeof(char) * (seqlen+1));
			memset(seqe,0,seqlen+1);
			if(pnode->c == 'R')
			{
				strcpy(seqe,seq);
			}
			else
			{
				if(pnode->dir == 'O')
				{
					strncpy(seqe,seq+1,seqlen-1);
					seqe[seqlen-1] = pnode->c;
				}
				else
				{
					seqe[0] = pnode->c;
					strncpy(seqe+1,seq,seqlen-1);
				}
			}
			cal_hash_value_directly_256bit(seqe,hashvalue_tmp,bit_para);
			uint64_t **bkmer_ptr;
			bkmer_ptr = Tgenerate_array<uint64_t>(sdBGidx.bkN);
			uint64_t **ukmer_ptr;
			ukmer_ptr = Tgenerate_array<uint64_t>(sdBGidx.ukN);
			uint64_t bkmerfindret,ukmerfindret;
			bkmerfindret = Tfind_arrindexN<uint64_t>(bkmer_ptr, sdBGidx.p_branchedkmer, hashvalue_tmp,unitperkmer);
			ukmerfindret = Tfind_arrindexN<uint64_t>(ukmer_ptr, sdBGidx.p_unbranchedkmer, hashvalue_tmp,unitperkmer);
			if(bkmerfindret != ULLONG_MAX)
			{
				uint8_t adinfo = sdBGidx.p_branchedkmerad[bkmerfindret];
				if(pnode->dir == 'I')
				{
					adinfo >>= 4;
				}
				for(uint32_t i = 0; i < 4; i++)
				{
					if(adinfo & 0x01)
					{
						pnode->p_child[i] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
						pnode->p_child[i]->dir = pnode->dir;
						pnode->p_child[i]->p_parent = pnode;
						pnode->p_child[i]->offset = 0;
						for(uint32_t ii = 0 ; ii < 4; ii++)
						{
							pnode->p_child[i]->p_child[ii] = NULL;
						}
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
						ext_treenode(bit_para, pnode->p_child[i], sdBGidx, extlen-1, seqe, unitperkmer);
					}
					adinfo >>= 1;
				}
			}
			if(ukmerfindret != ULLONG_MAX)//如果是unipath上的kmer  判断在unipath上的offset 根据是向出度方向还是入度方向扩展来找
			{
				uint32_t ukmerid = sdBGidx.p_unbranchedkmerid[ukmerfindret];
				char *pos = NULL;
				pos = strstr(sdBGidx.upath_arr[ukmerid],seqe);
				pnode->offset = pos - sdBGidx.upath_arr[ukmerid];
				char ch;
				if(pnode->dir == 'I')
				{
					ch = *(pos-1);
					pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
					pnode->p_child[0]->dir = pnode->dir;
					pnode->p_child[0]->p_parent = pnode;
					pnode->p_child[0]->offset = pnode->offset - 1;
					pnode->p_child[0]->c = ch;
					for(uint32_t ii = 0 ; ii < 4; ii++)
					{
						pnode->p_child[0]->p_child[ii] = NULL;
					}
				}
				else
				{
					ch = *(pos + strlen(seqe));
					pnode->p_child[0] = (struct TPTnode *)malloc(sizeof(struct TPTnode));
					pnode->p_child[0]->dir = pnode->dir;
					pnode->p_child[0]->p_parent = pnode;
					pnode->p_child[0]->offset = pnode->offset + 1;
					pnode->p_child[0]->c = ch;
					for(uint32_t ii = 0 ; ii < 4; ii++)
					{
						pnode->p_child[0]->p_child[ii] = NULL;
					}
				}
				ext_treenode(bit_para, pnode->p_child[0], sdBGidx, extlen-1, sdBGidx.upath_arr[ukmerid], unitperkmer);
			}
			free(seqe);
			free(hashvalue_tmp);
		}
	}
}

void count_extnum(struct TPTnode node, struct seedext *p_seedext)
{
	if(node.p_child[0] == NULL && node.p_child[1] == NULL && node.p_child[2] == NULL && node.p_child[3] == NULL)
	{
		p_seedext->num++;
	}
	else
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			if(node.p_child[i] != NULL)
			{
				count_extnum(*node.p_child[i], p_seedext);
			}
		}
	}
}

void init_seedext(struct TPTnode node, struct seedext *p_seedext, char *seed)
{
	p_seedext->num = 0;
	count_extnum(node, p_seedext);
	p_seedext->dir = node.dir;
	p_seedext->seqext = (char**)malloc(sizeof(char*)*p_seedext->num);
	p_seedext->num = 0;
	p_seedext->seed = (char*)malloc(sizeof(char)*(strlen(seed)+1));
	memset(p_seedext->seed,0,strlen(seed)+1);
	strcpy(p_seedext->seed,seed);
}

void calc_seedextpara(struct seedext *p_seedext, char *calcseq, uint32_t tau, sFMindex n_index, sFMindex r_index, struct build_para para)  //如果dir为I 需不需要将calcseq reverse
{
	if(calcseq != NULL && tau > 0)
	{
		uint32_t seedlen = strlen(p_seedext->seed);
		uint32_t seqlen = strlen(calcseq);
		p_seedext->p3_extedarr = (uint32_t ***)malloc(sizeof(uint32_t **) * p_seedext->num);
		p_seedext->p2_extchcnt = (uint32_t **)malloc(sizeof(uint32_t*) * p_seedext->num);
		uint32_t ** edmatrix;
		for(uint32_t i = 0; i < p_seedext->num; i++)
		{
			uint32_t seedextlen = strlen(p_seedext->seqext[i]);
			edmatrix = originalEd(calcseq,p_seedext->seqext[i],seqlen,seedextlen);
			p_seedext->p3_extedarr[i] = (uint32_t **)malloc(sizeof(uint32_t *) * seedextlen);
			char *seqtmp = (char *)malloc(sizeof(char) * (seedlen+seedextlen+1));
			memset(seqtmp,0,seedlen+seedextlen+1);
			char *p_seqpos = seqtmp;
			p_seqpos += seedextlen;
			if(p_seedext->dir == 'I')
			{
				strcpy(p_seqpos,p_seedext->seed);
			}
			if(p_seedext->dir == 'O')
			{
				char *seqrev = (char *)malloc(sizeof(char) * (seedlen+1));
				strcpy(seqrev,p_seedext->seed);
				reverseq(seqrev);
				strcpy(p_seqpos,seqrev);
				if(seqrev)
				{
					free(seqrev);
					seqrev = NULL;
				}
			}
			p_seqpos--;
			p_seedext->p2_extchcnt[i] = (uint32_t *)malloc(sizeof(uint32_t) * seedextlen);
			for(uint32_t j = 0; j < seedextlen; j++)
			{
				*p_seqpos = p_seedext->seqext[i][j];
				cout << p_seqpos << endl;
				uint32_t *p_sarry = NULL;
				if(p_seedext->dir == 'I')
				{
					p_sarry = calc_SArangeSeq(n_index, para, p_seqpos); //暂时将build_para设为相同  如果需要设置不同的build_para再调整
				}
				if(p_seedext->dir == 'O')
				{
					p_sarry = calc_SArangeSeq(r_index, para, p_seqpos);
				}
				p_seedext->p2_extchcnt[i][j] = p_sarry[1] - p_sarry[0] + 1;
				p_seqpos--;
				if(p_sarry)
				{
					free(p_sarry);
					p_sarry = NULL;
				}
				p_seedext->p3_extedarr[i][j] = (uint32_t *)malloc(sizeof(uint32_t) * (2*tau + 1));
				for(uint32_t k = 0; k < 2*tau+1; k++)
				{
					if(k+j >= seqlen)
					{
						p_seedext->p3_extedarr[i][j][k] = edmatrix[seqlen][j+1];
					}
					else
					{
						p_seedext->p3_extedarr[i][j][k] = edmatrix[k+j][j+1];
					}
				}
			}
			free_edmatrix(&edmatrix, seqlen);
			if(seqtmp)
			{
				free(seqtmp);
				seqtmp = NULL;
			}
		}
	}
}

void free_seedext(struct seedext *p_seedext)
{
	if(p_seedext)
	{
		if(p_seedext->seed)
		{
			free(p_seedext->seed);
			p_seedext->seed = NULL;
		}
		if(p_seedext->p2_extchcnt)
		{
			for(uint32_t i = 0; i < p_seedext->num; i++)
			{
				free(p_seedext->p2_extchcnt[i]);
				p_seedext->p2_extchcnt[i] = NULL;
			}
			free(p_seedext->p2_extchcnt);
		}
		if(p_seedext->p3_extedarr)
		{
			for(uint32_t i = 0; i < p_seedext->num; i++)
			{
				uint32_t seedlen = strlen(p_seedext->seqext[i]);
				for(uint32_t j = 0; j < seedlen; j++)
				{
					free(p_seedext->p3_extedarr[i][j]);
					p_seedext->p3_extedarr[i][j] = NULL;
				}
				free(p_seedext->p3_extedarr[i]);
				p_seedext->p3_extedarr[i] = NULL;
			}
			free(p_seedext->p3_extedarr);
			p_seedext->p3_extedarr = NULL;
		}
		if(p_seedext->seqext)
		{
			for(uint32_t i = 0; i < p_seedext->num; i++)
			{
				free(p_seedext->seqext[i]);
				p_seedext->seqext[i] = NULL;

			}
			p_seedext->num = 0;
			free(p_seedext->seqext);
			p_seedext->seqext = NULL;
		}
	}
}

void print_extree(struct TPTnode node,char *seq, struct seedext *p_seedext)  //如果第三个参数为NULL 只打印extree
{
	if(node.p_child[0] == NULL && node.p_child[1] == NULL && node.p_child[2] == NULL && node.p_child[3] == NULL)
	{
		seq[strlen(seq)] = node.c;
		printf("%s",seq);
		if(p_seedext != NULL)
		{
			p_seedext->seqext[p_seedext->num] = (char*)malloc(sizeof(char)*strlen(seq));
			memset(p_seedext->seqext[p_seedext->num],0,strlen(seq));
			strncpy(p_seedext->seqext[p_seedext->num++],seq+1,strlen(seq));
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
				print_extree(*node.p_child[i],seq,p_seedext);
				seq[strlen(seq)-1] = '\0';
			}
		}
	}
}

void destory_extree(struct TPTnode *pnode)
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
}

