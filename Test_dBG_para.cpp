/*
 * Test_dBG_para.cpp
 *
 *  Created on: Dec 25, 2019
 *      Author: pluto
 */

#include "Test_dBG_para.h"
#include "load_DBG.h"
#include "BplusTreeBit.h"
#include "BPlusTree_full.h"

void calc_onesinchar(uint8_t in,uint64_t &num)
{
	for(uint32_t i = 0; i < 8; i++)
	{
		num += (in & 1);
		in >>= 1;
	}

}

void test_2MethodRate(struct dBG * p_dBG, char * p_file_path)	//测试B+树的生成速度和查询速度都比**慢
{
	struct bit256KmerPara bit_para;
	get_para(&bit_para,p_dBG->L);
	uint32_t unitperkmer;
	unitperkmer = (uint32_t)ceil((double)p_dBG->L / 32);
	cout << "unitperkmer : " << unitperkmer << endl;

	uint64_t kN;
	FILE * fp_kmerpath;
	fp_kmerpath = fopen(p_file_path,"rb");
	fseek(fp_kmerpath,0,2);
	kN = ftell(fp_kmerpath)/(unitperkmer*sizeof(uint64_t));
	rewind(fp_kmerpath);

	uint64_t * const hashvalue_tmp = (uint64_t *)malloc(sizeof(uint64_t) * unitperkmer);
	uint64_t * kmer_array;
	uint64_t unitnum = kN * unitperkmer;
	kmer_array = (uint64_t *)malloc(sizeof(uint64_t)*unitnum);
	fread(kmer_array,sizeof(uint64_t),unitnum,fp_kmerpath);
	fclose(fp_kmerpath);
	struct NodeBit * p_root;
	uint64_t *p_kmer_array = kmer_array;
	struct nodeBit node_tmp;
	p_root = (struct NodeBit*)malloc(sizeof(struct NodeBit));
	Node_initial_bit(p_root);
	for(uint32_t i = 0; i < kN; i++)
	{
		node_tmp.hashValue = p_kmer_array;
		Insert_Value_bit(&p_root, node_tmp, bit_para);
		p_kmer_array += unitperkmer;
	}
	struct timeval tvs,tve;
	p_kmer_array = kmer_array;
	gettimeofday(&tvs,NULL);
	for(uint32_t i = 0; i < kN; i++)
	{
		MappingHashValueToID_bit(p_root,p_kmer_array,bit_para);
		p_kmer_array += unitperkmer;
	}
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout << "B+Tree Search time is: "<< span << endl;
	destory_tree_bit(p_root,bit_para);

	p_kmer_array = kmer_array;
	uint32_t **bkmer_ptr;
	bkmer_ptr = generate_array(kN, 0);
	gettimeofday(&tvs,NULL);
	for(uint32_t i = 0; i < kN; i++)
	{
		find_arrindex(bkmer_ptr, kmer_array, *p_kmer_array);   //为什么输入指向数组的指针会有段错误
		p_kmer_array += unitperkmer;
	}
	gettimeofday(&tve,NULL);
	span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout << "**array Search time is: "<< span << endl;
	free_genarray(&bkmer_ptr);
	free(kmer_array);
	free(hashvalue_tmp);
}

void test_bkmer_index(struct dBG * p_dBG, char * p_dbg_path)
{
	struct bit256KmerPara bit_para;
	get_para(&bit_para,p_dBG->L);

	//1)解析文件名：从p_dbg_path文件中读入4个文件名，分别为：
	//kmer026, kmer026ad, unipath026, nomLeu2.fa
	char ** p_dbg_file;
	uint32_t fN;
	fN=get_Dbg_file_name(p_dbg_path,&p_dbg_file);

	cout << fN << endl;
	cout << p_dbg_file[0] << endl;
	cout << p_dbg_file[1] << endl;
	cout << p_dbg_file[2] << endl;
	cout << p_dbg_file[3] << endl;

	//2)计算kmer的总数量
	uint64_t tN;
	uint64_t bkN,ukN;
	tN=get_total_kmer_number(p_dbg_file,&bkN,&ukN);
	cout << tN << endl;

	//3)装载分叉kmer
	uint64_t kbN;
	FILE* fp_branched_kmer_file;
	fp_branched_kmer_file=fopen(p_dbg_file[0],"rb");
	fseek(fp_branched_kmer_file,0,2);
	kbN=ftell(fp_branched_kmer_file)/(sizeof(uint64_t)*bit_para.kmer64Len);
	fseek(fp_branched_kmer_file,0,0);

	FILE* fp_branched_kmer_file_ad;
	fp_branched_kmer_file_ad=fopen(p_dbg_file[1],"rb");

	uint64_t * p_branched_kmer;
	p_branched_kmer=(uint64_t *)malloc(sizeof(uint64_t)*kbN*bit_para.kmer64Len);
	uint8_t * p_branched_kmer_ad;
	p_branched_kmer_ad=(uint8_t *)malloc(sizeof(uint8_t)*kbN);
	fread(p_branched_kmer,sizeof(uint64_t),kbN*bit_para.kmer64Len,fp_branched_kmer_file);
	fread(p_branched_kmer_ad,sizeof(uint8_t),kbN,fp_branched_kmer_file_ad);

	uint32_t ** p_branched_kmer_index;
	p_branched_kmer_index=generate_array(kbN,0);

	cout << "generate branched kmer index is over!" << endl;

	//4)测试新的索引是不是好用
	for(uint32_t i=0;i<kbN;i++)
	{
		if(find_arrindex(p_branched_kmer_index,\
				p_branched_kmer,\
				p_branched_kmer[i])!=i)
		{
			cout << "error: from branched kmer index!" << endl;
		}
	}
	fclose(fp_branched_kmer_file);
}

void test_kmernum(struct dBG * p_dBG, char * p_dbg_path)	//output unipath num & unbranched kmer num
{
	struct bit256KmerPara bit_para;
	get_para(&bit_para,p_dBG->L);

	//1)解析文件名：从p_dbg_path文件中读入4个文件名，分别为：
	//kmer026, kmer026ad, unipath026, nomLeu2.fa
	char ** p_dbg_file;
	uint32_t fN;
	fN=get_Dbg_file_name(p_dbg_path,&p_dbg_file);

	cout << fN << endl;
	cout << p_dbg_file[0] << endl;
	cout << p_dbg_file[1] << endl;
	cout << p_dbg_file[2] << endl;
	cout << p_dbg_file[3] << endl;
	cout << p_dbg_file[4] << endl;

	//5)计算unipath数量
	uint32_t uN;
	uN = get_total_unipath_number(p_dbg_file);
	cout <<"uN : "<< uN << endl;

	//7）读unipath
	FILE * fp_unipath;
	fp_unipath=fopen(p_dbg_file[2],"rb");
	uint64_t *p_unipath;
	p_unipath=(uint64_t*)malloc(sizeof(uint64_t)*uN);
	fread(p_unipath,sizeof(uint64_t),uN,fp_unipath);
	fclose(fp_unipath);

	FILE * fp_unipathcnt;
	fp_unipathcnt=fopen(p_dbg_file[4],"rb");
	uint64_t *p_unipathcnt;
	p_unipathcnt=(uint64_t*)malloc(sizeof(uint64_t)*uN);
	fread(p_unipathcnt,sizeof(uint64_t),uN,fp_unipathcnt);
	fclose(fp_unipathcnt);

	//8)计算super unipath总长度
	struct unipath cur_struct_unipath_tmp;
	uint64_t ukN = 0;
	uint32_t calc_tmp;
	for(uint32_t i=0;i<uN;i++)
	{
		get_unipath_struct(p_unipath[i],&cur_struct_unipath_tmp);
		ukN  += cur_struct_unipath_tmp.len * (p_unipathcnt[i] + 1);

	}
	cout << "ukN : " << ukN << endl;
	free(p_unipath);
	free(p_unipathcnt);

}

void test_lendis(struct dBG * p_dBG, char * p_dbg_path, int unidivnum)	//test unipath length distribution  unidivnum-->divide minlen~maxlen to unidivnum block
{
	struct bit256KmerPara bit_para;
	get_para(&bit_para,p_dBG->L);

	//1)解析文件名：从p_dbg_path文件中读入4个文件名，分别为：
	//kmer026, kmer026ad, unipath026, nomLeu2.fa
	char ** p_dbg_file;
	uint32_t fN;
	fN=get_Dbg_file_name(p_dbg_path,&p_dbg_file);

	cout << fN << endl;
	cout << p_dbg_file[0] << endl;
	cout << p_dbg_file[1] << endl;
	cout << p_dbg_file[2] << endl;
	cout << p_dbg_file[3] << endl;

	//3)装载分叉kmer
	uint64_t kbN;
	FILE* fp_branched_kmer_file;
	fp_branched_kmer_file=fopen(p_dbg_file[0],"rb");
	fseek(fp_branched_kmer_file,0,2);
	kbN = ftell(fp_branched_kmer_file)/(sizeof(uint64_t)*bit_para.kmer64Len);
	cout << "kbN : " << kbN << endl;
	fclose(fp_branched_kmer_file);

	//4）计算 branched kmer ad 中的边的个数
	uint64_t addegree = 0;
	FILE* fp_branched_kmer_file_ad;
	fp_branched_kmer_file_ad = fopen(p_dbg_file[1],"rb");
	uint8_t adinf;
	for(uint32_t i = 0; i < kbN; i++)
	{
		fread(&adinf,sizeof(uint8_t),1,fp_branched_kmer_file_ad);
		calc_onesinchar(adinf,addegree);
	}
	fclose(fp_branched_kmer_file_ad);


	//5)计算unipath数量
	uint32_t uN;
	uN = get_total_unipath_number(p_dbg_file);
	cout <<"uN : "<< uN << endl;

	//6)计算边的数量
	uint64_t edgenum;
	edgenum = uN + addegree / 2;
	cout << "Total edge num is: " << edgenum << endl;

	//7）读unipath
	FILE * fp_unipath;
	fp_unipath=fopen(p_dbg_file[2],"rb");
	uint64_t *p_unipath;
	p_unipath=(uint64_t*)malloc(sizeof(uint64_t)*uN);
	fread(p_unipath,sizeof(uint64_t),uN,fp_unipath);
	fclose(fp_unipath);

	FILE * fp_unipathcnt;
	fp_unipathcnt=fopen(p_dbg_file[4],"rb");
	uint64_t *p_unipathcnt;
	p_unipathcnt=(uint64_t*)malloc(sizeof(uint64_t)*uN);
	fread(p_unipathcnt,sizeof(uint64_t),uN,fp_unipathcnt);
	fclose(fp_unipath);

	//8)计算super unipath总长度
	struct unipath cur_struct_unipath_tmp;
	uint64_t superunipath_len = 0l;
	uint64_t ukN = 0;
	uint32_t unilenmin,unilenmax;
	uint32_t unicurlen;
	get_unipath_struct(p_unipath[0],&cur_struct_unipath_tmp);
	unicurlen = cur_struct_unipath_tmp.len+p_dBG->L - 1;
	unilenmin = unicurlen;
	unilenmax = unicurlen;
	for(uint32_t i=0;i<uN;i++)
	{
		get_unipath_struct(p_unipath[i],&cur_struct_unipath_tmp);
		ukN += cur_struct_unipath_tmp.len;
		superunipath_len += cur_struct_unipath_tmp.len+p_dBG->L;
		unicurlen = cur_struct_unipath_tmp.len+p_dBG->L - 1;
		if(unicurlen < unilenmin)
		{
			unilenmin = unicurlen;
		}
		if(unicurlen > unilenmax)
		{
			unilenmax = unicurlen;
		}
	}
	cout << "ukN : " << ukN << endl;
	cout << "SUL : " << superunipath_len << endl;

	//9) 计算unipath的长度分布区间  其中unidivnum是函数的输入参数  将unipath按长度分配至不同区间
	uint32_t interval = (uint32_t)ceil((double)(unilenmax - unilenmin) / unidivnum);
	uint32_t lenarray[unidivnum] = {0};
	uint32_t pos;
	for(uint32_t i=0;i<uN;i++)
	{
		get_unipath_struct(p_unipath[i],&cur_struct_unipath_tmp);
		unicurlen = cur_struct_unipath_tmp.len+p_dBG->L - 1;
		pos = (unicurlen - unilenmin) / interval;
		lenarray[pos] += p_unipathcnt[i] + 1;

	}
	printf("shortest unipath length : %d\n",unilenmin);
	printf("longest unipath length : %d\n",unilenmax);
	for(uint32_t i=0;i<unidivnum;i++)
	{
		printf("%d  ",lenarray[i]);
	}
	printf("\n");
	free(p_unipath);
	free(p_unipathcnt);

}

void Gen_unipathLenInf(struct dBG * p_dBG, char * p_dbg_path)	//利用CDBGO生成的counter文件  生成unipath长度和频率分布文件
{
	struct bit256KmerPara bit_para;
	get_para(&bit_para,p_dBG->L);
	char ** p_dbg_file;
	uint32_t fN;
	fN = get_Dbg_file_name(p_dbg_path,&p_dbg_file);

	uint32_t uN;
	uN = get_total_unipath_number(p_dbg_file);
	cout <<"uN : "<< uN << endl;

	uint32_t unicurlen;
	struct unipath cur_struct_unipath_tmp;
	FILE * fp_unipath;
	fp_unipath = fopen(p_dbg_file[2],"rb");
	uint64_t *p_unipath;
	p_unipath = (uint64_t*)malloc(sizeof(uint64_t)*uN);
	fread(p_unipath,sizeof(uint64_t),uN,fp_unipath);

	FILE * fp_unipathcnt;
	fp_unipathcnt=fopen(p_dbg_file[4],"rb");
	uint64_t *p_unipathcnt;
	uint64_t unipathcnt_tmp;
	p_unipathcnt=(uint64_t*)malloc(sizeof(uint64_t)*uN);
	fread(p_unipathcnt,sizeof(uint64_t),uN,fp_unipathcnt);
	fclose(fp_unipath);

	char fileoutname[32] = {0};
	sprintf(fileoutname,"%s%s",p_dbg_file[2],"_leninf");
	FILE *fp_unipathout;
	fp_unipathout = fopen(fileoutname,"wb");
	for(uint32_t i = 0; i < uN; i++)
	{
		unipathcnt_tmp = p_unipathcnt[i] + 1;
		fwrite(&unipathcnt_tmp,sizeof(uint64_t),1,fp_unipathout);
		get_unipath_struct(p_unipath[i],&cur_struct_unipath_tmp);
		unicurlen = cur_struct_unipath_tmp.len+p_dBG->L - 1;
		fwrite(&unicurlen,sizeof(uint32_t),1,fp_unipathout);
		fwrite("\n",sizeof(char),1,fp_unipathout);
	}
}

void merge_superunipath(struct dBG * p_dBG, char * p_dbg_path,uint8_t ** p_uni_adinf,uint32_t **p_pre_arrlen)	//连接unipath成为superunipath
{
	struct bit256KmerPara bit_para;
	get_para(&bit_para,p_dBG->L);

	//1)解析文件名：从p_dbg_path文件中读入4个文件名，分别为：
	//kmer026, kmer026ad, unipath026, nomLeu2.fa
	char ** p_dbg_file;
	uint32_t fN;
	fN = get_Dbg_file_name(p_dbg_path,&p_dbg_file);

	cout << fN << endl;
	cout << p_dbg_file[0] << endl;
	cout << p_dbg_file[1] << endl;
	cout << p_dbg_file[2] << endl;
	cout << p_dbg_file[3] << endl;

	//2)计算unipath的数量
	uint32_t uN;
	uN=get_total_unipath_number(p_dbg_file);
	cout <<"the total unipath number is : "<< uN << endl;

	//3）读unipath
	FILE * fp_unipath;
	fp_unipath=fopen(p_dbg_file[2],"rb");
	uint64_t *p_unipath;
	p_unipath=(uint64_t*)malloc(sizeof(uint64_t)*uN);
	fread(p_unipath,sizeof(uint64_t),uN,fp_unipath);
	fclose(fp_unipath);

	//4)计算super unipath总长度
	char * p_ref;
	uint32_t ref_len;
	uint32_t cur_ref_id;
	struct unipath cur_struct_unipath_tmp;
	cur_ref_id=0;
	char *p_superunipath = NULL;
	uint64_t superunipath_len = 0l;
	for(uint32_t i=0;i<uN;i++)
	{
		get_unipath_struct(p_unipath[i],&cur_struct_unipath_tmp);
		superunipath_len += cur_struct_unipath_tmp.len+p_dBG->L;
	}
	p_superunipath = (char *)malloc(sizeof(char) * superunipath_len);
	cout << "the calculated length of super unipath is : " << superunipath_len << endl;

	//5）生成super unipath
	char * p_char_tmp;
	cur_ref_id = 0;
	ReadSeq(&p_ref,&ref_len,p_dbg_file[cur_ref_id+3]);
	uint32_t cur_unipth_pos_on_super=0;

	uint8_t *uni_adinf;
	uint32_t *pre_arrlen;
	uni_adinf = (uint8_t *)malloc(sizeof(uint8_t)*uN);
	pre_arrlen = (uint32_t *)malloc(sizeof(uint32_t)*uN);

	for(uint32_t i=0;i<uN;i++)
	{
		uint32_t cur_unipth_len;
		get_unipath_struct(p_unipath[i],&cur_struct_unipath_tmp);
		if(cur_struct_unipath_tmp.ref_id!=cur_ref_id)
		{
			free(p_ref);
			cur_ref_id = cur_struct_unipath_tmp.ref_id;
			ReadSeq(&p_ref,&ref_len,p_dbg_file[cur_ref_id+3]);
		}
		p_char_tmp = p_ref + cur_struct_unipath_tmp.start;

		cur_unipth_len=cur_struct_unipath_tmp.len + p_dBG->L - 1;
		get_unipath_kmer_ad(p_char_tmp,cur_unipth_len,&uni_adinf[i]);
		pre_arrlen[i] = (i == 0) ? cur_unipth_len : (pre_arrlen[i-1] + cur_unipth_len);
		strncpy(p_superunipath+cur_unipth_pos_on_super,p_char_tmp,cur_unipth_len);
		cur_unipth_pos_on_super+=cur_unipth_len;
		p_superunipath[cur_unipth_pos_on_super]='U';
		cur_unipth_pos_on_super++;
	}
	*p_uni_adinf = uni_adinf;
	*p_pre_arrlen = pre_arrlen;
	free(p_ref);
	free(p_unipath);
	cout << "the length of generated super unipath is :"  << cur_unipth_pos_on_super << endl;

	//6）输出unipath
	FILE * out_superunipath;
	out_superunipath=fopen("superUnipath","wb");
	fwrite(p_superunipath,sizeof(char),superunipath_len,out_superunipath);
	fclose(out_superunipath);
	free(p_superunipath);
	cout << "merge superunipath over !" << endl;
}

void Test_dBG_Attribute(struct dBG * p_dBG, char * p_dbg_path)
{
	struct bit256KmerPara bit_para;
	get_para(&bit_para,p_dBG->L);

	char ** p_dbg_file;
	uint32_t fN;
	fN=get_Dbg_file_name(p_dbg_path,&p_dbg_file);
	//1)总的kmer个数
	uint64_t bkN, ukN,tN;
	tN = get_total_kmer_number(p_dbg_file, &bkN, &ukN);
    cout << "the branched kmer number is:" << bkN <<endl;
    cout << "the un-branched kmer number is:" << ukN << endl;
    cout <<"the total kmer number is:" << tN << endl;

	//2)unipath个数
	uint32_t uN;
	uN = get_total_unipath_number(p_dbg_file);
	cout <<"the total unipath number is : "<< uN << endl;

	//3）计算 branched kmer ad 中的边的个数
	uint64_t addegree = 0;
	FILE* fp_branched_kmer_file_ad;
	fp_branched_kmer_file_ad = fopen(p_dbg_file[1],"rb");
	uint8_t adinf;
	for(uint32_t i = 0; i < bkN; i++)
	{
		fread(&adinf,sizeof(uint8_t),1,fp_branched_kmer_file_ad);
		calc_onesinchar(adinf,addegree);
	}
	fclose(fp_branched_kmer_file_ad);

	//4) 边的数量
	uint64_t edgenum;
	edgenum = uN + addegree / 2;
	cout << "Total edge num is: " << edgenum << endl;

	//5）读unipath
	FILE * fp_unipath;
	fp_unipath=fopen(p_dbg_file[2],"rb");
	uint64_t *p_unipath;
	p_unipath=(uint64_t*)malloc(sizeof(uint64_t)*uN);
	fread(p_unipath,sizeof(uint64_t),uN,fp_unipath);
	fclose(fp_unipath);

	char * p_ref;
	uint32_t ref_len;
	uint32_t cur_ref_id;
	struct unipath cur_struct_unipath_tmp;
	cur_ref_id=0;
	uint64_t superunipath_len = 0l;
	for(uint32_t i=0;i<uN;i++)
	{
		get_unipath_struct(p_unipath[i],&cur_struct_unipath_tmp);
		superunipath_len += cur_struct_unipath_tmp.len+p_dBG->L-1;
	}
	cout << "the calculated length of super unipath is : " << superunipath_len << endl;
}
