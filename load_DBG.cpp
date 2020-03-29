/*
 * load_DBG.cpp
 *
 *  Created on: Sep 16, 2019
 *      Author: bio
 */

//装载dBG图
//输入数据为：kmer026, kmer026ad, unipath026, nomLeu2.fa
//装载的数据结构为：
//1）hash表，实现kmer与id的对应关系
//2）kmer_detail数组，存储kmer的具体信息，数组的下标与kmer的id对应
//3）unipath数组，存储unipath字符串，数组下标对应unipath的id

#include "load_DBG.h"
#include "Hash.h"

uint32_t get_Dbg_file_name(char* p_dbg_path, char *** p_dbg_file)
{
	uint32_t fN = 0;
	FILE * fp = NULL; //文件指针。
    int32_t c, lc=0; //c为文件当前字符，lc为上一个字符，供结尾判断用。
    uint32_t line_num = 0; //行数统计
    fp = fopen(p_dbg_path, "r");//以只读方式打开文件。
    if(fp == NULL)
    {
    	printf("%s is not exist!\n",p_dbg_path);
    	return 0;
    }
    while((c = fgetc(fp)) != EOF) //逐个读入字符直到文件结尾
    {
        if(c == '\n') //统计行数。
        {
        	line_num++;
        }
        lc = c; //保存上一字符。
    }
    if(lc != '\n')//处理末行
    {
    	line_num++;
    }
    rewind(fp); //重新指向文件开始
    char ** p_tmp = (char **)malloc(sizeof(char *) * line_num);
    char path_tmp[32] = {0};
    for( ; fN < line_num; fN++)
    {
    	char *name_tmp = (char *)malloc(sizeof(char) * 32);
		fgets(path_tmp,32,fp);
		strcat(name_tmp,path_tmp);
		name_tmp[strlen(name_tmp)-1] = '\0'; // replace '\n' with '\0'
		*(p_tmp+fN) = name_tmp;
    }
    fclose(fp);
    *p_dbg_file = p_tmp;
	return fN;
}

uint64_t get_total_kmer_number(char ** p_dbg_file, uint64_t *kN, uint64_t *ukN)
{
	uint64_t tN = 0;
	uint64_t kN_b = 0;
	FILE * fp0 = NULL;
	FILE * fp2 = NULL;
    uint64_t len = 0; //c为文件当前字符，lc为上一个字符，供结尾判断用。
	fp0 = fopen(p_dbg_file[0], "rb");
	fp2 = fopen(p_dbg_file[2], "rb");
	if(fp0 == NULL || fp2 == NULL)
	{
		printf("%s or %s not exist!\n",p_dbg_file[0],p_dbg_file[2]);
		return 0;
	}
	fseek(fp0,0,SEEK_END);//将文件内部的指针指向文件末尾
	len = ftell(fp0);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
	kN_b = len / sizeof(uint64_t);
    fclose(fp0);

    fseek(fp2,0,SEEK_END);//将文件内部的指针指向文件末尾
    len = ftell(fp2);
    fseek(fp2,0,0);
//    cout << len << endl;
    uint64_t *ptmp = NULL;
    ptmp = (uint64_t *)malloc(len);
    memset(ptmp,0,len);
    fread(ptmp,1,len,fp2);
    len = len / sizeof(uint64_t);
    uint32_t kN_unipath;
    kN_unipath=0;
    for(uint32_t i = 0; i < len; i++)
    {
    	uint64_t tmp;
    	tmp=(*(ptmp+i) >> 8 ) & 0x1ffff;
    	kN_unipath = kN_unipath + tmp;
    }

    tN = kN_b+kN_unipath;
    *kN = kN_b;
    *ukN = kN_unipath;

    free(ptmp);
    fclose(fp2);
	return tN;
}

uint32_t get_total_unipath_number(char ** p_dbg_file)
{
	uint32_t uN = 0;
	FILE * fp = NULL;
    uint32_t len = 0;
	fp = fopen(p_dbg_file[2],"r");
	if(fp == NULL)
	{
		printf("%s is not exist!\n",p_dbg_file[2]);
		return 0;
	}
	fseek(fp,0,SEEK_END);//将文件内部的指针指向文件末尾
	len = ftell(fp);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
	uN = len / sizeof(uint64_t);
    fclose(fp);
	return uN;
}

void get_unipath_struct(uint64_t a,struct unipath *b)
{
	b->start = a >> 25;
	b->len = (a >> 8) & 0x1ffff;
	b->ref_id = a & 0xff;
}

void get_unipath_kmer_ad(char* p_cur,uint32_t kmer_len,uint8_t * ad)
{
	*ad = 0;
	switch(*(p_cur-1))
	{
		case 'A':
			*ad += 0x80;
			break;
		case 'C':
			*ad += 0x40;
			break;
		case 'G':
			*ad += 0x20;
			break;
		case 'T':
			*ad += 0x10;
			break;
		default:
			break;
	}
	switch(*(p_cur+kmer_len))
	{
		case 'A':
			*ad += 0x08;
			break;
		case 'C':
			*ad += 0x04;
			break;
		case 'G':
			*ad += 0x02;
			break;
		case 'T':
			*ad += 0x01;
			break;
		default:
			break;
	}
}

void save_kmer_details(struct kmer_detail * a,\
		uint8_t is_branched,\
		uint8_t ad,\
		uint32_t unipath_id,\
		uint16_t unipath_offset)
{
	a->is_branched = is_branched;
	a->ad = ad;
	a->unipath_id = unipath_id;
	a->unipath_offset = unipath_offset;
}

/*
void loadDbg(struct dBG * p_dBG, char * p_dbg_path)
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

	struct NodeBit** p_kmer_root_tmp;
	p_kmer_root_tmp=bit256initialHashFTable();
	struct kmer_detail *p_kmer_detail_tmp;
	char ** p_unipath_tmp;

	//2)计算kmer的总数量
	uint64_t tN;
	tN=get_total_kmer_number(p_dbg_file);
	p_kmer_detail_tmp=(struct kmer_detail *)malloc(sizeof(struct kmer_detail )*tN);

	cout << tN << endl;

	//3)计算unipath的数量
	uint32_t uN;
	uN=get_total_unipath_number(p_dbg_file);
	p_unipath_tmp=(char**)malloc(sizeof(char*)*uN);

	cout << uN << endl;

	//4)装载分叉kmer
	uint64_t kbN;
	FILE* fp_branched_kmer_file;
	fp_branched_kmer_file=fopen(p_dbg_file[0],"rb");
	fseek(fp_branched_kmer_file,0,2);
	kbN=ftell(fp_branched_kmer_file)/(sizeof(uint64_t)*bit_para.kmer64Len);
	fseek(fp_branched_kmer_file,0,0);

	FILE* fp_branched_kmer_file_ad;
	fp_branched_kmer_file_ad=fopen(p_dbg_file[1],"rb");

	struct nodeBit node_tmp;
	node_tmp.hashValue=(uint64_t*)malloc(sizeof(uint64_t)*4);
	for(uint64_t i=0;i<kbN;i++)
	{
		uint8_t  kmer_ad_tmp;
		fread(node_tmp.hashValue,sizeof(uint64_t),bit_para.kmer64Len,fp_branched_kmer_file);
		fread(&kmer_ad_tmp,sizeof(uint8_t),1,fp_branched_kmer_file_ad);
		node_tmp.arrayID=i;
		bit256insertHashFTable(p_kmer_root_tmp,node_tmp,bit_para);
		//add the kmer information here
		save_kmer_details(p_kmer_detail_tmp+i,1,kmer_ad_tmp,0,0);
	}
	fclose(fp_branched_kmer_file);
	fclose(fp_branched_kmer_file_ad);

	/*-------------------  2019.10.9   merge unipath to super unipath sting  -------------------*/
//	uint64_t super_unipath_len = 0;
//	for(uint32_t i=0;i<uN;i++)
//	{
//		fread(&cur_unipath,sizeof(uint64_t),1,fp_unipath);
//		get_unipath_struct(cur_unipath,&cur_struct_unipath_tmp);
//		if(cur_struct_unipath_tmp.ref_id!=cur_ref_id)
//		{
//			free(p_ref);
//			cur_ref_id=cur_struct_unipath_tmp.ref_id;
//			ReadSeq(&p_ref,&ref_len,p_dbg_file[cur_ref_id+3]);
//		}
//		//存储当前unipath到unipath数组上，unipath长度为含有的不分叉的kmer的数量
//		char * p_char_tmp;
//		p_char_tmp=p_ref+cur_struct_unipath_tmp.start;
//		p_unipath_tmp[i]=(char*)malloc(sizeof(char)*(cur_struct_unipath_tmp.len+p_dBG->L-1));
//		for(uint32_t j=0;j<cur_struct_unipath_tmp.len+p_dBG->L-1;j++)
//		{
//			p_unipath_tmp[i][j]=p_char_tmp[j];
//			super_unipath_len++;
//		}
//	}
//	char * p_super_unipath = NULL;
//	p_super_unipath = (char *)malloc(sizeof(char) * super_unipath_len);
//	strcpy(p_super_unipath,p_unipath_tmp[0]);
//	strcat(p_super_unipath,"$");
//	for(uint32_t i=1;i<uN;i++)
//	{
//		strcat(p_super_unipath,p_unipath_tmp[i]);
//		strcat(p_super_unipath,"$");
//	}
	/*-------------------  2019.10.9   merge unipath to super unipath sting  -------------------*/


//	//5)装载unipath上的kmer
//	char * p_ref;
//	uint64_t cur_kmer_id;
//	cur_kmer_id=kbN;
//	uint64_t cur_unipath;
//	uint32_t ref_len;
//	uint32_t cur_ref_id;
//	uint32_t total_ref=fN-3;
//	struct unipath cur_struct_unipath_tmp;
//	//其中包含kmer026, kmer026ad, unipath026三个文件
//	//5.1)读reference序列
//	//5,2)一条一条处理unipath，当unipath在下一条reference上时，
//	//free当前的p_ref,病装载新的ref；
//	FILE * fp_unipath;
//	fp_unipath=fopen(p_dbg_file[2],"rb");
//	cur_ref_id=0;
//	ReadSeq(&p_ref,&ref_len,p_dbg_file[cur_ref_id+3]);
//	for(uint32_t i=0;i<uN;i++)
//	{
//		fread(&cur_unipath,sizeof(uint64_t),1,fp_unipath);
//		get_unipath_struct(cur_unipath,&cur_struct_unipath_tmp);
//		if(cur_struct_unipath_tmp.ref_id!=cur_ref_id)
//		{
//			free(p_ref);
//			cur_ref_id=cur_struct_unipath_tmp.ref_id;
//			ReadSeq(&p_ref,&ref_len,p_dbg_file[cur_ref_id+3]);
//		}
//		//存储当前unipath到unipath数组上，unipath长度为含有的不分叉的kmer的数量
//		char * p_char_tmp;
//		p_char_tmp=p_ref+cur_struct_unipath_tmp.start;
//		p_unipath_tmp[i]=(char*)malloc(sizeof(char)*(cur_struct_unipath_tmp.len+p_dBG->L-1));
//		for(uint32_t j=0;j<cur_struct_unipath_tmp.len+p_dBG->L-1;j++)
//		{
//			p_unipath_tmp[i][j]=p_char_tmp[j];
//		}
//		//解析当前unipath上的kmer并存储
//		uint8_t  kmer_ad_tmp;
//		cal_hash_value_directly_256bit(p_char_tmp,node_tmp.hashValue,bit_para);
//		node_tmp.arrayID=cur_kmer_id;
//		bit256insertHashFTable(p_kmer_root_tmp,node_tmp,bit_para);
//		get_unipath_kmer_ad(p_char_tmp,p_dBG->L,&kmer_ad_tmp);
//		save_kmer_details(p_kmer_detail_tmp+cur_kmer_id,0,kmer_ad_tmp,i,0);
//		cur_kmer_id++;
//		for(uint32_t j=1;j<cur_struct_unipath_tmp.len;j++)
//		{
//			cal_hash_value_indirectly_256bit(p_char_tmp+j,node_tmp.hashValue,node_tmp.hashValue,bit_para);
//			node_tmp.arrayID=cur_kmer_id;
//			bit256insertHashFTable(p_kmer_root_tmp,node_tmp,bit_para);
//			get_unipath_kmer_ad(p_char_tmp+j,p_dBG->L,&kmer_ad_tmp);
//			save_kmer_details(p_kmer_detail_tmp+cur_kmer_id,0,kmer_ad_tmp,i,j);
//			cur_kmer_id++;
//		}
//	}
//	free(p_ref);
//
//	//final
//	p_dBG->p_kmer_root=p_kmer_root_tmp;
//	p_dBG->p_kmer_detail=p_kmer_detail_tmp;
//	p_dBG->p_unipath=p_unipath_tmp;
//	p_dBG->kN=tN;
//	p_dBG->uN=uN;




