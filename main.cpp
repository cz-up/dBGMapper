/*
 * main.c
 *
 *  Created on: Sep 12, 2019
 *      Author: pluto
 */

#include "basic.h"
#include "read.h"
#include "load_DBG.h"
#include "load_DBG_full.h"
#include "Test_dBG_para.h"
#include "seeding.h"
#include "BPlusTree_full.h"

int main(int argc, char** argv)
{

	char *p_dbg_path;
	char *kmertest,*alignseq;
	struct dBG p_test;
	uint32_t extlen;
	char dir;
	for(uint32_t i=1;i<argc;i=i+2)
	{
		if(argv[i][0]=='-'&&argv[i][1]=='r')//reference
		{
			p_dbg_path=argv[i+1];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='L')//reference
		{
			p_test.L=atoi(argv[i+1]);
		}
		if(argv[i][0]=='-'&&argv[i][1]=='e')//reference
		{
			extlen=atoi(argv[i+1]);
		}
		if(argv[i][0]=='-'&&argv[i][1]=='d')//reference
		{
			dir = *argv[i+1];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='s')//reference
		{
			kmertest = argv[i+1];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='a')//reference
		{
			alignseq = argv[i+1];
		}

	}
//	struct timeval tvs,tve;
//	gettimeofday(&tvs,NULL);
	char nindex[] = "./nindex";
	char rindex[] = "./rindex";
	struct sFMindex nFMidx,rFMidx,FMtmp;
	struct build_para bd_para;
	read_bfile2mem(nindex,&nFMidx,0);
	read_bfile2mem(rindex,&rFMidx,0);

	cout <<"start..."<<endl;

	struct para_dBGindex sdBGindex;
	gen_dBG_index(&p_test, &sdBGindex, p_dbg_path,1);

	struct bit256KmerPara bit_para;
	get_para(&bit_para,p_test.L);
	uint32_t unitperkmer;
	unitperkmer = bit_para.kmer64Len;
//	char *kmer = "CCATGGCTGCTTTTCG";
//	char *kmer = "CAGGCAGGGGCAGGTG";
	struct TPTnode rootnode;
	rootnode.c = 'R';
	rootnode.offset = 0;
	rootnode.edarry = (uint32_t*)malloc(sizeof(uint32_t)*(2*1+1));
	rootnode.edarry[0] = 0;
	for(int i = 1; i < 2*1+1; i++)
	{
		rootnode.edarry[i] = i-1;
	}
	for(uint32_t ii = 0 ; ii < 4; ii++)
	{
		rootnode.p_child[ii] = NULL;
	}
	if(dir == 'O')
	{
		FMtmp = rFMidx;
	}
	if(dir == 'I')
	{
		FMtmp = nFMidx;
	}
	ext_treenode(bit_para, &rootnode, sdBGindex, FMtmp, dir,
			extlen, kmertest, alignseq, 1);
	cout << "ext_treenode finished!" << endl;
	char *extseq = new char[6]();
	struct seedext sedextest;
	init_seedext(rootnode, &sedextest, kmertest);
	cout << "init_seedext" << endl;
	print_extree(rootnode, extseq, NULL);
	cout << endl;
	cout << "print_extree done" << endl;
//	calc_seedextpara(&sedextest, "CCAC" , 1, nFMidx, rFMidx, bd_para);
//	for(uint32_t i = 0; i < strlen(sedextest.seqext[0]); i++)
//	{
//		for(uint32_t j = 0; j < 3; j++)
//		{
//			cout << sedextest.p3_extedarr[0][i][j] << " ";
//		}
//		cout << endl;
//	}
//	cout << "print the extend char frequency." << endl;
//	for(uint32_t i = 0; i < sedextest.num; i++)
//	{
//		for(uint32_t j = 0; j < strlen(sedextest.seqext[i]); j++)
//		{
//			cout << sedextest.p2_extchcnt[i][j] << " ";
//		}
//		cout << endl;
//	}
//	cout << "calc_seedext" << endl;
//	free_seedext(&sedextest);
//	cout << "free_seedext" << endl;
//	cout << "print_extree finished" << endl;
//	destory_extree(&rootnode);
//	cout << "destory_extree finished!" << endl;
	cout << "gen dBG over" << endl;
	free_dBGindex(&sdBGindex);
//	Test_dBG_Attribute(&p_test, p_dbg_path);

	char *funname = "gen_dBG_index";

	cout << "end..."<< endl;
//	gettimeofday(&tve,NULL);
//	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
//	cout << funname  << " time is: "<< span << endl;

	printf("%s is over!\n",funname);
//	getchar();
	return 0;
}
