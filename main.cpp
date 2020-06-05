/*
 * main.c
 *
 *  Created on: Sep 12, 2019
 *      Author: pluto
 */

#include "analyse_dBG.h"
#include "basic.h"
#include "Binary_Search.h"
#include "read.h"
#include "load_DBG_full.h"
#include "seeding.h"
#include "mapperStruct.h"
#include "prealignment.h"

int main(int argc, char** argv)
{

	char *p_dbg_path;
	char *kmertest,*alignseq;
	struct dBG p_test;
	uint32_t extlen;
	char dir;
	uint32_t tau;
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
		if(argv[i][0]=='-'&&argv[i][1]=='t')//reference
		{
			tau = atoi(argv[i+1]);
		}

	}
	for(int i = 0; i < 1; /*strlen(kmertest)-7;*/ ++i)
	{
		char readt[16] = {0};
//		strncpy(readt,kmertest+i+2,4);
		int ret = pre_alignment(alignseq, kmertest, 2, strlen(kmertest)-7, tau);
		if(!ret)
		{
			cout << "ret != 1" << "i = " << i << endl;
		}
	}
	cout << "test done\n";
	getchar();
	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);
	char nindex[] = "./nindex";
	char rindex[] = "./rindex";
	struct sFMindex nFMidx,rFMidx;
	struct build_para bd_para;
	read_bfile2index(nindex,&nFMidx,0);
	read_bfile2index(rindex,&rFMidx,0);

	cout <<"start..."<<endl;

	struct bit256KmerPara bit_para;
	get_para(&bit_para,p_test.L);
	struct para_dBGindex sdBGindex;
	char *ref;
	gen_dBG_index(bit_para, &sdBGindex, p_dbg_path,1,&ref);

	uint32_t unitperkmer;
	unitperkmer = bit_para.kmer64Len;
//	char *kmer = "CCATGGCTGCTTTTCG";
//	char *kmer = "CAGGCAGGGGCAGGTG";
	struct TPTnode rootnode;
	struct seed_extpara ext_set;
	bool unipathflag = false;
	ext_set.dir = dir;
	ext_set.orignseq = kmertest;
	ext_set.alignseq = alignseq;
	ext_set.tau = tau;
	ext_set.bit_para = bit_para;
	ext_set.sdBGidx = sdBGindex;
	if(ext_set.dir == 'I')
	{
		ext_set.FMidx = nFMidx;
	}
	if(ext_set.dir == 'O')
	{
		ext_set.FMidx = rFMidx;
	}
//	init_rootnode(&rootnode, ext_set, kmertest);
//	ext_treenode(ext_set, &rootnode, kmertest,extlen+tau, &unipathflag);
//	cout << "ext_treenode finished!" << endl;
//	char *extseq = new char[extlen+2]();
//	print_extree(rootnode, extseq);
//	cout << "print_extree done!\n";
//
//	struct seed_segment p_seedsegtmp;
//	char *seedseq = new char[extlen+2]();
//	print_specificlen(rootnode, ext_set, extlen, seedseq,&p_seedsegtmp);
//	cout << "print_specificlen done!\n";
//	cout << "seed segment seq:" << p_seedsegtmp.seedseq << endl;
//	destory_extree(&rootnode);
//	cout << "gen dBG over" << endl;
//	free_dBGindex(&sdBGindex);
	char *funname = "generate_seeds";

	//test mapperStruct
	vector<ms_seed>  vseed;
	vector<ms_candidate>  vcand;
	vector<ms_result>  vrslt;
	generate_seeds(kmertest,ext_set.tau,nFMidx,vseed);
	generate_candidate(vseed, nFMidx, vcand);
	printvec_seed(vseed);
	printvec_candidate(vcand);
	verification(vseed, vcand, vrslt, kmertest, ref, ext_set.tau);
	printvec_result(vrslt);
	//
	cout << "end..."<< endl;
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout << "time is: "<< span << endl;

	printf("%s is over!\n",funname);
//	getchar();
	return 0;
}
