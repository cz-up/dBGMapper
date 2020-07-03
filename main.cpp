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
#include "method.h"

int main(int argc, char** argv)
{

	char *p_dbg_path;
	char *kmertest,*alignseq;
	struct dBG p_test;
	uint32_t extlen;
	char dir;
	uint32_t tau;
	uint32_t loop;
	uint32_t poset = 0;
	uint32_t read_s;
	for(uint32_t i=1;i<argc;i=i+2)
	{
		if(argv[i][0]=='-'&&argv[i][1]=='r')//
		{
			p_dbg_path=argv[i+1];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='L')//
		{
			p_test.L=atoi(argv[i+1]);
		}
		if(argv[i][0]=='-'&&argv[i][1]=='e')//
		{
			extlen=atoi(argv[i+1]);
		}
		if(argv[i][0]=='-'&&argv[i][1]=='d')//
		{
			dir = *argv[i+1];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='s')//
		{
			kmertest = argv[i+1];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='a')//
		{
			alignseq = argv[i+1];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='t')//
		{
			tau = atoi(argv[i+1]);
		}
		if(argv[i][0]=='-'&&argv[i][1]=='l')//
		{
			loop = atoi(argv[i+1]);
		}
		if(argv[i][0]=='-'&&argv[i][1]=='c')//
		{
			read_s = atoi(argv[i+1]);
		}

	}
//	struct PH_Node *PH_NodeA;
//	generate_PHNArray(&PH_NodeA, tau);
//	struct PH_Node NodeT;
//	NodeT.start_seed_id = p_test.L;
//	NodeT.end_seed_id = loop;
//	NodeT.tau = read_s;
//	int phnindex = PHNode_index(PH_NodeA, NodeT);
//	cout << "phnode index:" << phnindex << endl;
//	getchar();

	struct timeval tvs,tve;
//	gettimeofday(&tvs,NULL);
	char nindex[] = "./nindex";
	char rindex[] = "./rindex";
	struct sFMindex nFMidx,rFMidx;
	struct build_para bd_para;
	read_bfile2index(nindex,&nFMidx,0);
	read_bfile2index(rindex,&rFMidx,0);
//
//	cout <<"start..."<<endl;
//
	struct bit256KmerPara bit_para;
	get_para(&bit_para,p_test.L);
	struct para_dBGindex sdBGindex;
	char *ref;
	gen_dBG_index(bit_para, &sdBGindex, p_dbg_path,1,&ref);
//
	uint32_t unitperkmer;
	unitperkmer = bit_para.kmer64Len;
	struct TPTnode rootnode;
	struct seed_extpara ext_set;
	bool unipathflag = false;
	ext_set.dir = dir;
	ext_set.orignseq = kmertest;
	ext_set.alignseq = alignseq;
	ext_set.tau = tau;
	ext_set.bit_para = bit_para;
	ext_set.sdBGidx = sdBGindex;
	ext_set.ponunipath = false;
	if(ext_set.dir == 'I')
	{
		ext_set.FMidx = nFMidx;
	}
	if(ext_set.dir == 'O')
	{
		ext_set.FMidx = rFMidx;
	}
	gettimeofday(&tvs,NULL);
	init_rootnode(&rootnode, ext_set, kmertest);
	for(uint32_t i = 0; i < loop; ++i)
	{
		ext_treenode(ext_set, &rootnode, kmertest,extlen+tau);
	}
	cout << "ext_treenode finished!" << endl;
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout << "time is: "<< span << endl;
	char *extseq = new char[extlen+2]();
	print_extree(rootnode, extseq);
	cout << "print_extree done!\n";
	/*
//
//	struct seed_segment p_seedsegtmp;
//	char *seedseq = new char[extlen+2]();
//	print_specificlen(rootnode, ext_set, extlen, seedseq,&p_seedsegtmp);
//	cout << "print_specificlen done!\n";
//	cout << "seed segment seq:" << p_seedsegtmp.seedseq << endl;
//	destory_extree(&rootnode);
//	cout << "gen dBG over" << endl;

//	char *readtmp = "CAACTTGACAGGGGCATGCA";
//	char *reftmp = "CAACTTGACTGGGGCATGCTGG";
//	unsigned int _tautmp = 2;
//	int ret = min_lastlineEd(readtmp, reftmp, strlen(readtmp), strlen(readtmp)+_tautmp, _tautmp);
//	cout << ret << endl;
//	getchar();
//	free_dBGindex(&sdBGindex);

	char *funname = "generate_seeds";
	cout <<"start..."<<endl;
	int ref_len = strlen(ref) - read_s - 1;
	double timetmp = 0;
	//test mapperStruct
	int _cntcand = 0;
	gettimeofday(&tvs,NULL);
	struct timeval tvs1,tve1;
	for(int i = 0; i < loop; ++i)
	{
		char *str2save = new char[read_s+1];
		srand((unsigned)time( NULL));
		int pos = rand() % ref_len;
		pos = poset ? poset : pos;
		strncpy(str2save,ref+pos,read_s);
//		for(int i = 0; i < 2; ++i)
//		{
//			int pos = rand() % strlen(str2save);
////			int pos = rand() % (i+1) * 10;
//			switch(str2save[pos])
//			{
//			case 'A':
//				str2save[pos] = 'C';
//				break;
//			case 'C':
//				str2save[pos] = 'G';
//				break;
//			case 'G':
//				str2save[pos] = 'T';
//				break;
//			default:
//				str2save[pos] = 'A';
//				break;
//			}
//		}
		vector<ms_seed>  vseed;
		vector<ms_candidate>  vcand;
		vector<ms_result>  vrslt;
		generate_seeds(str2save,tau,nFMidx,vseed);
//		printvec_seed(vseed);
		generate_candidate(vseed, nFMidx, vcand);
		_cntcand += vcand.size();
		gettimeofday(&tvs1,NULL);
		verification(vseed, vcand, vrslt, str2save, ref, tau);
		gettimeofday(&tve1,NULL);
		timetmp += tve1.tv_sec-tvs1.tv_sec + (tve1.tv_usec-tvs1.tv_usec) / 1000000.0;
//		printvec_result(vrslt);
		free_seed(vseed);
		delete []str2save;
//		vector<ms_seed>().swap(vseed);

	}
	cout << "average size of candidate:" << _cntcand / loop << endl;
//	printvec_candidate(vcand);
//	printvec_result(vrslt);
	//
	cout << "end..."<< endl;
	gettimeofday(&tve,NULL);
	free_FMindex(&nFMidx);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout << "time is: "<< span << endl;
	cout << timetmp << endl;
	cout << "verification cost " << timetmp / span * 100 << "%" << endl;

	printf("%s is over!\n",funname);
	*/
	return 0;
}
