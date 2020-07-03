#include "basic.h"

void reverseq(char *seq)
{
	char *front = seq;
	char *rear = seq + strlen(seq) - 1;
	char tmp;
	while(front < rear)
	{
		tmp = *front;
		*front = *rear;
		*rear = tmp;
		front++;
		rear--;
	}
}

void kmercpy(uint64_t *dest, const uint64_t * from, uint32_t unitper)
{
	if(dest == NULL || from == NULL)
	{
		printf("kmercpy error\n");
		exit(1);
	}
	for(uint32_t i = 0; i < unitper; i++)
	{
		dest[i] = from[i];
	}
}

void get_para(struct bit256KmerPara *para1,uint64_t kmer_length)
{
	struct bit256KmerPara *p=para1;
	(*p).kmer1Len=kmer_length*2;
	(*p).remainer1to64=(*p).kmer1Len%64;
	(*p).kmer64Len=(*p).kmer1Len/64+((*p).remainer1to64?1:0);
	if((*p).remainer1to64==0)
	{
		(*p).codefor1=0;
	}
	else
	{
		(*p).codefor1=0;
		for(uint32_t i=0;i<(*p).remainer1to64-1;i++)
		{
			(*p).codefor1=(*p).codefor1|1;
			(*p).codefor1=(*p).codefor1<<1;
		}
		(*p).codefor1=(*p).codefor1|1;
	}
}


void fill_char_with_four_char(uint8_t* current,char*p)
{
	for(uint32_t i=0;i<3;i++)
	{
		switch(p[i])
		{
			case 'A':
				*current=*current<<2;
				break;
			case 'C':
				*current=*current|1;
				*current=*current<<2;
				break;
			case 'G':
				*current=*current|2;
				*current=*current<<2;
				break;
			case 'T':
				*current=*current|3;
				*current=*current<<2;
				break;
			default:
				*current=*current<<2;
				break;
		}
	}
	switch(p[3])
	{
		case 'A':
			break;
		case 'C':
			*current=*current|1;
			break;
		case 'G':
			*current=*current|2;
			break;
		case 'T':
			*current=*current|3;
			break;
		default:
			break;
	}
}
void ReadSeq_bit(uint8_t **seq1,uint64_t *seq_length,char* p_ref)
{
	uint32_t cursize;
	uint32_t maxsize;
	uint32_t addsize;

	maxsize=pow(2,20);
	addsize=pow(2,20);

	uint8_t *seq;//=seq1;
	seq=(uint8_t*) malloc (sizeof(uint8_t)*maxsize);
	cursize=maxsize;

	uint32_t buffer_size=_BufferSize_;
	char buffer_line[_BufferSize_];
	memset(buffer_line,0,buffer_size);

	FILE *fp;
	fp = fopen(p_ref,"r+");
	if(fp==NULL)
	{
		cout <<"file can not be open!" << endl;
	}

	uint64_t len=0;
	uint32_t cycle=0;
	uint32_t number_left=0;
	uint32_t number_buffer_char;
	char array_left[4];
	while (1)
	{
		if(fgets(buffer_line,buffer_size-number_left,fp)==NULL)
		{
			break;
		}
		if(buffer_line[0]=='>')
		{
			continue;
		}
		else
		{
			char buffer_line_tmp[_BufferSize_];
			if(number_left!=0)
			{
				strcpy(buffer_line_tmp,array_left);
				strcpy(buffer_line_tmp+number_left,buffer_line);
			}
			else
			{
				strcpy(buffer_line_tmp,buffer_line);
			}

			uint32_t buffer_line_char_number=strlen(buffer_line_tmp);
			if(buffer_line_tmp[buffer_line_char_number-1]=='\n')
			{
				buffer_line_char_number--;
			}
			number_left=buffer_line_char_number%4;
			number_buffer_char=buffer_line_char_number/4;
			if(number_left!=0)
			{
				for(uint32_t i=0;i<number_left;i++)
				{
					array_left[i]=buffer_line_tmp[buffer_line_char_number-(number_left-i)];
				}
				array_left[number_left]='\0';
			}

			if(len+number_buffer_char<cursize)
			{
				for(uint32_t i=0;i<number_buffer_char;i++)
				{
					uint8_t tmp=0;

					for(uint32_t j=0;j<4;j++)
					{
						if(buffer_line_tmp[4*i+j]>='a')
						{
							buffer_line_tmp[4*i+j]-=32;
						}
					}
					fill_char_with_four_char(&tmp,buffer_line_tmp+4*i);
					seq[len]=tmp;
					len++;
				}
			}
			else
			{
				seq=(uint8_t*) realloc (seq,sizeof(uint8_t)*(cursize+addsize));
				cursize=cursize+addsize;
				for(uint32_t i=0;i<buffer_line_char_number/4;i++)
				{
					uint8_t tmp=0;

					for(uint32_t j=0;j<4;j++)
					{
						if(buffer_line_tmp[4*i+j]>='a')
						{
							buffer_line_tmp[4*i+j]-=32;
						}
					}
					fill_char_with_four_char(&tmp,buffer_line_tmp+4*i);
					seq[len]=tmp;
					len++;
				}
				cout <<"add 1024*1024 byte for seq: " << cycle++ <<endl;
			}
		}
		memset(buffer_line,0,buffer_size);
	}
	*seq_length=len*4+number_left;
	if(number_left!=0)
	{
		for(uint32_t i=number_left;i<4;i++)
		array_left[i]='A';
		uint8_t tmp;
		fill_char_with_four_char(&tmp,array_left);
		seq[len]=tmp;
		len++;
	}

	seq1[0]=seq;

	uint8_t* seq_shift1=(uint8_t*) malloc (sizeof(uint8_t)*len);
	for(uint32_t i=0;i<len-1;i++)
	{
		seq_shift1[i]=seq[i]<<2;
		seq_shift1[i]=seq_shift1[i]|(seq[i+1]>>6);
	}
	seq_shift1[len-1]=seq[len-1]<<2;
	seq1[1]=seq_shift1;

	uint8_t* seq_shift2=(uint8_t*) malloc (sizeof(uint8_t)*len);
	for(uint32_t i=0;i<len-1;i++)
	{
		seq_shift2[i]=seq_shift1[i]<<2;
		seq_shift2[i]=seq_shift2[i]|(seq_shift1[i+1]>>6);
	}
	seq_shift2[len-1]=seq_shift1[len-1]<<2;
	seq1[2]=seq_shift2;

	uint8_t* seq_shift3=(uint8_t*) malloc (sizeof(uint8_t)*len);
	for(uint32_t i=0;i<len-1;i++)
	{
		seq_shift3[i]=seq_shift2[i]<<2;
		seq_shift3[i]=seq_shift3[i]|(seq_shift2[i+1]>>6);
	}
	seq_shift3[len-1]=seq_shift2[len-1]<<2;
	seq1[3]=seq_shift3;

	cout << "the length of seq is: " << *seq_length << endl;
}

void ReadSeq(char **seq1,uint32_t *seq_length,char* p_ref)
{
	uint32_t cursize;
	uint32_t maxsize;
	uint32_t addsize;

	maxsize=pow(2,20);
	addsize=pow(2,20);

	char *seq;
	seq=(char*) malloc (sizeof(char)*maxsize);
	cursize=maxsize;

	uint32_t buffer_size=256;
	char buffer_line[256] = {0};

	FILE *fp;
	fp = fopen(p_ref,"r+");
	if(fp==NULL)
	{
		cout <<"file can not be open!" << endl;
	}

	uint32_t len=0;
	while (fgets(buffer_line,buffer_size,fp)!=NULL)
	{
		if(buffer_line[0]=='>')
			continue;
		else
		{
			if(len+buffer_size<cursize)
			{
				for(uint32_t i=0;i<buffer_size;i++)
				{
					if(buffer_line[i]=='\n'||buffer_line[i]=='\0')
					{
						break;
					}
					if(buffer_line[i]>='a')
					{
						buffer_line[i]-=32;
					}
					if(buffer_line[i]!='A'&&buffer_line[i]!='C'&&buffer_line[i]!='G'&&buffer_line[i]!='T')
					{
						buffer_line[i]='A';
					}
					seq[len]=buffer_line[i];
					len++;
				}
			}
			else
			{
				seq=(char*) realloc (seq,sizeof(char)*(cursize+addsize));
				cursize=cursize+addsize;
				for(uint32_t i=0;i<buffer_size;i++)
				{
					if(buffer_line[i]=='\n'||buffer_line[i]=='\0')
					{
						break;
					}
					if(buffer_line[i]>='a')
					{
						buffer_line[i]-=32;
					}
					if(buffer_line[i]!='A'&&buffer_line[i]!='C'&&buffer_line[i]!='G'&&buffer_line[i]!='T')
					{
						buffer_line[i]='A';
					}
					seq[len]=buffer_line[i];
					len++;
				}
			}
		}
		memset(buffer_line,0,buffer_size*sizeof(char));
	}
	*seq_length=len;
	*seq1=seq;
	cout << "the length of seq is: " << len << endl;
}
uint32_t cmp256BitKmer(uint64_t*a,uint64_t*b,uint32_t len)
{
	uint32_t r=2;
	for(uint32_t i=0;i<len;i++)
	{
		if(a[i]<b[i])
		{
			r=0;
			break;
		}
		else
		{
			if(a[i]>b[i])
			{
				r=1;
				break;
			}
		}
	}
	return r;
}
void cal_hash_value_directly_256bit(char *seq,uint64_t * current,\
		struct bit256KmerPara para)
{
	uint64_t tmp;
	char *k_mer_temp=seq;
	for(uint32_t i=0;i<para.kmer64Len;i++)
	{
		char* loop_tmp=k_mer_temp+32*i;
		if(i==para.kmer64Len-1&&para.remainer1to64!=0)
		{
			tmp=0;
			for(uint32_t j=0;j<para.remainer1to64/2-1;j++)
			{
				switch(loop_tmp[j])
				{
					case 'A':
						tmp=tmp<<2;
						break;
					case 'C':
						tmp=tmp|1;
						tmp=tmp<<2;
						break;
					case 'G':
						tmp=tmp|2;
						tmp=tmp<<2;
						break;
					case 'T':
						tmp=tmp|3;
						tmp=tmp<<2;
						break;
					default:
						tmp=tmp<<2;
						break;
				}
			}
			switch(loop_tmp[para.remainer1to64/2-1])
			{
				case 'A':
					break;
				case 'C':
					tmp=tmp|1;
					break;
				case 'G':
					tmp=tmp|2;
					break;
				case 'T':
					tmp=tmp|3;
					break;
				default:
					break;
			}
			current[i]=tmp;
		}
		else
		{
			tmp=0;
			for(uint32_t j=0;j<31;j++)
			{
				switch(loop_tmp[j])
				{
					case 'A':
						tmp=tmp<<2;
						break;
					case 'C':
						tmp=tmp|1;
						tmp=tmp<<2;
						break;
					case 'G':
						tmp=tmp|2;
						tmp=tmp<<2;
						break;
					case 'T':
						tmp=tmp|3;
						tmp=tmp<<2;
						break;
					default:
						tmp=tmp<<2;
						break;
				}
			}
			switch(loop_tmp[31])
			{
				case 'A':
					break;
				case 'C':
					tmp=tmp|1;
					break;
				case 'G':
					tmp=tmp|2;
					break;
				case 'T':
					tmp=tmp|3;
					break;
				default:
					break;
			}
			current[i]=tmp;
		}
	}
}
void cal_hash_value_indirectly_256bit(char *seq,uint64_t* original,uint64_t* current,\
		struct bit256KmerPara para)
{
	char *k_mer_temp=seq;
	for(uint32_t i=0;i<para.kmer64Len-1;i++)
	{
		if(i==para.kmer64Len-2&&para.remainer1to64!=0)
		{
			current[i]=original[i]<<2;
			current[i]=current[i]|(original[i+1]>>(para.remainer1to64-2));
		}
		else
		{
			current[i]=original[i]<<2;
			current[i]=current[i]|(original[i+1]>>62);
		}
	}
	if(para.remainer1to64==0)
	{
		current[para.kmer64Len-1]=original[para.kmer64Len-1]<<2;
		switch(k_mer_temp[para.kmer1Len/2-1])
		{
			case 'A':
				//current=current|0;
				break;
			case 'C':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|1;
				//current=current<<2;
				break;
			case 'G':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|2;
				//current=current<<2;
				break;
			case 'T':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|3 ;
				//current=current<<2;
				break;
			default:
				break;
		}
	}
	else
	{
		current[para.kmer64Len-1]=original[para.kmer64Len-1]<<2;
		current[para.kmer64Len-1]=current[para.kmer64Len-1]&para.codefor1;
		switch(k_mer_temp[para.kmer1Len/2-1])
		{
			case 'A':
				//current=current|0;
				break;
			case 'C':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|1;
				//current=current<<2;
				break;
			case 'G':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|2;
				//current=current<<2;
				break;
			case 'T':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|3 ;
				//current=current<<2;
				break;
			default:
				break;
		}
	}
}
