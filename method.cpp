/*
 * method.cpp
 *
 *  Created on: Feb 16, 2019
 *      Author: bio
 */
#include "method.h"

uint32_t min_lastlineEd(char *a,char * b,int32_t x,int32_t y, uint32_t &tau)  //b is ref   1。read 长度不同时间   2.逻辑是否有问题  3.candidate数量
{
	uint32_t r;
	uint32_t ** m;
	m=(uint32_t **)malloc(sizeof(uint32_t *)*(x+1));
	for(int32_t i=0;i<x+1;i++)
	{
		m[i]=(uint32_t *)malloc(sizeof(uint32_t)*(y+1));
	}

	for(int32_t i=0;i<x+1;++i)
	{
		m[i][0]=i;
	}
	for(int32_t i=0;i<y+1;++i)
	{
		m[0][i]=i;
	}

	for(int32_t i=1;i<x+1;++i)
	{
		for(int32_t j=1;j<y+1;++j)
		{
			int32_t tmp=0;
			if(a[i-1]!=b[j-1])
			{
				tmp=1;
			}
			m[i][j]=min(m[i][j-1]+1,m[i-1][j-1]+tmp);
			m[i][j]=min(m[i][j],m[i-1][j]+1);
		}
	}
	r = m[x][0];
	uint32_t col = 0;
	for(int32_t i=1;i<x+1;++i)
	{
		if(m[x][i] <= tau)
		{
			r = m[x][i];
			col = i;
		}
	}
	for(int32_t i=0;i<x+1;++i)
	{
		free(m[i]);
		m[i] = NULL;
	}
	free(m);
	m = NULL;
	if(r <= tau)
	{
		tau -= r;
		return col;
	}
	else
	{
		return -1;
	}
}


uint32_t **originalEd(char *a,char * b,int32_t x,int32_t y)
{
	uint32_t r;
	uint32_t ** m;
	m=(uint32_t **)malloc(sizeof(uint32_t *)*(x+1));
	for(int32_t i=0;i<x+1;i++)
	{
		m[i]=(uint32_t *)malloc(sizeof(uint32_t)*(y+1));
	}

	for(int32_t i=0;i<x+1;i++)
	{
		m[i][0]=i;
	}
	for(int32_t i=0;i<y+1;i++)
	{
		m[0][i]=i;
	}

	for(int32_t i=1;i<x+1;i++)
	{
		for(int32_t j=1;j<y+1;j++)
		{
			int32_t tmp=0;
			if(a[i-1]!=b[j-1])
			{
				tmp=1;
			}
			m[i][j]=min(m[i][j-1]+1,m[i-1][j-1]+tmp);
			m[i][j]=min(m[i][j],m[i-1][j]+1);
		}
	}
	return m;
}

void free_edmatrix(uint32_t ***m, int32_t x)
{
	if(*m != NULL)
	{
		for(int32_t i=0;i<x+1;i++)
		{
			free((*m)[i]);
			(*m)[i] = NULL;
		}
		free(*m);
		*m = NULL;
	}
}

uint32_t rangeEd(char *a,char * b,int32_t x,int32_t y,int32_t ed)
{
	if((x-y > ed) || (x-y < -ed))
	{
		return (uint32_t)ed+1;
	}
	uint32_t r;
	uint32_t ** m;
	m=(uint32_t **)malloc(sizeof(uint32_t *)*(x+1));
	for(int32_t i=0;i<x+1;i++)
	{
		m[i]=(uint32_t *)malloc(sizeof(uint32_t)*(y+1));
	}

	for(int32_t i=0;i<x+1;i++)
	{
		m[i][0]=i;
	}
	for(int32_t i=0;i<y+1;i++)
	{
		m[0][i]=i;
	}

	int32_t y_start;
	int32_t y_end;
	int32_t mid_ed=ed/2;
	int32_t l_bound;
	int32_t r_bound;
	int32_t y1=1;
	if(x<=y)
	{
		l_bound=mid_ed;
		r_bound=(y-x)+mid_ed;
	}
	else
	{
		l_bound=(x-y)+mid_ed;
		r_bound=mid_ed;
	}
	for(int32_t i=1;i<x+1;i++)
	{
		y_start=max(y1,i-l_bound);
		y_end=min(i+r_bound,y);
		int32_t line_check = ed+1;
		for(int32_t j=y_start;j<=y_end;j++)
		{
			uint32_t tmp=0;
			if(j==i-l_bound && j==i+r_bound)
			{
				if(a[i-1]!=b[j-1])
				{
					tmp=1;
				}
				m[i][j]=m[i-1][j-1]+tmp;
			}
			else if(j==i-l_bound)
			{
				if(a[i-1]!=b[j-1])
				{
					tmp=1;
				}
				m[i][j]=min(m[i-1][j]+1,m[i-1][j-1]+tmp);
			}
			else if(j==i+r_bound)
			{
				if(a[i-1]!=b[j-1])
				{
					tmp=1;
				}
				m[i][j]=min(m[i][j-1]+1,m[i-1][j-1]+tmp);
			}
			else
			{
				if(a[i-1]!=b[j-1])
				{
					tmp=1;
				}
				m[i][j]=min(m[i][j-1]+1,m[i-1][j-1]+tmp);
				m[i][j]=min(m[i][j],m[i-1][j]+1);
			}
			if(m[i][j] < line_check)
			{
				line_check = m[i][j];
			}
		}
		if(line_check == (ed+1))
		{
			for(int32_t i=0;i<x+1;i++)
			{
				free(m[i]);
			}
			free(m);
			return ed+1;
		}
	}
	r = m[x][y];
	for(int32_t i=0;i<x+1;i++)
	{
		free(m[i]);
	}
	free(m);
	return r;
}

struct cell
{
	int32_t a;
	int32_t b;
};
uint32_t ** levelEd(char *a,char * b,int32_t x,int32_t y,int32_t ed)
{
	if((x-y > ed) || (x-y < -ed))
	{
		return NULL;
	}
	uint32_t r = ed+1;
	uint32_t ** m;
	m=(uint32_t **)malloc(sizeof(uint32_t *)*(x+1));
	for(int32_t i=0;i<x+1;i++)
	{
		m[i]=(uint32_t *)malloc(sizeof(uint32_t)*(y+1));
		memset(m[i],127,sizeof(uint32_t)*(y+1));
	}

//	for(uint32_t i=0;i<x+1;i++)
//	{
//		m[i][0]=i;
//	}
//	for(uint32_t i=0;i<y+1;i++)
//	{
//		m[0][i]=i;
//	}
	m[0][0] = 0;

	struct cell * p_last;
	struct cell * p_cur;
	struct cell * H_last;
	struct cell * H_cur;
	uint32_t n_last=0;
	uint32_t n_cur=0;
	H_last=(struct cell*)malloc(sizeof(struct cell)*x*y);
	H_cur=(struct cell*)malloc(sizeof(struct cell)*x*y);
	p_last = H_last;
	p_cur = H_cur;

	struct cell tmp;

	tmp.a=0;
	tmp.b=0;
	H_last[n_last]=tmp;
	n_last++;
	uint32_t len=1;
	while(1)
	{
		if(tmp.a+len>x||tmp.b+len>y)
		{
			break;
		}
		if(a[tmp.a+len-1]==b[tmp.b+len-1])
		{
			H_last[n_last].a=tmp.a+len;
			H_last[n_last].b=tmp.b+len;
			m[tmp.a+len][tmp.a+len]=0;
			n_last++;
			len++;
		}
		else
		{
			break;
		}
	}

	for(int32_t i=1;i<=ed;i++)
	{
		for(int32_t j=0;j<n_last;j++)
		{
			if(H_last[j].a+1<x+1)
			{
				if(m[H_last[j].a+1][H_last[j].b]>i)
				{
					H_cur[n_cur].a=H_last[j].a+1;
					H_cur[n_cur].b=H_last[j].b;
					m[H_last[j].a+1][H_last[j].b]=i;
					n_cur++;
					len=1;
					tmp=H_cur[n_cur-1];
					while(1)
					{
						if(tmp.a+len>x||tmp.b+len>y)
						{
							break;
						}
						if(a[tmp.a+len-1]==b[tmp.b+len-1]&&m[tmp.a+len][tmp.b+len]>i)
						{
							H_cur[n_cur].a=tmp.a+len;
							H_cur[n_cur].b=tmp.b+len;
							m[tmp.a+len][tmp.b+len]=i;
							n_cur++;
							len++;
						}
						else
						{
							break;
						}
					}
				}
			}

			if(H_last[j].b+1<y+1)
			{
				if(m[H_last[j].a][H_last[j].b+1]>i)
				{
					H_cur[n_cur].a=H_last[j].a;
					H_cur[n_cur].b=H_last[j].b+1;
					m[H_last[j].a][H_last[j].b+1]=i;
					n_cur++;
					len=1;
					tmp=H_cur[n_cur-1];
					while(1)
					{
						if(tmp.a+len>x||tmp.b+len>y)
						{
							break;
						}
						if(a[tmp.a+len-1]==b[tmp.b+len-1]&&m[tmp.a+len][tmp.b+len]>i)
						{
							H_cur[n_cur].a=tmp.a+len;
							H_cur[n_cur].b=tmp.b+len;
							m[tmp.a+len][tmp.b+len]=i;
							n_cur++;
							len++;
						}
						else
						{
							break;
						}
					}
				}
			}

			if(H_last[j].b+1<y+1&&H_last[j].a+1<x+1)
			{
				if(m[H_last[j].a+1][H_last[j].b+1]>i)
				{
					H_cur[n_cur].a=H_last[j].a+1;
					H_cur[n_cur].b=H_last[j].b+1;
					m[H_last[j].a+1][H_last[j].b+1]=i;
					n_cur++;
					len=1;
					tmp=H_cur[n_cur-1];
					while(1)
					{
						if(tmp.a+len>x||tmp.b+len>y)
						{
							break;
						}
						if(a[tmp.a+len-1]==b[tmp.b+len-1]&&m[tmp.a+len][tmp.b+len]>i)
						{
							H_cur[n_cur].a=tmp.a+len;
							H_cur[n_cur].b=tmp.b+len;
							m[tmp.a+len][tmp.b+len]=i;
							n_cur++;
							len++;
						}
						else
						{
							break;
						}
					}
				}
			}
		}
		struct cell * cell_tmp;
		cell_tmp = H_last;
		H_last = H_cur;
		H_cur = cell_tmp ;
		n_last=n_cur;
		n_cur=0;

		if(m[x][y]<=ed)
		{
			r = m[x][y];
			break;
		}

	}

//	for(uint32_t i=0;i<x+1;i++)
//	{
//		free(m[i]);
//	}
//	free(m);

	free(p_last);
	free(p_cur);
	return m;
}

typedef struct range
{
	int start;
	int end;
	int len;

}RANGE;
int edit_range(string str1, string str2,int tor)
{
    int max1 = str1.size();
    int max2 = str2.size();
	if(max1==0 || max2==0)
	{
		return max1+max2;
	}

	if(abs(max1-max2)>tor)
	{
		return tor+1;
	}

    int **ptr = new int*[max1 + 1];
    for(int i = 0; i < max1 + 1 ;i++)
    {
        ptr[i] = new int[max2 + 1];
    }

    for(int i = 0 ;i < max1 + 1 ;i++)
    {
        ptr[i][0] = i;
    }

    for(int i = 0 ;i < max2 + 1;i++)
    {
        ptr[0][i] = i;
    }

    for(int i = 1 ;i < max1 + 1 ;i++)
    {
		int b=1;
		int lb=max(1,i-tor);
		int ub=min(max2,i+tor);
        for(int j = lb ;j<= ub; j++)
        {
			int temp;
			if(tor>0)
			{
				if(i-j==tor)
				{
					temp = ptr[i-1][j] + 1;
				}
				else
				{
					if(i-j==-tor)
					{
						temp = ptr[i][j-1] + 1;
					}
					else
					{
						temp = min(ptr[i-1][j] + 1, ptr[i][j-1] + 1);
					}
				}
			}
			else
			{
				temp=tor+1;
			}

            int d;
            if(str1[i-1] == str2[j-1])
            {
                d = 0 ;
            }
            else
            {
                d = 1 ;
            }
            ptr[i][j] = min(temp, ptr[i-1][j-1] + d);
			if(ptr[i][j]<=tor)
			{
				b=0;
			}
        }
		if(b==1)
		{
			for(int i = 0; i < max1 + 1; i++)
			{
				delete[] ptr[i];
				ptr[i] = NULL;
			}

			delete[] ptr;
			ptr = NULL;

			return tor+1;
		}
    }
	/*cout << "**************************" << endl;
    for(int i = 0 ;i < max1 + 1 ;i++)
    {
        for(int j = 0; j< max2 + 1; j++)
        {
            cout << ptr[i][j] << "	" ;
        }
        cout << endl;
    }
    cout << "**************************" << endl;*/

 	int dis = ptr[max1][max2];

	if(dis<0)
	{
		cout << "dis xiaoyu 0!" << endl;
		for(int i = 0 ;i < max1 + 1 ;i++)
		{
			for(int j = 0; j< max2 + 1; j++)
			{
				cout << ptr[i][j] << "	" ;
			}
			cout << endl;
		}
	}

    for(int i = 0; i < max1 + 1; i++)
    {
        delete[] ptr[i];
        ptr[i] = NULL;
    }

    delete[] ptr;
    ptr = NULL;


    return dis;
}
bool checkPos(string s, int pos, string seg)
{
	int l=1;
	for(int i=0;i<seg.size();i++)
	{
		if(seg[i]!=s[pos+i])
		{
			l=0;
		}
		if(i+pos==s.size())
		{
			l=0;
		}
	}
	return l;
}
void recur_ed_divide(string q,int tor,vector <RANGE> &ranges,int start)
{
	int N_tor=tor+1;
	int q_a=q.size()/N_tor;
	int q_b=q.size()%N_tor;

	vector <int> echSegLen;

	RANGE range_tmp;

	for(int i=0;i<q_b;i++)
	{
		echSegLen.push_back(q_a+1);
		range_tmp.len=q_a+1;
		ranges.push_back(range_tmp);
	}
	for(int i=q_b;i<N_tor;i++)
	{
		echSegLen.push_back(q_a);
		range_tmp.len=q_a;
		ranges.push_back(range_tmp);
	}

	for(int i=1;i<N_tor;i++)
	{
		echSegLen[i]=echSegLen[i]+echSegLen[i-1];
	}


	ranges[0].start=0+start;
	ranges[0].end=echSegLen[0]-1+start;

	for(int i=1;i<N_tor;i++)
	{
		ranges[i].start=echSegLen[i-1]+start;
		ranges[i].end=echSegLen[i]-1+start;
	}
}
int recur_ed_without_index(string s, string r, RANGE range_s, RANGE range_r,int tor,int lb, int ub,int level)
{
	level+=1;
	if(abs(range_s.len-range_r.len)>tor)
	{
		return tor+1;
	}

	vector<RANGE> RR;
	recur_ed_divide(r.substr(range_r.start,range_r.len),tor,RR,range_r.start);

	int ed=tor+1;

	for(int i=0;i<tor+1;i++)
	{
		int lowerBound,upperBound;
		lowerBound=max(lb+RR[i].start,range_s.start);
		upperBound=min(ub+RR[i].start, range_s.end-RR[i].len+1);

		for(int j=lowerBound;j<=upperBound;j++)
		{
			if(checkPos(s,j,r.substr(RR[i].start,RR[i].len)))
			{
				int edmin;

				RANGE sl,sr;
				RANGE rl,rr;

				sl.start=range_s.start;
				sl.end=j-1;
				sl.len=sl.end-sl.start+1;

				rl.start=range_r.start;
				rl.end=RR[i].start-1;
				rl.len=rl.end-rl.start+1;

				sr.start=j+RR[i].len;
				sr.end=range_s.end;
				sr.len=sr.end-sr.start+1;

				rr.start=RR[i].end+1;
				rr.end=range_r.end;
				rr.len=rr.end-rr.start+1;

				if(sl.len<max(range_s.end/2,50) || rl.len<max(range_s.end/2,50))
				{
					edmin=edit_range(s.substr(sl.start,sl.len),r.substr(rl.start,rl.len),tor-abs(sr.len-rr.len));
				}
				else
				{
					edmin=recur_ed_without_index(s,r,sl,rl,tor-abs(sr.len-rr.len),lb,ub,level);
				}

				if(edmin<=tor-abs(sr.len-rr.len))
				{
					if(sr.len<max(range_s.end/2,50) || rr.len<max(range_s.end/2,50))
					{
						edmin+=edit_range(s.substr(sr.start,sr.len),r.substr(rr.start,rr.len),tor-abs(sl.len-rl.len));
					}
					else
					{
						edmin+=recur_ed_without_index(s,r,sr,rr,tor-abs(sl.len-rl.len),lb,ub,level);
					}

					if(edmin<=tor && level==1)
					{
						return edmin;
					}
					if(edmin<ed)
					{
						ed=edmin;
					}
				}
			}
		}
	}
	return ed;
}
int recur_ed(string s,string r,int tor)
{
	int diff=s.size()-r.size();
	if((diff > tor) || (diff < -tor))
	{
		return tor+1;
	}
	int deta=(tor-abs(diff))/2;
	int lowerBound, upperBound;
	RANGE rangeS,rangeR;
	if(diff>=0)
	{
		lowerBound=-deta;
		upperBound=diff+deta;
	}
	else
	{
		lowerBound=diff-deta;
		upperBound=deta;
	}
	rangeS.start=0;
	rangeS.end=s.size()-1;
	rangeS.len=s.size();

	rangeR.start=0;
	rangeR.end=r.size()-1;
	rangeR.len=r.size();

	int ed=recur_ed_without_index(s,r,rangeS,rangeR,tor,lowerBound,upperBound,0);

	return ed;
}

