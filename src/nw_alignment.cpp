#include "util.h"


#define MINUS_INF -0x3fffffff


typedef struct {
	int32_t h, e;
} eh_t;

#define LIKELY(x) __builtin_expect((x), 1) //gcc内置函数, 帮助编译器分支优化








const float MaxPenalty = -65536;
const float OPEN_GAP = -1;
const float EXTEND_GAP = -0.5;
const float NEW_GAP = -1.5;

double max(float x, float y)
{
	return x > y ? x : y;
}

double max(float x, float y, float z)
{
	return x > y ? max(x, z) : max(y, z);
}
//m是read
//n是ref
void nw_alignment(int m, string& s1, int n, string& s2)
{

			printf("%s\n",s1.c_str());

	printf("%s\n",s2.c_str());
	int i, j;

	m = m + 1, n = n + 1;

	float** r = new float*[m];
	float** t = new float*[m];
	float** s = new float*[m];

	for (i = 0; i < m; i++)
	{
		r[i] = new float[n];
		t[i] = new float[n];
		s[i] = new float[n];
	}

	// initialization
	r[0][0] = t[0][0] = s[0][0] = 0;
	for (i = 1; i < m; i++)
	{
		r[i][0] = MaxPenalty;
		s[i][0] = t[i][0] = OPEN_GAP + i*EXTEND_GAP;
	}

	for (j = 1; j < n; j++)
	{
		t[0][j] = MaxPenalty;
		s[0][j] = r[0][j] = OPEN_GAP + j*EXTEND_GAP;
	}

	for (i = 1; i < m; i++)
	{
		for (j = 1; j < n; j++)
		{
			r[i][j] = max(r[i][j - 1] + EXTEND_GAP, s[i][j - 1] + NEW_GAP);
			t[i][j] = max(t[i - 1][j] + EXTEND_GAP, s[i - 1][j] + NEW_GAP);
			s[i][j] = max(s[i - 1][j - 1] + (nst_nt4_table[(uint8_t)s1[i - 1]] == nst_nt4_table[(uint8_t)s2[j - 1]] ? 1.5 : -1.5), r[i][j], t[i][j]);
		}
	}
	// back tracking




	i = m - 1, j = n - 1;
	while (i > 0 || j > 0) {
		//来自左边
		if (s[i][j] == r[i][j]) {
			s1.insert(i, 1, '-');
			j--;
		}
		//来自上方
		else if (s[i][j] == t[i][j]) {
			s2.insert(j, 1, '-');
			i--;
		}
		else {
			i--, j--;
		}
	}

		printf("%s\n",s1.c_str());

	printf("%s\n",s2.c_str());
	for (i = 0; i < m; i++)
	{
		delete[] r[i]; 
		delete[] t[i]; 
		delete[] s[i];
	}
	delete[] r; delete[] t; delete[] s;

	
}

// //这里的m是ref，s1是参考序列的
// //这里的n的read，s2是read的
// void nw_alignment_band(int m, string& s1, int n, string& s2)
// {

// 	printf("%s\n",s1.c_str());

// 	printf("%s\n",s2.c_str());
	

// 	int i, j,w=5,qlen=n;
// 	float h1;

// 	m = m + 1, n = n + 1;


// 	float* t = new float[n];
// 	float* s = new float[n];

// 	int n_col = qlen < 2*w+1? qlen : 2*w+1;
// 	int** diraction = new int*[m];

// 	for (i = 0; i < m; i++){
// 		diraction[i] = new int[n_col];
// 	}

// 	// initialization
// 	float r = MINUS_INF;
// 	t[0] = MaxPenalty;
// 	s[0] = 0;


// 	for (j = 1; j < n && j <= w; j++){
// 		t[j] = MaxPenalty;
// 		s[j] = OPEN_GAP + j*EXTEND_GAP;
// 	}
// 	for (; j <= n; ++j){
// 		s[j] =t[j] = MaxPenalty;
// 	}


// 	for(int x =0;x<n;x++){
// 		printf("%f ",x,s[x]);
// 	}
// 	printf("\n");

// 	for (i = 1; i < m; i++){
// 		int32_t  beg, end;
// 		float f = MINUS_INF ,h1, temp;
// 		beg = i > w? i - w : 0;
// 		end = i + w + 1 < qlen? i + w + 1 : qlen; // only loop through [beg,end) of the query sequence
// 		h1 = beg == 0? (OPEN_GAP + EXTEND_GAP * (i + 1)) : MINUS_INF;
// 		for(int y=0;y<beg;y++){
// 			printf(" ");
// 		}

// 		printf("%d-%d\t", beg,end);
// 		for (j = beg; j < end; ++j){
// 				// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
// 				// Cells are computed in the following order:
// 				//   M(i,j)   = H(i-1,j-1) + S(i,j)
// 				//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
// 				//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
// 				//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
// 				float h, M = s[j], E = t[j];
// 				uint8_t D; // direction
// 				s[j] = h1;
// 				//   M(i,j)   = H(i-1,j-1) + S(i,j)
// 				M += (nst_nt4_table[(uint8_t)s1[i - 1]] == nst_nt4_table[(uint8_t)s2[j - 1]] ? 1.5 : -1.5);
// 				//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}				
// 				D = M >= E? 0 : 1;
// 				h = M >= E? M : E;
// 				D = h >= f? D : 2;
// 				h = h >= f? h : f;
// 				h1 = h;
// 				// E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
// 				temp = M - OPEN_GAP;
// 				E -= EXTEND_GAP;
// 				E  = E > temp? E    : temp;
// 				t[j] = E;
// 				temp = m - OPEN_GAP;
// 				f -= EXTEND_GAP;
// 				f  = f > temp? f    : temp;
// 				diraction[i][j - beg] = D; // z[i,j] keeps h for the current cell and e/f for the next cell
// 				printf("%d ",D);
// 		}
// 		printf("\n");
// 		s[end] = h1; t[end] = MaxPenalty;
// 	}

// 	float score = s[n-1];
// 	// back tracking



// 	i = m - 1; 

// 	int k = n_col - 1;

// 	printf("%d-%d\n",i,k);
// 	int which =0;
// 	while (i >= 0 && k >= 0) {
// 		which = diraction[i][(k)];
// 		//左上回溯
// 		if (which == 0)  {
// 			--i, --k;
// 		}    
// 		//向上回溯
// 		else if (which == 1){
// 			s2.insert(j, 1, '-');
// 			i--;
// 		} else  {//向左回溯
// 			--k;
// 			s1.insert(i, 1, '-');
// 		}               
// 	}

// 	for (i = 0; i < m; i++)
// 	{
// 		delete[] diraction[i]; 
// 	}
// 	 delete[] t; delete[] s;

// 	printf("%s\n",s1.c_str());

// 	printf("%s\n",s2.c_str());

// }

// void nw_alignment_band(int m, string& s1, int n, string& s2)
// {
// 	int i, j,w=20,qlen =n,beg,end;

// 	m = m + 1, n = n + 1;

// 	float** f = new float*[m];
// 	float** e = new float*[m];
// 	float** h = new float*[m];

// 	for (i = 0; i < m; i++)
// 	{
// 		f[i] = new float[n];
// 		e[i] = new float[n];
// 		h[i] = new float[n];
// 	}

// 	// initialization
// 	f[0][0] = e[0][0] = h[0][0] = 0;
// 	for (i = 1; i < m; i++){
// 		f[i][0] = MaxPenalty;
// 		h[i][0] = e[i][0] = OPEN_GAP + i*EXTEND_GAP;
// 	}

// 	for (j = 1; j <= qlen && j <= w; j++)
// 	{
// 		e[0][j] = MaxPenalty;
// 		h[0][j] = f[0][j] = OPEN_GAP + j*EXTEND_GAP;
// 	}

// // 	for (j = 1; j <= qlen && j <= w; ++j)
// // 		eh[j].h = -(o_ins + e_ins * j), eh[j].e = MINUS_INF;
// 	for (; j <= qlen; ++j) {
// 		h[0][j] = e[0][j] = MaxPenalty; // everything is -inf outside the band
// 	}


// 	for (i = 1; i < m; i++)
// 	{
// 		beg = i > w? i - w : 0;
// 		end = i + w + 1 < qlen? i + w + 1 : qlen; // only loop through [beg,end) of the query sequence
// 		for (j = beg; j < end; j++)
// 		{
// 			f[i][j] = max(f[i][j - 1] + EXTEND_GAP, h[i][j - 1] + NEW_GAP);
// 			e[i][j] = max(e[i - 1][j] + EXTEND_GAP, h[i - 1][j] + NEW_GAP);
// 			h[i][j] = max(h[i - 1][j - 1] + (nst_nt4_table[(uint8_t)s1[i - 1]] == nst_nt4_table[(uint8_t)s2[j - 1]] ? 1.5 : -1.5), f[i][j], e[i][j]);
// 		}
// 	}
// 	// back tracking

	
// 	i = m - 1, j = (i + w + 1 < qlen? i + w + 1 : qlen) - 1;
// 	while (i > 0 || j > 0) {
// 		//来自左边
// 		if (h[i][j] == f[i][j]) {
// 			s1.insert(i, 1, '-');
// 			j--;
// 		}
// 		//来自上方
// 		else if (h[i][j] == e[i][j]) {
// 			s2.insert(j, 1, '-');
// 			i--;
// 		}
// 		else {
// 			i--, j--;
// 		}
// 	}
// 	for (i = 0; i < m; i++)
// 	{
// 		delete[] f[i]; 
// 		delete[] e[i]; 
// 		delete[] h[i];
// 	}
// 	delete[] f; delete[] e; delete[] h;
// }

void fillScoreMatrix(int a, int b, int8_t mat[25])
{
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j){
			mat[k++] = i == j? a : -b;
		}
			
		mat[k++] = -1; // ambiguous base
	}
	for (j = 0; j < 5; ++j) {
		mat[k++] = -1;
	}
}

// static inline uint32_t *push_cigar(int *n_cigar, int *m_cigar, uint32_t *cigar, int op, int len)
// {
// 	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
// 		if (*n_cigar == *m_cigar) {
// 			*m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
// 			cigar = (uint32_t *)realloc(cigar, (*m_cigar) << 2);
// 		}
// 		cigar[(*n_cigar)++] = len<<4 | op;
// 	} else cigar[(*n_cigar)-1] += len<<4;
// 	return cigar;
// }

void nw_alignment_band(int qlen, string& query, int tlen, string& target)
{

	// printf("%s\n",query.c_str());

	// printf("%s\n",target.c_str());
	int o_del =2,  e_del=1,  o_ins=2,  e_ins=1,  w=10;
	eh_t *eh;

	int i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, score, n_col;
	uint8_t *z; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be a little more complex

	// allocate memory
	n_col = qlen < 2*w+1? qlen : 2*w+1; // maximum #columns of the backtrack matrix
	z = (uint8_t *) malloc((long)n_col * tlen) ;
	eh = (eh_t *)calloc(qlen + 1, 8);
	// generate the query profile

	// fill the first row
	eh[0].h = 0; 
	eh[0].e = MINUS_INF;
	for (j = 1; j <= qlen && j <= w; ++j)
		eh[j].h = -(o_ins + e_ins * j), eh[j].e = MINUS_INF;
	for (; j <= qlen; ++j) eh[j].h = eh[j].e = MINUS_INF; // everything is -inf outside the band
	// DP loop
	for (i = 0; LIKELY(i < tlen); ++i) { // target sequence is in the outer loop
		int32_t f = MINUS_INF, h1, beg, end, t;

		beg = i > w? i - w : 0;
		end = i + w + 1 < qlen? i + w + 1 : qlen; // only loop through [beg,end) of the query sequence
		h1 = beg == 0? -(o_del + e_del * (i + 1)) : MINUS_INF;
		
			uint8_t *zi = &z[(long)i * n_col];
			for (j = beg; LIKELY(j < end); ++j) {
				// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
				// Cells are computed in the following order:
				//   M(i,j)   = H(i-1,j-1) + S(i,j)
				//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
				//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
				//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
				// We have to separate M(i,j); otherwise the direction may not be recorded correctly.
				// However, a CIGAR like "10M3I3D10M" allowed by local() is disallowed by global().
				// Such a CIGAR may occur, in theory, if mismatch_penalty > 2*gap_ext_penalty + 2*gap_open_penalty/k.
				// In practice, this should happen very rarely given a reasonable scoring system.
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				uint8_t d; // direction
				p->h = h1;
				m += (nst_nt4_table[(uint8_t)target[i - 1]] == nst_nt4_table[(uint8_t)query[j - 1]] ? 1 : -2);
				d = m >= e? 0 : 1;
				h = m >= e? m : e;
				d = h >= f? d : 2;
				h = h >= f? h : f;
				h1 = h;
				t = m - oe_del;
				e -= e_del;
				d |= e > t? 1<<2 : 0;
				e  = e > t? e    : t;
				p->e = e;
				t = m - oe_ins;
				f -= e_ins;
				d |= f > t? 2<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f > t? f    : t;
				zi[j - beg] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}

		
		eh[end].h = h1; eh[end].e = MINUS_INF;
	}
	score = eh[qlen].h;
	 // backtrack
		int n_cigar = 0, m_cigar = 0, which = 0,insert_pos;
		uint32_t *cigar = 0, tmp;
		i = tlen - 1; k = (i + w + 1 < qlen? i + w + 1 : qlen) - 1; // (i,k) points to the last cell
		while (i >= 0 && k >= 0) {
			which = z[(long)i * n_col + (k - (i > w? i - w : 0))] >> (which<<1) & 3;
			if (which == 0)  {
				--i, --k;
			}    
			else if (which == 1){
				insert_pos =i>qlen ? qlen :i;
				// printf("insert_pos:%d qlen:%d\n",insert_pos,qlen);
				query.insert(insert_pos, 1, '-');
				 --i;
			} 
			else  {
				insert_pos =k>tlen ? tlen :k;
				// printf("insert_pos:%d tlen:%d\n",insert_pos,tlen);
				target.insert(insert_pos, 1, '-');
				 --k;
			}              
		}
		if (i >= 0){
			query.insert(0, i, '-');
		} 
		if (k >= 0){
			target.insert(0, k, '-');
		} 

	free(eh); 
	free(z);
	// return score;
}