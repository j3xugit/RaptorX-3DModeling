#include "cdhit_seqdb.h"


//=========== inner function ============//
////For smiple len1 <= len2, len2 is for existing representative
////walk along all diag path of two sequences,
////find the diags with most aap
////return top n diags
////added on 2006 11 13
////band 0                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                            XXXXXXXXXXXXXXX                  seq1
////band 1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                             XXXXXXXXXXXXXXX                 seq1
////extreme right (+)           XXXXXXXXXXXXXXXXXX               seq2, rep seq
////    band = len2-1                            XXXXXXXXXXXXXXX seq1
////band-1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                           XXXXXXXXXXXXXXX                   seq1
////extreme left (-)            XXXXXXXXXXXXXXXXXX               seq2, rep seq
////              XXXXXXXXXXXXXXX   band = -(len1-1)             seq1
////index of diag_score = band+len1-1;
int SequenceDB::diag_test_aapn(int NAA1, char iseq2[], int len1, int len2, WorkingBuffer & buffer, const Options & options,
		int &best_sum, int band_width, int &band_left, int &band_center, int &band_right, int required_aa1)
{
	int i, i1, j, k;
	int *pp;
	int nall = len1+len2-1;
	Vector<int> & taap = buffer.taap;
	Vector<INTs> & aap_begin = buffer.aap_begin;
	Vector<INTs> & aap_list = buffer.aap_list;
	Vector<int> & diag_score = buffer.diag_score;
	Vector<int> & diag_score2 = buffer.diag_score2;

	if (nall > MAX_DIAG) CDHIT_bomb_error("in diag_test_aapn, MAX_DIAG reached");
	for (pp=&diag_score[0], i=nall; i; i--, pp++) *pp=0;
	for (pp=&diag_score2[0], i=nall; i; i--, pp++) *pp=0;

	int c22, cpx;
	int len11 = len1-1;
	int len22 = len2-1;
	i1 = len11;
	for (i=0; i<len22; i++,i1++) {
		c22 = iseq2[i]*NAA1+ iseq2[i+1];
		cpx = 1 + (iseq2[i] != iseq2[i+1]);
		if ( (j=taap[c22]) == 0) continue;
		int m = aap_begin[c22];
		for(int k=0; k<j; k++){
			diag_score[ i1 - aap_list[m+k] ] ++;
			diag_score2[ i1 - aap_list[m+k] ] += cpx;
		}
	}

	//find the best band range
	//  int band_b = required_aa1;
	int band_b = required_aa1-1 >= 0 ? required_aa1-1:0;  // on dec 21 2001
	int band_e = nall - required_aa1;
	int band_m = ( band_b+band_width-1 < band_e ) ? band_b+band_width-1 : band_e;
	int best_score=0;
	int best_score2=0;
	int max_diag = 0;
	int max_diag2 = 0;
	int imax_diag = 0;
	for (i=band_b; i<=band_m; i++){
		best_score += diag_score[i];
		best_score2 += diag_score2[i];
		if( diag_score2[i] > max_diag2 ){
			max_diag2 = diag_score2[i];
			max_diag = diag_score[i];
			imax_diag = i;
		}
	}
	int from=band_b;
	int end =band_m;
	int score = best_score;
	int score2 = best_score2;
	for (k=from, j=band_m+1; j<band_e; j++, k++) {
		score -= diag_score[k]; 
		score += diag_score[j]; 
		score2 -= diag_score2[k]; 
		score2 += diag_score2[j]; 
		if ( score2 > best_score2 ) {
			from = k + 1;
			end  = j;
			best_score = score;
			best_score2 = score2;
			if( diag_score2[j] > max_diag2 ){
				max_diag2 = diag_score2[j];
				max_diag = diag_score[j];
				imax_diag = j;
			}
		}
	}
	int mlen = imax_diag;
	if( imax_diag > len1 ) mlen = nall - imax_diag;
	int emax = int((1.0 - options.cluster_thd) * mlen) + 1;
	for (j=from; j<imax_diag; j++) { // if aap pairs fail to open gap
		if ( (imax_diag - j) > emax || diag_score[j] < 1 ) {
			best_score -= diag_score[j]; from++;
		} else break;
	}
	for (j=end; j>imax_diag; j--) { // if aap pairs fail to open gap
		if ( (j - imax_diag) > emax || diag_score[j] < 1 ) {
			best_score -= diag_score[j]; end--;
		} else break;
	}

	//  delete [] diag_score;
	band_left = from - len1 + 1; 
	band_right= end - len1 + 1;
	band_center = imax_diag - len1 + 1;
	best_sum = best_score;
	return OK_FUNC;
}
// END diag_test_aapn
 

int SequenceDB::diag_test_aapn_est(int NAA1, char iseq2[], int len1, int len2, WorkingBuffer & buffer, const Options & options,
        int &best_sum, int band_width, int &band_left, int &band_center, int &band_right, int required_aa1)
{
	int i, i1, j, k;
	int *pp, *pp2;
	int nall = len1+len2-1;
	int NAA2 = NAA1 * NAA1;
	int NAA3 = NAA2 * NAA1;
	Vector<int> & taap = buffer.taap;
	Vector<INTs> & aap_begin = buffer.aap_begin;
	Vector<INTs> & aap_list = buffer.aap_list;
	Vector<int> & diag_score = buffer.diag_score;
	Vector<int> & diag_score2 = buffer.diag_score2;

	if (nall > MAX_DIAG) CDHIT_bomb_error("in diag_test_aapn_est, MAX_DIAG reached");
	pp = & diag_score[0];
	pp2 = & diag_score2[0];
	for (i=nall; i; i--, pp++, pp2++) *pp = *pp2 =0;

	INTs *bip;
	int c22, cpx;
	int len22 = len2-3;
	i1 = len1-1;
	for (i=0; i<len22; i++,i1++,iseq2++) {
		unsigned char c0 = iseq2[0];
		unsigned char c1 = iseq2[1];
		unsigned char c2 = iseq2[2];
		unsigned char c3 = iseq2[3];
		if ((c0==4) || (c1==4) || (c2==4) || (c3==4)) continue; //skip N

		c22 = c0*NAA3+ c1*NAA2 + c2*NAA1 + c3;
		if ( (j=taap[c22]) == 0) continue;
		cpx = 1 + (c0 != c1) + (c1 != c2) + (c2 != c3);
		bip = & aap_list[ aap_begin[c22] ];     //    bi = aap_begin[c22];
		for (; j; j--, bip++) { 
			diag_score[i1 - *bip]++;
			diag_score2[i1 - *bip] += cpx;
		}
	}
#if 0
	int mmax = 0;
	int immax = 0;
	for(i=0; i<=nall; i++){
		if( i%len2 ==0 or i == nall ) printf( "\n" );
		printf( "%3i ", diag_score[i] );
		if( diag_score[i] > mmax ){
			mmax = diag_score[i];
			immax = i;
		}
	}
#endif

	//find the best band range
	//  int band_b = required_aa1;
	int band_b = required_aa1-1 >= 0 ? required_aa1-1:0;  // on dec 21 2001
	int band_e = nall - required_aa1;
	int band_m = ( band_b+band_width-1 < band_e ) ? band_b+band_width-1 : band_e;
	int best_score=0;
	int best_score2=0;
	int max_diag = 0;
	int max_diag2 = 0;
	int imax_diag = 0;
	for (i=band_b; i<=band_m; i++){
		best_score += diag_score[i];
		best_score2 += diag_score2[i];
		if( diag_score2[i] > max_diag2 ){
			max_diag2 = diag_score2[i];
			max_diag = diag_score[i];
			imax_diag = i;
		}
	}
	int from=band_b;
	int end =band_m;
	int score = best_score;  
	int score2 = best_score2;  

	for (k=from, j=band_m+1; j<band_e; j++, k++) {
		score -= diag_score[k]; 
		score += diag_score[j]; 
		score2 -= diag_score2[k]; 
		score2 += diag_score2[j]; 
		if ( score2 > best_score2 ) {
			from = k + 1;
			end  = j;
			best_score = score;
			best_score2 = score2;
			if( diag_score2[j] > max_diag2 ){
				max_diag2 = diag_score2[j];
				max_diag = diag_score[j];
				imax_diag = j;
			}
		}
	}
#if 0
	printf( "%i\n", required_aa1 );
	printf( "max=%3i  imax=%3i; band:  %3i  %3i  %i\n", max_diag, imax_diag, band_b, band_e, band_m );
	printf( "best: %i\n", best_score );
	printf( "from: %i, end: %i,  best: %i\n", from, end, best_score );
#endif
	int mlen = imax_diag;
	if( imax_diag > len1 ) mlen = nall - imax_diag;
	int emax = int((1.0 - options.cluster_thd) * mlen) + 1;
	for (j=from; j<imax_diag; j++) { // if aap pairs fail to open gap
		if ( (imax_diag - j) > emax || diag_score[j] < 1 ) {
			best_score -= diag_score[j]; from++;
		} else break;
	}
	for (j=end; j>imax_diag; j--) { // if aap pairs fail to open gap
		if ( (j - imax_diag) > emax || diag_score[j] < 1 ) {
			best_score -= diag_score[j]; end--;
		} else break;
	}
	//printf( "best: %i\n", best_score );

	band_left = from-len1+1; 
	band_right= end-len1+1;
	band_center = imax_diag - len1 + 1;
	best_sum = best_score;
#if 0
	printf( "max=%3i  imax=%3i; band:  %3i  %3i  %i\n", mmax, immax, band_b, band_e, band_m );
	printf( "%3i:  best: %i,  %i  %i  %i\n", required_aa1, best_score, band_left, band_right, band_width );
#endif
	return OK_FUNC;
}
// END diag_test_aapn_est



/*
local alignment of two sequence within a diag band
for band 0 means direction (0,0) -> (1,1)
         1 means direction (0,1) -> (1,2)
        -1 means direction (1,0) -> (2,1)
added on 2006 11 13
band 0                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                            XXXXXXXXXXXXXXX                  seq1
band 1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                             XXXXXXXXXXXXXXX                 seq1
extreme right (+)           XXXXXXXXXXXXXXXXXX               seq2, rep seq
    band = len2-1                            XXXXXXXXXXXXXXX seq1
band-1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                           XXXXXXXXXXXXXXX                   seq1
extreme left (-)            XXXXXXXXXXXXXXXXXX               seq2, rep seq
              XXXXXXXXXXXXXXX   band = -(len1-1)             seq1
iseq len are integer sequence and its length,
mat is matrix, return ALN_PAIR class

       band:  -101   seq2 len2 = 17
                \\\1234567890123456
              0  \xxxxxxxxxxxxxxxxx
              1   xxxxxxxxxxxxxxxxx\ most right band = len2-1
              2   xxxxxxxxxxxxxxxxx
    seq1      3   xxxxxxxxxxxxxxxxx
    len1 = 11 4   xxxxxxxxxxxxxxxxx
              5   xxxxxxxxxxxxxxxxx
              6   xxxxxxxxxxxxxxxxx
              7   xxxxxxxxxxxxxxxxx
              8   xxxxxxxxxxxxxxxxx
              9   xxxxxxxxxxxxxxxxx
              0   xxxxxxxxxxxxxxxxx
                  \
                   most left band = -(len1-1)

*/

int SequenceDB::local_band_align( char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix &mat, 
		int &best_score, int &iden_no, int &alnln, float &dist, int *alninfo,
		int band_left, int band_center, int band_right, WorkingBuffer & buffer, const Options & options)
{
	int i, j, k, j1;
	int64_t best_score1;
	iden_no = 0;

	if ( (band_right >= len2 ) ||
			(band_left  <= -len1) ||
			(band_left  > band_right) ) return FAILED_FUNC;

	// allocate mem for score_mat[len1][len2] etc
	int band_width = band_right - band_left + 1;
	int band_width1 = band_width + 1;

	MatrixInt64 & score_mat = buffer.score_mat;
	MatrixInt   & back_mat = buffer.back_mat;

	//printf( "%i  %i\n", band_right, band_left );

	if( (int)score_mat.size() <= len1 ){
		VectorInt   row( band_width1, 0 );
		VectorInt64 row2( band_width1, 0 );
		while( (int)score_mat.size() <= len1 ){
			score_mat.Append( row2 );
			back_mat.Append( row );
		}
	}
	for(i=0; i<=len1; i++){
		if( score_mat[i].Size() < band_width1 ) score_mat[i].Resize( band_width1 );
		if( back_mat[i].Size() < band_width1 ) back_mat[i].Resize( band_width1 );
	}

	best_score = 0;
	/*
	   seq2 len2 = 17            seq2 len2 = 17      seq2 len2 = 17
	   01234567890123456       01234567890123456    01234567890123456
	   0     xxxxxxxxxxxxxxxxx \\\\\\XXXxxxxxxxxxxxxxx    xXXXXXXXxxxxxxxxx
	   1\\\\\Xxxxxxxxxxxxxxxxx  \\\\\Xxx\xxxxxxxxxxxxx    xx\xxxxx\xxxxxxxx
	   2 \\\\X\xxxxxxxxxxxxxxx   \\\\Xxxx\xxxxxxxxxxxx    xxx\xxxxx\xxxxxxx
	   seq1 3  \\\Xx\xxxxxxxxxxxxxx    \\\Xxxxx\xxxxxxxxxxx    xxxx\xxxxx\xxxxxx
	   len1 4   \\Xxx\xxxxxxxxxxxxx     \\Xxxxxx\xxxxxxxxxx    xxxxx\xxxxx\xxxxx
	   = 11 5    \Xxxx\xxxxxxxxxxxx      \Xxxxxxx\xxxxxxxxx    xxxxxx\xxxxx\xxxx
	   6     Xxxxx\xxxxxxxxxxx       Xxxxxxxx\xxxxxxxx    xxxxxxx\xxxxx\xxx
	   7     x\xxxx\xxxxxxxxxx       x\xxxxxxx\xxxxxxx    xxxxxxxx\xxxxx\xx
	   8     xx\xxxx\xxxxxxxxx       xx\xxxxxxx\xxxxxx    xxxxxxxxx\xxxxx\x
	   9     xxx\xxxx\xxxxxxxx       xxx\xxxxxxx\xxxxx    xxxxxxxxxx\xxxxx\
	   0     xxxx\xxxx\xxxxxxx       xxxx\xxxxxxx\xxxx    xxxxxxxxxxx\xxxxx
	   band_left < 0           band_left < 0        band_left >=0
	   band_right < 0          band_right >=0       band_right >=0
	   init score_mat, and iden_mat (place with upper 'X')
	 */

	if (band_left < 0) {  //set score to left border of the matrix within band
		int tband = (band_right < 0) ? band_right : 0;
		//for (k=band_left; k<tband; k++)
		for (k=band_left; k<=tband; k++) { // fixed on 2006 11 14
			i = -k;
			j1 = k-band_left;
			// penalty for leading gap opening = penalty for gap extension
			score_mat[i][j1] =  mat.ext_gap * i;
			back_mat[i][j1] = DP_BACK_TOP;
		}
		back_mat[-tband][tband-band_left] = DP_BACK_NONE;
	}

	if (band_right >=0) { //set score to top border of the matrix within band
		int tband = (band_left > 0) ? band_left : 0;
		for (j=tband; j<=band_right; j++) {
			j1 = j-band_left;
			score_mat[0][j1] = mat.ext_gap * j;
			back_mat[0][j1] = DP_BACK_LEFT;
		}
		back_mat[0][tband-band_left] = DP_BACK_NONE;
	}

	int gap_open[2] = { mat.gap, mat.ext_gap };
	int max_diag = band_center - band_left;
	int extra_score[4] = { 4, 3, 2, 1 };
	for (i=1; i<=len1; i++) {
		int J0 = 1 - band_left - i;
		int J1 = len2 - band_left - i;
		if( J0 < 0 ) J0 = 0;
		if( J1 >= band_width ) J1 = band_width;
		for (j1=J0; j1<=J1; j1++){
			j = j1+i+band_left;

			int ci = iseq1[i-1];
			int cj = iseq2[j-1];
			int sij = mat.matrix[ci][cj];
			//int iden_ij = (ci == cj);
			int back;

			/* extra score according to the distance to the best diagonal */
			int extra = extra_score[ abs(j1 - max_diag) & 3 ]; // max distance 3
			sij += extra * (sij>0);

			back = DP_BACK_LEFT_TOP;
			best_score1 = score_mat[i-1][j1] + sij;
			int gap0 = gap_open[ (i == len1) | (j == len2) ];
			int gap = 0;
			int64_t score;

			if( j1 > 0 ){
				gap = gap0;
				if( back_mat[i][j1-1] == DP_BACK_LEFT ) gap = mat.ext_gap;
				if( (score = score_mat[i][j1-1] + gap) > best_score1 ){
					back = DP_BACK_LEFT;
					best_score1 = score;
				}
			}
			if(j1+1<band_width){
				gap = gap0;
				if( back_mat[i-1][j1+1] == DP_BACK_TOP ) gap = mat.ext_gap;
				if( (score = score_mat[i-1][j1+1] + gap) > best_score1 ){
					back = DP_BACK_TOP;
					best_score1 = score;
				}
			}
			score_mat[i][j1] = best_score1;
			back_mat[i][j1]  = back;
			//printf( "%2i(%2i) ", best_score1, iden_no1 );

		}
		//printf( "\n" );
	}
	i = j = 0;
	if( len2 - band_left < len1 ){
		i = len2 - band_left;
		j = len2;
	}else if( len1 + band_right < len2 ){
		i = len1;
		j = len1 + band_right;
	}else{
		i = len1;
		j = len2;
	}
	j1 = j - i - band_left;
	best_score = (int)score_mat[i][j1];
	best_score1 = score_mat[i][j1];

#if 1
//	const char *letters = "acgtn";
//	const char *letters2 = "ACGTN";
#else
//	const char *letters = "arndcqeghilkmfpstwyvbzx";
//	const char *letters2 = "ARNDCQEGHILKMFPSTWYVBZX";
#endif
	int back = back_mat[i][j1];
	int last = back;
	int count = 0, count2 = 0, count3 = 0;
	int match, begin1, begin2, end1, end2;
	int gbegin1, gbegin2, gend1, gend2;
	int64_t score, smin = best_score1, smax = best_score1 - 1;
	int posmin, posmax, pos = 0;
	int bl, dlen = 0, dcount = 0;
	posmin = posmax = 0;
	begin1 = begin2 = end1 = end2 = 0;
	gbegin1 = gbegin2 = gend1 = gend2 = 0;

#ifdef PRINT
#define PRINT
	if( options.verbose)printf( "%i %i\n", best_score, score_mat[i][j1] );
	if( options.verbose)printf( "%i %i %i\n", band_left, band_center, band_right );
	if( options.verbose)printf( "%i %i %i %i\n", i, j, j1, len2 );
#endif
#ifdef MAKEALIGN
#define MAKEALIGN
	char AA[ MAX_SEQ ], BB[ MAX_SEQ ];
	int NN = 0;
	int IA, IB;
	for(IA=len1;IA>i;IA--){
		AA[NN] = letters[ iseq1[IA-1] ];
		BB[NN++] = '-';
	}
	for(IB=len2;IB>j;IB--){
		AA[NN] = '-';
		BB[NN++] = letters[ iseq2[IB-1] ];
	}
#endif

	int indels = 0;
	int max_indels = 0;
	while( back != DP_BACK_NONE ){
		switch( back ){
		case DP_BACK_TOP  :
#ifdef PRINT
			if( options.verbose)printf( "%5i: %c %c %9i\n", pos, letters[ iseq1[i-1] ], '|', score_mat[i][j1] );
#endif
#ifdef MAKEALIGN
			AA[NN] = letters[ iseq1[i-1] ];
			BB[NN++] = '-';
#endif
			bl = (last != back) & (j != 1) & (j != len2);
			dlen += bl;
			dcount += bl;
			score = score_mat[i][j1];
			if( score < smin ){
				count2 = 0;
				smin = score;
				posmin = pos - 1;
				begin1 = i;
				begin2 = j;
			}
			i -= 1;
			j1 += 1;
			break;
		case DP_BACK_LEFT :
#ifdef PRINT
			if( options.verbose)printf( "%5i: %c %c %9i\n", pos, '|', letters[ iseq2[j-1] ], score_mat[i][j1] );
#endif
#ifdef MAKEALIGN
			AA[NN] = '-';
			BB[NN++] = letters[ iseq2[j-1] ];
#endif
			bl = (last != back) & (i != 1) & (i != len1);
			dlen += bl;
			dcount += bl;
			score = score_mat[i][j1];
			if( score < smin ){
				count2 = 0;
				smin = score;
				posmin = pos - 1;
				begin1 = i;
				begin2 = j;
			}
			j1 -= 1;
			j -= 1;
			break;
		case DP_BACK_LEFT_TOP :
#ifdef PRINT
			if( iseq1[i-1] == iseq2[j-1] ){
				if( options.verbose)printf( "%5i: %c %c %9i\n", pos, letters2[ iseq1[i-1] ], letters2[ iseq2[j-1] ], score_mat[i][j1] );
			}else{
				if( options.verbose)printf( "%5i: %c %c %9i\n", pos, letters[ iseq1[i-1] ], letters[ iseq2[j-1] ], score_mat[i][j1] );
			}
#endif
#ifdef MAKEALIGN
			if( iseq1[i-1] == iseq2[j-1] ){
				AA[NN] = letters2[ iseq1[i-1] ];
				BB[NN++] = letters2[ iseq2[j-1] ];
			}else{
				AA[NN] = letters[ iseq1[i-1] ];
				BB[NN++] = letters[ iseq2[j-1] ];
			}
#endif
			if( alninfo && options.global_identity ){
				if( i == 1 || j == 1 ){
					gbegin1 = i-1;
					gbegin2 = j-1;
				}else if( i == len1 || j == len2 ){
					gend1 = i-1;
					gend2 = j-1;
				}
			}
			score = score_mat[i][j1];
			i -= 1;
			j -= 1;
			match = iseq1[i] == iseq2[j];
			dlen += 1;
			dcount += ! match;
			if( score > smax ){
				count = 0;
				smax = score;
				posmax = pos;
				end1 = i;
				end2 = j;
			}
			count += match;
			count2 += match;
			count3 += match;
			if( score < smin ){
				int mm = match == 0;
				count2 = 0;
				smin = score;
				posmin = pos - mm;
				begin1 = i + mm;
				begin2 = j + mm;
			}
			break;
		default : 
			if( options.verbose)printf( "%i\n", back ); 
			break;
		}
		if( options.is454 ){
			if( back == DP_BACK_LEFT_TOP ){
				if( indels > max_indels ) max_indels = indels;
				indels = 0;
			}else{
				if( last == DP_BACK_LEFT_TOP ){
					indels = 1;
				}else if( indels ){
					indels += 1;
				}
			}
		}
		pos += 1;
		last = back;
		back = back_mat[i][j1];
	}
	if( options.is454 && max_indels > options.max_indel ) return FAILED_FUNC;
	iden_no = options.global_identity ? count3 : count - count2;
	alnln = posmin - posmax + 1;
	dist = dcount/(float)dlen;
	//dist = - 0.75 * log( 1.0 - dist * 4.0 / 3.0 );
	int umtail1 = len1 - 1 - end1;
	int umtail2 = len2 - 1 - end2;
	int umhead = begin1 < begin2 ? begin1 : begin2;
	int umtail = umtail1 < umtail2 ? umtail1 : umtail2;
	int umlen = umhead + umtail;
	if( umlen > options.unmatch_len ) return FAILED_FUNC;
	if( umlen > len1 * options.short_unmatch_per ) return FAILED_FUNC;
	if( umlen > len2 * options.long_unmatch_per ) return FAILED_FUNC;
	if( alninfo ){
		alninfo[0] = begin1;
		alninfo[1] = end1;
		alninfo[2] = begin2;
		alninfo[3] = end2;
		if( options.global_identity ){
			alninfo[0] = gbegin1;
			alninfo[1] = gend1;
			alninfo[2] = gbegin2;
			alninfo[3] = gend2;
		}
	}
#ifdef PRINT
	if( options.verbose)printf( "%6i %6i:  %4i %4i %4i %4i\n", alnln, iden_no, begin1, end1, begin2, end2 );
	if( options.verbose)printf( "%6i %6i:  %4i %4i\n", posmin, posmax, posmin - posmax, count - count2 );
	if( options.verbose)printf( "smin = %9i, smax = %9i\n", smin, smax );
	if( options.verbose)printf( "dlen = %5i, dcount = %5i, dist = %.3f\n", dlen, dcount, dcount/(float)dlen );
#endif
#ifdef MAKEALIGN
	float identity = iden_no / (float)( options.global_identity ? len1 : alnln);
	if( identity < options.cluster_thd ) return OK_FUNC;
	while(i--){
		AA[NN] = letters[ iseq1[i-1] ];
		BB[NN++] = '-';
	}
	while(j--){
		AA[NN] = '-';
		BB[NN++] = letters[ iseq2[j-1] ];
	}
	AA[NN] = '\0';
	BB[NN] = '\0';
	for(i=0; i<NN/2; i++){
		char aa = AA[i], bb = BB[i];
		AA[i] = AA[NN-i-1];
		BB[i] = BB[NN-i-1];
		AA[NN-i-1] = aa;
		BB[NN-i-1] = bb;
	}
	static int fcount = 0; 
	fcount += 1;
	FILE *fout = fopen( "alignments.txt", "a" );
	if( fout == NULL ){
		if( fcount <= 1 ) fprintf( stderr, "alignment files open failed\n" );
		return OK_FUNC;
	}
	fprintf( fout, "\n\n######################################################\n" );
	fprintf( fout, "# length X = %i\n", len2 );
	fprintf( fout, "# length Y = %i\n", len1 );
	fprintf( fout, "# best align X: %i-%i\n", begin2+1, end2+1 );
	fprintf( fout, "# best align Y: %i-%i\n", begin1+1, end1+1 );
	if( alninfo ){
		fprintf( fout, "# align X: %i-%i\n", alninfo[2]+1, alninfo[3]+1 );
		fprintf( fout, "# align Y: %i-%i\n", alninfo[0]+1, alninfo[1]+1 );
	}
	fprintf( fout, "# alignment length: %i\n", alnln );
	fprintf( fout, "# identity count: %i\n", iden_no );
	fprintf( fout, "# identity: %g\n", identity );
	fprintf( fout, "# distance: %g\n", dist );
	if( options.is454 ) fprintf( fout, "# max indel: %i\n", max_indels );
#if 0
	fprintf( fout, "%i %s\n", seq1->index, AA );
	fprintf( fout, "%i %s\n", seq2->index, BB );
#else
	bool printaa = true;
	IB = IA = 0;
	fprintf( fout, "\n\nX " );
	while( IA < NN ){
		if( printaa ){
			fprintf( fout, "%c", BB[IB] );
			IB += 1;
			if( IB % 75 ==0 or IB == NN ) printaa = false, fprintf( fout, "\nY " );
		}else{
			fprintf( fout, "%c", AA[IA] );
			IA += 1;
			if( IA % 75 ==0 ) printaa = true, fprintf( fout, "\n\nX " );
		}
	}
#endif
	fclose( fout );
#endif

	return OK_FUNC;
} // END int local_band_align


//--------- make Comp_AAN_idx -----//
void SequenceDB::make_comp_iseq(int len, char *iseq_comp, char *iseq)
{
	int i, c[5] = {3,2,1,0,4};
	for (i=0; i<len; i++) iseq_comp[i] = c[ (int)iseq[len-i-1] ];
} // make_comp_iseq

void SequenceDB::make_comp_short_word_index(int NAA, int *NAAN_array, Vector<int> &Comp_AAN_idx) {
	int i, j, k, icomp, k1;
	int c[4] = {3,2,1,0};
	unsigned char short_word[32]; //short_word[12] is enough

	int NAA1 = NAAN_array[1];
	int NAAN = NAAN_array[NAA];

	for (i=0; i<NAAN; i++) {
		// decompose i back to short_word
		for (k=i, j=0; j<NAA; j++) {
			short_word[j] = (unsigned char) (k % NAA1);
			k  = k / NAA1;
		}

		// calc_comp_aan_list
		icomp=0;
		for (k=0, k1=NAA-1; k<NAA; k++, k1--) icomp += c[short_word[k1]] * NAAN_array[k];

		Comp_AAN_idx[i] = icomp;
	}
} // make_comp_short_word_index


//----------- upper bound related operation -----------//
// when alignment coverage such as -aL is specified
// if a existing rep is too long, it won't be qulified 
int SequenceDB::upper_bound_length_rep(int len, double opt_s, int opt_S, double opt_aL, int opt_AL )
{
	int len_upper_bound = 99999999;
	double r1 = (opt_s > opt_aL) ? opt_s : opt_aL;
	int    a2 = (opt_S < opt_AL) ? opt_S : opt_AL;
	if (r1 > 0.0) len_upper_bound = (int) ( ((float) len)  / r1);
	if ((len+a2) < len_upper_bound)  len_upper_bound = len+a2;

	return len_upper_bound;
} // END upper_bound_length_rep
int SequenceDB::upper_bound_length_rep(int len, const Options & options )
{
	double opt_s = options.diff_cutoff;
	int    opt_S = options.diff_cutoff_aa;
	double opt_aL = options.long_coverage;
	int    opt_AL = options.long_control;
	return upper_bound_length_rep( len, opt_s, opt_S, opt_aL, opt_AL );
}

//------------- statistics related operation ---------//
void SequenceDB::cal_aax_cutoff(double &aa1_cutoff, double &aa2_cutoff, double &aan_cutoff,
		double cluster_thd, int tolerance, int naa_stat_start_percent,
		int naa_stat[5][61][4], int NAA)
{
	aa1_cutoff = cluster_thd;
	aa2_cutoff = 1 - (1-cluster_thd)*2;
	aan_cutoff = 1 - (1-cluster_thd)*NAA;
	if (tolerance==0) return; 

	int clstr_idx = (int) (cluster_thd * 100) - naa_stat_start_percent;
	if (clstr_idx <0) clstr_idx = 0;
	double d2  = ((double) (naa_stat[tolerance-1][clstr_idx][3]     )) / 100;
	double dn  = ((double) (naa_stat[tolerance-1][clstr_idx][5-NAA] )) / 100;
	aa2_cutoff = d2 > aa2_cutoff ? d2 : aa2_cutoff;
	aan_cutoff = dn > aan_cutoff ? dn : aan_cutoff;
	return;
} // END cal_aax_cutoff


void SequenceDB::update_aax_cutoff(double &aa1_cutoff, double &aa2_cutoff, double &aan_cutoff,
		int tolerance, int naa_stat_start_percent,
		int naa_stat[5][61][4], int NAA, double cluster_thd)
{
	if (cluster_thd > 1.0) cluster_thd = 1.00;

	double aa1_t, aa2_t, aan_t;
	cal_aax_cutoff(aa1_t, aa2_t, aan_t, cluster_thd, tolerance, naa_stat_start_percent,
			naa_stat, NAA);
	if (aa1_t > aa1_cutoff) aa1_cutoff = aa1_t;
	if (aa2_t > aa2_cutoff) aa2_cutoff = aa2_t;
	if (aan_t > aan_cutoff) aan_cutoff = aan_t;
	return;  
} // END update_aax_cutoff


//------------ calc_ann_list -----------//
int SequenceDB::calc_ann_list(int len, char *seqi, int NAA, int& aan_no, Vector<int> & aan_list, Vector<INTs> & aan_list_no, bool est) 
{
	int i, j, k, i0, i1, k1;

	// check_aan_list 
	aan_no = len - NAA + 1;
	for (j=0; j<aan_no; j++) {
		aan_list[j] = 0;
		for (k=0, k1=NAA-1; k<NAA; k++, k1--) aan_list[j] += seqi[j+k] * NAAN_array[k1];
	}
	if( est ){
		// for the short word containing 'N', mask it to '-1'
		for (j=0; j<len; j++){
			if ( seqi[j] == 4 ) {                      // here N is 4
				i0 = (j-NAA+1 > 0) ? j-NAA+1 : 0;
				i1 = j < aan_no ? j : aan_no - 1;
				for (i=i0; i<=i1; i++) aan_list[i]=-1;
			}
		}
	}

	sort(aan_list.begin(), aan_list.begin() + aan_no);
	for(j=0; j<aan_no; j++) aan_list_no[j]=1;
	for(j=aan_no-1; j; j--) {
		if (aan_list[j] == aan_list[j-1]) {
			aan_list_no[j-1] += aan_list_no[j];
			aan_list_no[j]=0;
		}
	}
	return OK_FUNC;
} // END calc_ann_list


//==================================== major SequenceDB class content ===============================//

//----------- class SequenceDB ---------//
void SequenceDB::Read( const char *file, const Options & options )
{
	Sequence one;
	Sequence dummy;
	Sequence des;
	FILE *swap = NULL;
	FILE *fin = fopen( file, "r" );
	char *buffer = NULL;
	char *res = NULL;
	size_t swap_size = 0;
	int option_l = options.min_length;
	if( fin == NULL ) CDHIT_bomb_error( "Failed to open the database file" );
	if( options.store_disk ) swap = CDHIT_OpenTempFile( CDHIT_temp_dir );
	Clear();
	dummy.swap = swap;
	buffer = new char[ MAX_LINE_SIZE+1 ];

	while (! feof( fin ) || one.size) { /* do not break when the last sequence is not handled */
		buffer[0] = '>';
		if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL && one.size == 0) break;
		if( buffer[0] == '+' ){
			int len = strlen( buffer );
			int len2 = len;
			while( len2 && buffer[len2-1] != '\n' ){
				if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
				len2 = strlen( buffer );
				len += len2;
			}
			one.des_length2 = len;
			dummy.des_length2 = len;
			fseek( fin, one.size, SEEK_CUR );
		}else if (buffer[0] == '>' || buffer[0] == '@' || (res==NULL && one.size)) {
			if ( one.size ) { // write previous record
				one.dat_length = dummy.dat_length = one.size;
				one.Format();
				one.index = dummy.index = sequences.size();
				if ( one.size > option_l ) {
					if ( swap ) {
						swap_size += one.size;
						// so that size of file < MAX_BIN_SWAP about 2GB
						if ( swap_size >= MAX_BIN_SWAP) {
							dummy.swap = swap = CDHIT_OpenTempFile( CDHIT_temp_dir );
							swap_size = one.size;
						}
						dummy.size = one.size;
						dummy.offset = ftell( swap );
						dummy.des_length = one.des_length;
						sequences.Append( new Sequence( dummy ) ); 
						one.ConvertBases();
						fwrite( one.data, 1, one.size, swap );
					}else{
						//printf( "==================\n" );
						sequences.Append( new Sequence( one ) ); 
						//printf( "------------------\n" );
						//if( sequences.size() > 10 ) break;
					}
					//if( sequences.size() >= 10000 ) break;
				}
			}
			one.size = 0;
			one.des_length2 = 0;

			int len = strlen( buffer );
			int len2 = len;
			des.size = 0;
			des += buffer;
			while( len2 && buffer[len2-1] != '\n' ){
				if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
				des += buffer;
				len2 = strlen( buffer );
				len += len2;
			}
			size_t offset = ftell( fin );
			one.des_begin = dummy.des_begin = offset - len;
			one.des_length = dummy.des_length = len;

			int i = 0;
			if( des.data[i] == '>' || des.data[i] == '@' || des.data[i] == '+' ) i += 1;
			if( des.data[i] == ' ' || des.data[i] == '\t' ) i += 1;
			if( options.des_len && options.des_len < des.size ) des.size = options.des_len;
			while( i < des.size && ! isspace( des.data[i] ) ) i += 1;
			des.data[i] = 0;
			one.identifier = dummy.identifier = des.data;
		} else {
			one += buffer;
		}
	}
#if 0
	int i, n = 0;
	for(i=0; i<sequences.size(); i++) n += sequences[i].bufsize + 4;
	cout<<n<<"\t"<<sequences.capacity() * sizeof(Sequence)<<endl;
	int i;
	scanf( "%i", & i );
#endif
	one.identifier = dummy.identifier = NULL;
	delete[] buffer;
	fclose( fin );
}



//-> note that the size of in_head and in_seq must equal !
void SequenceDB::Read_From_String( const Options & options,
	vector <string> & in_head, vector <string> & in_seq )
{
	//initial check
	if(in_head.size() != in_seq.size())CDHIT_bomb_error( "Size of in_head and in_seq should be equal" );
	//real process
	Sequence one;
	Sequence dummy;
	Sequence des;
	FILE *swap = NULL;
//	FILE *fin = fopen( file, "r" );
	char *buffer = NULL;
	size_t swap_size = 0;
	int option_l = options.min_length;
//	if( fin == NULL ) CDHIT_bomb_error( "Failed to open the database file" );
	if( options.store_disk ) swap = CDHIT_OpenTempFile( CDHIT_temp_dir );
	Clear();
	dummy.swap = swap;
	buffer = new char[ MAX_LINE_SIZE+1 ];

	//real process
	int wsi;
	int size=(int)in_seq.size();
	for(wsi=0;wsi<size;wsi++)
	{
		//get des
		{
			strcpy(buffer,in_head[wsi].c_str());
			des.size = 0;
			des += buffer;
			int len = strlen(buffer);
			one.des_length = dummy.des_length = len;
			int i = 0;
			if( des.data[i] == '>' || des.data[i] == '@' || des.data[i] == '+' ) i += 1;
			if( des.data[i] == ' ' || des.data[i] == '\t' ) i += 1;
			if( options.des_len && options.des_len < des.size ) des.size = options.des_len;
			while( i < des.size && ! isspace( des.data[i] ) ) i += 1;
			des.data[i] = 0;
			one.identifier = dummy.identifier = des.data;
		}//end des
		//get seq
		{
			strcpy(buffer,in_seq[wsi].c_str());
			one.size = 0;
			one += buffer;
			//appedn to data
			if ( one.size ) { // write previous record
				one.dat_length = dummy.dat_length = one.size;
				one.Format();
				one.index = dummy.index = sequences.size();
				if ( one.size > option_l ) {
					if ( swap ) {
						swap_size += one.size;
						// so that size of file < MAX_BIN_SWAP about 2GB
						if ( swap_size >= MAX_BIN_SWAP) {
							dummy.swap = swap = CDHIT_OpenTempFile( CDHIT_temp_dir );
							swap_size = one.size;
						}
						dummy.size = one.size;
						dummy.offset = ftell( swap );
						dummy.des_length = one.des_length;
						sequences.Append( new Sequence( dummy ) ); 
						one.ConvertBases();
						fwrite( one.data, 1, one.size, swap );
					}else{
						//printf( "==================\n" );
						sequences.Append( new Sequence( one ) ); 
						//printf( "------------------\n" );
						//if( sequences.size() > 10 ) break;
					}
					//if( sequences.size() >= 10000 ) break;
				}
			}
		}//end seq
	}
#if 0
	int i, n = 0;
	for(i=0; i<sequences.size(); i++) n += sequences[i].bufsize + 4;
	cout<<n<<"\t"<<sequences.capacity() * sizeof(Sequence)<<endl;
	int i;
	scanf( "%i", & i );
#endif
	one.identifier = dummy.identifier = NULL;
	delete[] buffer;
//	fclose( fin );
}


#if 0
void SequenceDB::Sort( int first, int last )
{
	int lower=first+1, upper=last;
	Sequence *pivot;
	Sequence *val;
	if( first >= last ) return;
	val = sequences[first];
	sequences[first] = sequences[ (first+last)/2 ];
	sequences[ (first+last)/2 ] = val;
	pivot = sequences[ first ];

	while( lower <= upper ){
		while( lower <= last && sequences[lower]->stats < pivot->stats ) lower ++;
		while( pivot->stats < sequences[upper]->stats ) upper --;
		if( lower < upper ){
			val = sequences[lower];
			sequences[lower] = sequences[upper];
			sequences[upper] = val;
			upper --;
		}
		lower ++;
	}
	val = sequences[first];
	sequences[first] = sequences[upper];
	sequences[upper] = val;
	if( first < upper-1 ) Sort( first, upper-1 );
	if( upper+1 < last ) Sort( upper+1, last );
}
#endif
void SequenceDB::SortDivide( Options & options, bool sort )
{
	int i, len;
	long long total_letter=0;
	long long total_desc=0;
	int max_len = 0, min_len = 99999999;
	int N = sequences.size();
	for (i=0; i<N; i++) {
		Sequence *seq = sequences[i];
		len = seq->size;
		total_letter += len;
		if (len > max_len) max_len = len;
		if (len < min_len) min_len = len;
		if (seq->swap == NULL) seq->ConvertBases();
		if( seq->identifier ) total_desc += strlen( seq->identifier );
	}
	options.total_letters = (size_t)total_letter;
	options.total_desc = (size_t)total_desc;
	options.max_length = max_len;
	options.max_entries = max_len * MAX_TABLE_SEQ;
	if (max_len >= 65536 && sizeof(INTs) <=2) 
		CDHIT_bomb_warning("Some seqs longer than 65536, you may define LONG_SEQ");
	if (max_len > MAX_SEQ ) 
		CDHIT_bomb_warning("Some seqs too long, you need to increase MAX_SEQ");
	if( options.verbose)cout << "longest and shortest : " << max_len << " and " << min_len << endl;
	if( options.verbose)cout << "Total letters: " << total_letter << endl;
	// END change all the NR_seq to iseq

	if( sort ){
		// **************************** Form NR_idx[], Sort them from Long to short
		int M = max_len - min_len + 1;
		Vector<int> count( M, 0 ); // count for each size = max_len - i
		Vector<int> accum( M, 0 ); // count for all size > max_len - i
		Vector<int> offset( M, 0 ); // offset from accum[i] when filling sorting
		Vector<Sequence*> sorting( N ); // TODO: use a smaller class if this consumes to much memory!

		for (i=0; i<N; i++) count[ max_len - sequences[i]->size ] ++;
		for (i=1; i<M; i++) accum[i] = accum[i-1] + count[i-1];
		for (i=0; i<N; i++){
			int len = max_len - sequences[i]->size;
			int id = accum[len] + offset[len];
			//sequences[i].index = id;
			sorting[id] = sequences[i];
			offset[len] ++;
		}
		options.max_entries = 0;
		for (i=0; i<N; i++){
			sequences[i] = sorting[i];
			if( i < MAX_TABLE_SEQ ) options.max_entries += sequences[i]->size;
		}
#if 0
		if( options.isEST ){
			int start = 0;
			for (i=0; i<M; i++){
				Sort( start, accum[i] );
				start = accum[i];
			}
		}
#endif
		if( options.verbose)cout << "Sequences have been sorted" << endl;
		// END sort them from long to short
	}
}// END sort_seqs_divide_segs

void SequenceDB::DivideSave( const char *db, const char *newdb, int n, const Options & options )
{
	if( n == 0 || sequences.size() ==0 ) return;

	size_t max_seg = options.total_letters / n + sequences[0]->size;
	if( max_seg >= MAX_BIN_SWAP ) max_seg = (size_t) MAX_BIN_SWAP;

	FILE *fin = fopen( db, "r" );
	char *buf = new char[MAX_LINE_SIZE+1];
	char outfile[512];
	size_t seg_size = 0;
	int i, j, count, rest, seg = 0;
	sprintf( outfile, "%s-%i", newdb, 0 );
	FILE *fout = fopen( outfile, "w+" );
	n = sequences.size();
	for (i=0; i<n; i++){
		Sequence *seq = sequences[i];
		int qs = seq->des_length2 ? seq->des_length2 + seq->dat_length : 0;
		fseek( fin, seq->des_begin, SEEK_SET );

		seg_size += seq->size;
		if( seg_size >= max_seg ){
			seg += 1;
			sprintf( outfile, "%s-%i", newdb, seg );
			fclose( fout );
			fout = fopen( outfile, "w+" );
			seg_size = seq->size;
		}

		count = (seq->des_length + seq->dat_length + qs) / MAX_LINE_SIZE;
		rest = (seq->des_length + seq->dat_length + qs) % MAX_LINE_SIZE;
		//printf( "count = %6i,  rest = %6i\n", count, rest );
		for (j=0; j<count; j++){
			if( fread( buf, 1, MAX_LINE_SIZE, fin ) ==0 ) CDHIT_bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, MAX_LINE_SIZE, fout );
		}
		if( rest ){
			if( fread( buf, 1, rest, fin ) ==0 ) CDHIT_bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, rest, fout );
		}
	}
	fclose( fin );
	fclose( fout );
	delete []buf;
}
void SequenceDB::WriteClusters( const char *db, const char *newdb, const Options & options )
{
	FILE *fin = fopen( db, "r" );
	FILE *fout = fopen( newdb, "w+" );
	int i, j, n = rep_seqs.size();
	int count, rest;
	char *buf = new char[MAX_LINE_SIZE+1];
	vector<uint64_t> sorting( n );
	if( fin == NULL || fout == NULL ) CDHIT_bomb_error( "file opening failed" );
	for (i=0; i<n; i++) sorting[i] = ((uint64_t)sequences[ rep_seqs[i] ]->index << 32) | rep_seqs[i];
	sort( sorting.begin(), sorting.end() );
	for (i=0; i<n; i++){
		Sequence *seq = sequences[ sorting[i] & 0xffffffff ];
		int qs = seq->des_length2 ? seq->des_length2 + seq->dat_length : 0;
		fseek( fin, seq->des_begin, SEEK_SET );

		count = (seq->des_length + seq->dat_length + qs) / MAX_LINE_SIZE;
		rest = (seq->des_length + seq->dat_length + qs) % MAX_LINE_SIZE;
		//printf( "count = %6i,  rest = %6i\n", count, rest );
		for (j=0; j<count; j++){
			if( fread( buf, 1, MAX_LINE_SIZE, fin ) ==0 ) CDHIT_bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, MAX_LINE_SIZE, fout );
		}
		if( rest ){
			if( fread( buf, 1, rest, fin ) ==0 ) CDHIT_bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, rest, fout );
		}
	}
	fclose( fin );
	fclose( fout );
	delete []buf;
}

//----- write cluster center and member to string -----//__140310__//
void SequenceDB::WriteExtra1D_To_String( const Options & options,
	vector <vector <string> > &member, vector < int >&center )
{
	//init
	center.clear();
	member.clear();
	//process
	int i, k, N = sequences.size();
	vector<long long> sorting( N );
	for (i=0; i<N; i++) sorting[i] = ((long long)sequences[i]->index << 32) | i;
	sort( sorting.begin(), sorting.end() );
	int M = rep_seqs.size();
	Vector<Vector<int> > clusters( M );
	for (i=0; i<N; i++){
		int k = sorting[i] & 0xffffffff;
		int id = sequences[k]->cluster_id;
		clusters[id].Append( k );
	}
	//output
	for (i=0; i<M; i++)
	{
		//init
		vector <string> memb;
		memb.clear();
		int cent=0;
		//proc
		for (k=0; k<(int)clusters[i].size(); k++)
		{
			string str=sequences[ clusters[i][k] ]->identifier;
			memb.push_back(str);
			if(!sequences[ clusters[i][k] ]->identity) //cent
			{
				if(cent==0)center.push_back(k);
				cent=1;
			}
		}
		//final
		member.push_back(memb);
	}
}

void SequenceDB::WriteExtra1D( const Options & options )
{
	string db_clstr = options.output + ".clstr";
	string db_clstr_bak = options.output + ".bak.clstr";
	int i, k, N = sequences.size();
	vector<long long> sorting( N );
	for (i=0; i<N; i++) sorting[i] = ((long long)sequences[i]->index << 32) | i;
	sort( sorting.begin(), sorting.end() );

	FILE *fin = fopen( options.input.c_str(), "r" );
	FILE *fout = fopen( db_clstr_bak.c_str(), "w+" );
	char *buf = new char[ MAX_DES + 1 ];
	for (i=0; i<N; i++) {
		Sequence *seq = sequences[ sorting[i] & 0xffffffff ];
		seq->PrintInfo( seq->cluster_id, fin, fout, options, buf );
	}
	fclose( fout );

	if( options.verbose)cout << "writing clustering information" << endl;
	int M = rep_seqs.size();
	Vector<Vector<int> > clusters( M );
	for (i=0; i<N; i++){
		int k = sorting[i] & 0xffffffff;
		int id = sequences[k]->cluster_id;
		clusters[id].Append( k );
	}

	fout = fopen( db_clstr.c_str(), "w+" );
	for (i=0; i<M; i++) {
		fprintf( fout, ">Cluster %i\n", i );
		for (k=0; k<(int)clusters[i].size(); k++)
			sequences[ clusters[i][k] ]->PrintInfo( k, fin, fout, options, buf );
	}
	delete []buf;
	fclose( fin );
}
void SequenceDB::WriteExtra2D( SequenceDB & other, const Options & options )
{
	string db_clstr = options.output + ".clstr";
	string db_clstr_bak = options.output + ".bak.clstr";
	int i, k, N = other.sequences.size();
	int N2 = sequences.size();
	vector<long long> sorting( N );
	for (i=0; i<N; i++) sorting[i] = ((long long)other.sequences[i]->index << 32) | i;
	sort( sorting.begin(), sorting.end() );

	FILE *fin = fopen( options.input.c_str(), "r" );
	FILE *fin2 = fopen( options.input2.c_str(), "r" );
	FILE *fout = fopen( db_clstr_bak.c_str(), "w+" );
	char *buf = new char[ MAX_DES + 1 ];
	for (i=0; i<N; i++) {
		Sequence *seq = other.sequences[ sorting[i] & 0xffffffff ];
		seq->PrintInfo( seq->cluster_id, fin, fout, options, buf );
	}
	for (i=0; i<N2; i++) {
		Sequence *seq = sequences[i];
		if( seq->state & IS_REDUNDANT ) seq->PrintInfo( seq->cluster_id, fin2, fout, options, buf );
	}
	fclose( fout );

	if( options.verbose)cout << "writing clustering information" << endl;
	Vector<Vector<int> > clusters( N );
	for (i=0; i<N2; i++){
		int id = sequences[i]->cluster_id;
		if( sequences[i]->state & IS_REDUNDANT ) clusters[id].Append( i );
	}

	fout = fopen( db_clstr.c_str(), "w+" );
	for (i=0; i<N; i++) {
		Sequence *seq = other.sequences[ i ];
		fprintf( fout, ">Cluster %i\n", i );
		seq->PrintInfo( 0, fin, fout, options, buf );
		for (k=0; k<(int)clusters[i].size(); k++)
			sequences[ clusters[i][k] ]->PrintInfo( k+1, fin2, fout, options, buf );
	}
	delete []buf;
	fclose( fin );
}
void SequenceDB::ClusterOne( Sequence *seq, int id, WordTable & table,
		WorkingParam & param, WorkingBuffer & buffer, const Options & options )
{
	if (seq->state & IS_REDUNDANT) return;
	int frag_size = options.frag_size;
	int NAA = options.NAA;
	int len = seq->size;
	int len_bound = upper_bound_length_rep(len, options);
	param.len_upper_bound = len_bound;
	int flag = CheckOne( seq, table, param, buffer, options );

	if( flag == 0 ){
		if ((seq->identity>0) && (options.cluster_best)) {
			// because of the -g option, this seq is similar to seqs in old SEGs
			seq->state |= IS_REDUNDANT ;
			seq->Clear();
		} else {                  // else add to NR90 db
			int aan_no = len - NAA + 1;
			int size = rep_seqs.size();
			rep_seqs.Append( id );
			seq->cluster_id = size;
			seq->identity = 0;
			seq->state |= IS_REP;
			if (frag_size){ /* not used for EST */
				int frg1 = (len - NAA ) / frag_size + 1;
				table.AddWordCountsFrag( aan_no, buffer.word_encodes_backup, 
						buffer.word_encodes_no, frg1, frag_size );
			}else{
				table.AddWordCounts(aan_no, buffer.word_encodes, buffer.word_encodes_no, table.sequences.size(), options.isEST);
			}
			table.sequences.Append( seq );
			if( frag_size ){
				while( (int)table.sequences.size() < table.frag_count )
					table.sequences.Append( seq );
			}
		}
	}
	if ( (id+1) % 1000 == 0 ) {
		int size = rep_seqs.size();
		if( options.verbose)
		{
			printf( "." );
			fflush( stdout );
		}
		if ( (id+1) % 10000 == 0 ) if( options.verbose) printf( "\r..........%9i  finished  %9i  clusters\n", id+1, size );
	}
}
size_t SequenceDB::MinimalMemory( int frag_no, int bsize, int T, const Options & options )
{
	int N = sequences.size();
	int F = frag_no < MAX_TABLE_SEQ ? frag_no : MAX_TABLE_SEQ;
	size_t mem_need = 0;
	size_t mem, mega = 1000000;
	int table = T > 1 ? 2 : 1;

	if( options.verbose)printf( "\nApproximated minimal memory consumption:\n" );
	mem = N*sizeof(Sequence) + options.total_desc + N;
	if( options.store_disk == false ) mem += options.total_letters + N;
	if( options.verbose)printf( "%-16s: %iM\n", "Sequence", (int)(mem/mega) );
	mem_need += mem;

	mem = bsize;
	if( options.verbose)printf( "%-16s: %i X %iM = %iM\n", "Buffer", T, (int)(mem/mega), (int)(T*mem/mega) );
	mem_need += T*mem;

	mem = F*(sizeof(Sequence*) + sizeof(IndexCount)) + NAAN*sizeof(NVector<IndexCount>);
	if( options.verbose)printf( "%-16s: %i X %iM = %iM\n", "Table", table, (int)(mem/mega), (int)(table*mem/mega) );
	mem_need += table*mem;

	mem = sequences.capacity()*sizeof(Sequence*) + N*sizeof(int);
	mem += Comp_AAN_idx.size()*sizeof(int);
	if( options.verbose)printf( "%-16s: %iM\n", "Miscellaneous", (int)(mem/mega) );
	mem_need += mem;

	if( options.verbose)printf( "%-16s: %iM\n\n", "Total", (int)(mem_need/mega) );

	if(options.max_memory && options.max_memory < mem_need + 50*table ){
		char msg[200];
		sprintf( msg, "not enough memory, please set -M option greater than %i\n", 
				(int)(50*table + mem_need/mega) );
		CDHIT_bomb_error(msg);
	}
	return mem_need;
}
void SequenceDB::DoClustering( int T, const Options & options )
{
	int i, k;
	int NAA = options.NAA;
	double aa1_cutoff = options.cluster_thd;
	double aas_cutoff = 1 - (1-options.cluster_thd)*4;
	double aan_cutoff = 1 - (1-options.cluster_thd)*options.NAA;
	int seq_no = sequences.size();
	int frag_no = seq_no;
	int frag_size = options.frag_size;
	valarray<size_t>  letters(T);

	//printf( "%li\n", options.mem_limit );

	if (frag_size){ 
		frag_no = 0;
		for (i=0; i<seq_no; i++) frag_no += (sequences[i]->size - NAA) / frag_size + 1;
	}

	if( ! options.isEST )
		cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd,
				options.tolerance, naa_stat_start_percent, naa_stat, NAA);

	Vector<WorkingParam> params(T);
	Vector<WorkingBuffer> buffers(T);
	for(i=0; i<T; i++){
		params[i].Set( aa1_cutoff, aas_cutoff, aan_cutoff );
		buffers[i].Set( frag_no, options );
	}

	// word_table as self comparing table and table buffer:
	WordTable word_table( options.NAA, NAAN );

	WordTable last_table( options.NAA, NAAN );

	int N = sequences.size();
	size_t mem_need = MinimalMemory( frag_no, buffers[0].total_bytes, T, options );
	size_t mem, mega = 1000000;

	size_t tabsize = 0;
	size_t mem_limit = (options.max_memory - mem_need) / sizeof(IndexCount);

	if( options.verbose)printf( "Table limit with the given memory limit:\n" );
	if( options.verbose)printf( "Max number of representatives: %i\n", MAX_TABLE_SEQ );
	if( options.max_memory ){
		if( options.verbose)printf( "Max number of word counting entries: %li\n", mem_limit );
	}else{
		mem_limit = options.max_entries;
		if( (int)mem_limit > T * MAX_TABLE_SIZE ) mem_limit = T * MAX_TABLE_SIZE;
	}
	if( options.verbose)printf( "\n" );

	size_t mem_limit2 = mem_limit / 64;
	if( mem_limit2 > 1E6 ) mem_limit2 = (size_t)1E6;

	omp_set_num_threads(T);
	for(i=0; i<N; ){
		int start = i;
		int m = i;
		size_t sum = 0;
		size_t lim = mem_limit / T;
		float redundancy = (rep_seqs.size() + 1.0) / (i + 1.0);
		if( i ==0 ) lim /= 8; // first SCB with small size
		if( lim < mem_limit2 ) lim = (lim + mem_limit2) / 2; // SCB size has lower limit
		if(lim==0)lim=1;
		while( m < N && sum < lim ){


			Sequence *seq = sequences[m];
			if( ! (seq->state & IS_REDUNDANT) ){
				if ( options.store_disk ) seq->SwapIn();
				sum += (size_t)(seq->size * redundancy);
			}
			if( i ==0 && (m > 1E6 || sum > 1E8) ) break;
			m ++;
		}
		if( m >= N ){
			m = N;
			if( m > i + 1E3 ) m = i + (N - i) / T;
		}
		for(k=i; k<m; k++){
			Sequence *seq = sequences[k];
			if( ! (seq->state & IS_REDUNDANT) ){
				mem_limit -= (size_t)(seq->size * redundancy);
			}
		}
		//printf( "m = %i  %i,  %i\n", i, m, m-i );
		if( options.verbose)printf( "\r# comparing sequences from  %9i  to  %9i\n", i, m );
		if( last_table.size ){
			int print = (m-i)/20 + 1;
			#pragma omp parallel for schedule( dynamic, 1 )
			for(int j=i; j<m; j++){
				Sequence *seq = sequences[j];
				if (seq->state & IS_REDUNDANT) continue;
				int tid = omp_get_thread_num();
				CheckOne( seq, last_table, params[tid], buffers[tid], options );
				if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
				if( j%print==0 ){
					if( options.verbose)
					{
						printf( "." ); 
						fflush( stdout );
					}
				}
			}
			int may_stop = 0;
			int self_stop = 0;
			int JN = N;
			float p0 = 0;
			int min = last_table.sequences[ last_table.sequences.size()-1 ]->size;
			int m0 = m;
			#pragma omp parallel for schedule( dynamic, 1 )
			for(int j=m-1; j<JN; j++){
				if( j+1 == JN ){
					//printf( "stoping\n" );
					//#pragma omp atomic
					may_stop = 1;
				}
				int tid = omp_get_thread_num();
				if( j == (m0-1) ){ // use m0 to avoid other iterations satisfying the condition:
					int ks0 = i;
					for(int ks=i; ks<m; ks++){
						Sequence *seq = sequences[ks];
						i = ks + 1;
						if (seq->state & IS_REDUNDANT) continue;
						ClusterOne( seq, ks, word_table, params[tid], buffers[tid], options );
						if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
						if( may_stop && ks > (ks0+1000) ) break;
						if( word_table.size >= mem_limit ) break;
						int tmax = MAX_TABLE_SEQ - (frag_size ? seq->size / frag_size + 1 : 0);
						if( (int)word_table.sequences.size() >= tmax || (int)word_table.frag_count >= tmax ) break;
					}
					self_stop = 1;
				}else{
					Sequence *seq = sequences[j];
					if (seq->state & IS_REDUNDANT) continue;
					if ( options.store_disk ){
						#pragma omp critical
						seq->SwapIn();
					}
					CheckOne( seq, last_table, params[tid], buffers[tid], options );
					if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
					if( min > params[tid].len_upper_bound ){
						may_stop = 1;
						#pragma omp critical
						JN = j;
						continue;
					}
					if( self_stop && tid ==1 ){
						float p = (100.0*j)/N;
						if( p > p0+1E-1 ){ // print only if the percentage changed
							if( options.verbose)
							{
								printf( "\r%4.1f%%", p );
								fflush( stdout );
							}
							p0 = p;
						}
					}

				}
			}
		}
		if( i == start || m == N ){//{|| word_table.size < mem_limit2 ){
			//printf( "comparing the first or last or very small group ...\n" ); fflush( stdout );
			for(k=i; k<m; ){
				int kk, mm = k, sum = 0;
				while( mm < m && sum < 1E5 ){
					if( ! (sequences[mm]->state & IS_REDUNDANT) ) sum += sequences[mm]->size;
					mm += 1;
				}
				if( mm < k + 1000 ) mm = k + 1000;
				if( mm > m ) mm = m;
				#pragma omp parallel for schedule( dynamic, 1 )
				for(kk=k; kk<mm; kk++){
					Sequence *seq = sequences[kk];
					if (seq->state & IS_REDUNDANT) continue;
					int tid = omp_get_thread_num();
					CheckOne( seq, word_table, params[tid], buffers[tid], options );
					if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
				}
				bool bk = false;
				for(int ks=k; ks<mm; ks++){
					Sequence *seq = sequences[ks];
					i = k = ks + 1;
					if (seq->state & IS_REDUNDANT) continue;
					ClusterOne( seq, ks, word_table, params[0], buffers[0], options );
					bk = true;
					if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
					if( word_table.size >= mem_limit ) break;
					int tmax = MAX_TABLE_SEQ - (frag_size ? seq->size / frag_size + 1 : 0);
					if( (int)word_table.sequences.size() >= tmax || (int)word_table.frag_count >= tmax ) break;
					bk = false;
				}
				if( bk ) break;
			}
		}else if( i < m ){
			if( options.verbose)printf( "\r---------- %6i remaining sequences to the next cycle\n", m-i );
		}
		if( options.verbose)printf( "---------- new table with %8i representatives\n", (int)word_table.sequences.size() );
		if( (last_table.size + word_table.size) > tabsize )
			tabsize = last_table.size + word_table.size;
		last_table.Clear();
		last_table.sequences.swap( word_table.sequences );
		last_table.indexCounts.swap( word_table.indexCounts );
		last_table.size = word_table.size;
		word_table.size = 0;
	}
	if( options.verbose)printf( "\n%9li  finished  %9li  clusters\n", (long)sequences.size(), (long)rep_seqs.size() );
	mem = (mem_need + tabsize*sizeof(IndexCount))/mega;
	if( options.verbose)printf( "\nApprixmated maximum memory consumption: %iM\n", (int)mem );
	last_table.Clear();
	word_table.Clear();
}

int SequenceDB::CheckOne( Sequence *seq, WordTable & table, WorkingParam & param, WorkingBuffer & buf, const Options & options )
{
	int len = seq->size;
	param.len_upper_bound = upper_bound_length_rep(len, options);
	if( options.isEST ) return CheckOneEST( seq, table, param, buf, options );
	return CheckOneAA( seq, table, param, buf, options );
}
int SequenceDB::CheckOneAA( Sequence *seq, WordTable & table, WorkingParam & param, WorkingBuffer & buf, const Options & options )
{
	NVector<IndexCount> & lookCounts = buf.lookCounts;
	NVector<uint32_t> & indexMapping = buf.indexMapping;
	Vector<INTs> & word_encodes_no = buf.word_encodes_no;
	Vector<int>  & word_encodes = buf.word_encodes;
	double aa1_cutoff = param.aa1_cutoff;
	double aa2_cutoff = param.aas_cutoff;
	double aan_cutoff = param.aan_cutoff;

	char *seqi = seq->data;
	int k, j1, len = seq->size;
	int flag = 0;
	int frag_size = options.frag_size;
	int & aln_cover_flag = param.aln_cover_flag;
	int & required_aa1 = param.required_aa1;
	int & required_aa2 = param.required_aas;
	int & required_aan = param.required_aan;
	int & min_aln_lenS = param.min_aln_lenS;
	int & min_aln_lenL = param.min_aln_lenL;

	int NAA = options.NAA;

	param.ControlShortCoverage( len, options );
	param.ComputeRequiredBases( options.NAA, 2, options );

	buf.EncodeWords( seq, options.NAA, false );

	// if minimal alignment length > len, return
	// I can not return earlier, because I need to calc the word_encodes etc
	if (options.min_control>len) return 0; // return flag=0

	// lookup_aan
	int aan_no = len - options.NAA + 1;
	table.CountWords(aan_no, word_encodes, word_encodes_no, lookCounts, indexMapping, false, required_aan);

	// contained_in_old_lib()
	int len_upper_bound = param.len_upper_bound;
	int len_lower_bound = param.len_lower_bound;
	int band_left, band_right, best_score, band_width1, best_sum, len2, alnln, len_eff1;
	int tiden_no, band_center;
	float tiden_pc, distance=0;
	int talign_info[5];
	int sum;
	char *seqj;
	int frg2 = frag_size ? (len - NAA + options.band_width ) / frag_size + 1 + 1 : 0;
	int lens;
	int has_aa2 = 0;

	IndexCount *ic = lookCounts.items;
	ic = lookCounts.items;
	for(; ic->count; ic++){
		if( ! frag_size ){
			indexMapping[ ic->index ] = 0;
			if ( ic->count < required_aan ) continue;
		}

		Sequence *rep = table.sequences[ ic->index ];
		len2 = rep->size;
		if (len2 > len_upper_bound ) continue;
		if (options.has2D && len2 < len_lower_bound ) continue;
		if( frag_size ){
			uint32_t count = ic->count;
			uint32_t *ims = & indexMapping[ ic->index ];
			k = (len2 - NAA) / frag_size + 1;
			sum = 0;
			for (j1=0; j1<frg2; j1++){
				uint32_t im = ims[j1];
				if( im ) sum += lookCounts[im-1].count;
			}
			count = sum;
			for (j1=frg2; j1<k; j1++) {
				uint32_t im1 = ims[j1];
				uint32_t im2 = ims[j1-frg2];
				if( im1 ) sum += lookCounts[im1-1].count;
				if( im2 ) sum -= lookCounts[im2-1].count;
				if ( sum > (int)count) count = sum;
			}
			if ( (int)count < required_aan ) continue;
		}

		param.ControlLongCoverage( len2, options );

		if ( has_aa2 == 0 )  { // calculate AAP array
			buf.ComputeAAP( seqi, seq->size );
			has_aa2 = 1;
		}
		seqj = rep->data; //NR_seq[NR90_idx[j]];

		band_width1 = (options.band_width < len+len2-2 ) ? options.band_width : len+len2-2;
		diag_test_aapn(NAA1, seqj, len, len2, buf, options, best_sum,
				band_width1, band_left, band_center, band_right, required_aa1);
		if ( best_sum < required_aa2 ) continue;

		int rc = FAILED_FUNC;
		if (options.print || aln_cover_flag) //return overlap region
			rc = local_band_align(seqi, seqj, len, len2, mat,
					best_score, tiden_no, alnln, distance, talign_info+1,
					band_left, band_center, band_right, buf, options);
		else
			rc = local_band_align(seqi, seqj, len, len2, mat,
					best_score, tiden_no, alnln, distance, NULL, 
					band_left, band_center, band_right, buf, options);
		if ( rc == FAILED_FUNC ) continue;
		if ( tiden_no < required_aa1 ) continue;
		lens = len;
		if( options.has2D && len > len2 ) lens = len2;
		len_eff1 = (options.global_identity == 0) ? alnln : lens;
		tiden_pc = tiden_no / (float) len_eff1;
		if( options.useDistance ){
			if (distance > options.distance_thd ) continue;
			if (distance >= seq->distance) continue; // existing distance
		}else{
			if (tiden_pc < options.cluster_thd) continue;
			if (tiden_pc <= seq->identity) continue; // existing iden_no
		}
		if (aln_cover_flag) {
			if ( talign_info[4]-talign_info[3]+1 < min_aln_lenL) continue;
			if ( talign_info[2]-talign_info[1]+1 < min_aln_lenS) continue;
		}
		if( options.has2D ) seq->state |= IS_REDUNDANT ;
		flag = 1; seq->identity = tiden_pc; seq->cluster_id = rep->cluster_id;
		seq->distance = distance;
		seq->coverage[0] = talign_info[1] +1;
		seq->coverage[1] = talign_info[2] +1;
		seq->coverage[2] = talign_info[3] +1;
		seq->coverage[3] = talign_info[4] +1;
		if (! options.cluster_best) break;
		update_aax_cutoff(aa1_cutoff, aa2_cutoff, aan_cutoff,
				options.tolerance, naa_stat_start_percent, naa_stat, NAA, tiden_pc);
		param.ComputeRequiredBases( options.NAA, 2, options );
	}
	if( frag_size ) ic = lookCounts.items;
	while( ic->count ){
		indexMapping[ ic->index ] = 0;
		ic += 1;
	}
	lookCounts.size = 0;
	if (flag == 1) { // if similar to old one delete it
		if (! options.cluster_best) {
			seq->Clear();
			seq->state |= IS_REDUNDANT ;
		}
	}
	return flag;

}
int SequenceDB::CheckOneEST( Sequence *seq, WordTable & table, WorkingParam & param, WorkingBuffer & buf, const Options & options )
{
	NVector<IndexCount> & lookCounts = buf.lookCounts;
	NVector<uint32_t> & indexMapping = buf.indexMapping;
	Vector<INTs> & word_encodes_no = buf.word_encodes_no;
	Vector<int>  & word_encodes = buf.word_encodes;
	Vector<int> & aan_list_comp = buf.aan_list_comp;
	char *seqi_comp = buf.seqi_comp;

	int & aln_cover_flag = param.aln_cover_flag;
	int & required_aa1 = param.required_aa1;
	int & required_aas = param.required_aas;
	int & required_aan = param.required_aan;
	int & min_aln_lenS = param.min_aln_lenS;
	int & min_aln_lenL = param.min_aln_lenL;

	char *seqi = seq->data;
	int j, len = seq->size;
	int flag = 0;

	param.ControlShortCoverage( len, options );
	param.ComputeRequiredBases( options.NAA, 4, options );
	int skip = buf.EncodeWords( seq, options.NAA, true );
	required_aan -= skip;
	if( required_aan <= 0 ) required_aan = 1;

	// if minimal alignment length > len, return
	// I can not return earlier, because I need to calc the word_encodes etc
	if (options.min_control>len) return 0; // return flag=0

	int aan_no = len - options.NAA + 1;

	// contained_in_old_lib()
	int len_upper_bound = param.len_upper_bound;
	int len_lower_bound = param.len_lower_bound;
	int band_left, band_right, best_score, band_width1, best_sum, len2, alnln, len_eff1;
	int tiden_no, band_center;
	float tiden_pc, distance=0;
	int talign_info[5];
	int j0, comp;
	char *seqj;

	for(comp=0; comp<2; comp++){
		if( comp ){
			for (j0=0; j0<aan_no; j0++) {
				j = word_encodes[j0];
				if ( j<0 ) aan_list_comp[j0] = j;
				else       aan_list_comp[j0] = Comp_AAN_idx[j];
			}
			make_comp_iseq(len, seqi_comp, seqi);
			seqi = seqi_comp;
		}
		int has_aas = 0;

		if( comp ){
			table.CountWords(aan_no, aan_list_comp, word_encodes_no, lookCounts, indexMapping, true, required_aan );
		}else{
			table.CountWords(aan_no, word_encodes, word_encodes_no, lookCounts, indexMapping, true, required_aan ); 
		}

		IndexCount *ic = lookCounts.items;
		ic = lookCounts.items;
		for(; ic->count; ic++){
			indexMapping[ ic->index ] = 0;
			if ( ic->count < required_aan ) continue;
			Sequence *rep = table.sequences[ic->index];

			len2 = rep->size;
			if (len2 > len_upper_bound ) continue;
			if (options.has2D && len2 < len_lower_bound ) continue;

			seqj = rep->data;

			param.ControlLongCoverage( len2, options );

			if ( has_aas == 0 )  { // calculate AAP array
				buf.ComputeAAP2( seqi, seq->size );
				has_aas = 1;
			}

			band_width1 = (options.band_width < len+len2-2 ) ? options.band_width : len+len2-2;
			diag_test_aapn_est(NAA1, seqj, len, len2, buf, options, best_sum,
					band_width1, band_left, band_center, band_right, required_aa1);
			//printf( "%i %i\n", best_sum, required_aas );
			if ( best_sum < required_aas ) continue;
			//if( comp and flag and (not options.cluster_best) and j > rep->cluster_id ) goto Break;

			int rc = FAILED_FUNC;
			if (options.print || aln_cover_flag){ //return overlap region
				rc = local_band_align(seqi, seqj, len, len2, mat,
						best_score, tiden_no, alnln, distance, talign_info+1,
						band_left, band_center, band_right, buf, options);
				if( comp ){
					talign_info[1] = len - talign_info[1] - 1;
					talign_info[2] = len - talign_info[2] - 1;
				}
			}else{
				//printf( "%5i %5i %5i %5i\n", band_width1, band_right-band_left, band_left, band_right );
				rc = local_band_align(seqi, seqj, len, len2, mat,
						best_score, tiden_no, alnln, distance, talign_info+1,
						band_left, band_center, band_right, buf, options);
			}
			if ( rc == FAILED_FUNC ) continue;
			//printf( "%i  %i  %i\n", best_score, tiden_no, required_aa1 );
			if ( tiden_no < required_aa1 ) continue;
			if ( options.is454 ){
				if (talign_info[3] != talign_info[1]) continue; // same start
				if (talign_info[1] > 1) continue; // one mismatch allowed at beginning
				if ((len-talign_info[2]) > 2) continue; // one mismatch allowed at end
			}

			len_eff1 = (options.global_identity == 0) ? alnln : len;
			tiden_pc = tiden_no / (float)len_eff1;
			//printf( "%i %i\n", tiden_no, options.cluster_thd100 );
			if( options.useDistance ){
				if (distance > options.distance_thd ) continue;
				if (options.cluster_best && distance >= seq->distance) continue; // existing distance
			}else{
				if (tiden_pc < options.cluster_thd) continue;
				if (options.cluster_best && tiden_pc < seq->identity) continue; // existing iden_no
			}
			if (aln_cover_flag) {
				if ( talign_info[4]-talign_info[3]+1 < min_aln_lenL) continue;
				if( comp ){
					if ( talign_info[1]-talign_info[2]+1 < min_aln_lenS) continue;
				}else{
					if ( talign_info[2]-talign_info[1]+1 < min_aln_lenS) continue;
				}
			}
			if( options.cluster_best && fabs(tiden_pc - seq->identity) < 1E-9 && rep->cluster_id >= seq->cluster_id ) continue;
			if( (! options.cluster_best) && flag !=0 && rep->cluster_id >= seq->cluster_id ) continue;
			flag = comp ? -1 : 1;
			seq->identity = tiden_pc;
			seq->distance = distance;
			seq->cluster_id = rep->cluster_id;
			seq->coverage[0] = talign_info[1] +1;
			seq->coverage[1] = talign_info[2] +1;
			seq->coverage[2] = talign_info[3] +1;
			seq->coverage[3] = talign_info[4] +1;
			if (! options.cluster_best) break;
		}
		while( ic->count ){
			indexMapping[ ic->index ] = 0;
			ic += 1;
		}
		lookCounts.size = 0;
		if (! options.option_r ) break;
	}
	if ((flag == 1) || (flag == -1)) { // if similar to old one delete it
		if (! options.cluster_best) {
			seq->Clear();
			seq->state |= IS_REDUNDANT ;
		}
		if( flag == -1 )
			seq->state |= IS_MINUS_STRAND; 
		else
			seq->state &= ~IS_MINUS_STRAND; 
	}
	return flag;
}
void SequenceDB::ComputeDistance( const Options & options )
{
	int i, j, N = sequences.size();
	int best_score, best_sum;
	int band_width1, band_left, band_center, band_right;
	int tiden_no, alnln;
	int talign_info[5];
	float distance;
	WorkingBuffer buf( N, options );

	Vector<NVector<float> > dists( N, NVector<float>(N) );

	Sequence comseq( *sequences[0] );

	for(i=0; i<N; i++){
		Sequence *seq = sequences[i];
		char *seqi = seq->data;
		int len = seq->size;
		buf.EncodeWords( seq, options.NAA, false );
		buf.ComputeAAP2( seqi, seq->size );
		dists[i][i] = 0.0;
		if((i+1)%1000 ==0) if( options.verbose)printf( "%9i\n", (i+1) );
		for(j=0; j<i; j++){
			Sequence *rep = sequences[j];
			char *seqj = rep->data;
			int len2 = rep->size;
			band_width1 = (options.band_width < len+len2-2 ) ? options.band_width : len+len2-2;
			diag_test_aapn_est(NAA1, seqj, len, len2, buf, options, best_sum,
					band_width1, band_left, band_center, band_right, 0);
			local_band_align(seqi, seqj, len, len2, mat,
					best_score, tiden_no, alnln, distance, talign_info+1,
					band_left, band_center, band_right, buf, options);
			dists[seq->index][rep->index] = dists[rep->index][seq->index] = distance;
		}
		if (! options.option_r ) break;
		comseq.index = seq->index;
		comseq.size = len;
		for(j=0; j<len; j++) comseq.data[i] = seq->data[len-i-1];
		seqi = comseq.data;
		buf.EncodeWords( &comseq, options.NAA, false );
		buf.ComputeAAP2( seqi, seq->size );
		for(j=0; j<i; j++){
			Sequence *rep = sequences[j];
			char *seqj = rep->data;
			int len2 = rep->size;
			band_width1 = (options.band_width < len+len2-2 ) ? options.band_width : len+len2-2;
			diag_test_aapn_est(NAA1, seqj, len, len2, buf, options, best_sum,
					band_width1, band_left, band_center, band_right, 0);
			local_band_align(seqi, seqj, len, len2, mat,
					best_score, tiden_no, alnln, distance, talign_info+1,
					band_left, band_center, band_right, buf, options);
			if( distance < dists[seq->index][rep->index] )
				dists[seq->index][rep->index] = dists[rep->index][seq->index] = distance;
		}
	}
	string output = options.output + ".dist";
	FILE *fout = fopen( output.c_str(), "w+" );
	fprintf( fout, "1" );
	for(i=1; i<N; i++) fprintf( fout, "\t%i", i+1 );
	fprintf( fout, "\n" );
	for(i=0; i<N; i++){
		fprintf( fout, "%g", dists[i][0] );
		for(j=1; j<N; j++) fprintf( fout, "\t%g", dists[i][j] );
		fprintf( fout, "\n" );
	}
	fclose( fout );
}
void SequenceDB::DoClustering( const Options & options )
{
	int i;
	int NAA = options.NAA;
	double aa1_cutoff = options.cluster_thd;
	double aas_cutoff = 1 - (1-options.cluster_thd)*4;
	double aan_cutoff = 1 - (1-options.cluster_thd)*options.NAA;
	int seq_no = sequences.size();
	int frag_no = seq_no;
	int frag_size = options.frag_size;

#if 0
	ComputeDistance( options );
	return;
#endif

	if( options.threads > 1 ){
		DoClustering( options.threads, options );
		CDHIT_CleanUpTempFiles();
		return;
	}

	if (frag_size){ 
		frag_no = 0;
		for (i=0; i<seq_no; i++) frag_no += (sequences[i]->size - NAA) / frag_size + 1;
	}

	if( ! options.isEST )
		cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd,
				options.tolerance, naa_stat_start_percent, naa_stat, NAA);

	WorkingParam param( aa1_cutoff, aas_cutoff, aan_cutoff );
	WorkingBuffer buffer( frag_no, options );

	WordTable word_table( options.NAA, NAAN );

	size_t mem_need = MinimalMemory( frag_no, buffer.total_bytes, 1, options );
	size_t mem, mega = 1000000;
	size_t mem_limit = (options.max_memory - mem_need) / sizeof(IndexCount);
	int N = sequences.size();

	size_t total_letters = options.total_letters;
	size_t tabsize = 0;

	if( options.verbose)printf( "Table limit with the given memory limit:\n" );
	if( options.verbose)printf( "Max number of representatives: %i\n", MAX_TABLE_SEQ );
	if( options.max_memory ){
		if( options.verbose)printf( "Max number of word counting entries: %li\n", mem_limit );
	}else{
		mem_limit = options.max_entries;
		if( mem_limit > MAX_TABLE_SIZE ) mem_limit = MAX_TABLE_SIZE;
	}
	if( options.verbose)printf( "\n" );

	for(i=0; i<N; ){
		float redundancy = (rep_seqs.size() + 1.0) / (i + 1.0);
		size_t sum = 0;
		int m = i;
		//printf( "\rdefining new sequence group\n", i, m );
		//fflush( stdout );
		if(mem_limit==0)mem_limit=1;
		while( m < N && sum < mem_limit ){
			Sequence *seq = sequences[m];
			if( ! (seq->state & IS_REDUNDANT) ){
				if ( options.store_disk ) seq->SwapIn();
				sum += (size_t)(seq->size * redundancy);
			}
			m ++;
		}
		if( m > N ) m = N;
		if( options.verbose)
		{
			printf( "\rcomparing sequences from  %9i  to  %9i\n", i, m );
			fflush( stdout );
		}
		for(int ks=i; ks<m; ks++){
			Sequence *seq = sequences[ks];
			i = ks + 1;
			if (seq->state & IS_REDUNDANT) continue;
			ClusterOne( seq, ks, word_table, param, buffer, options );
			total_letters -= seq->size;
			if( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
			if( word_table.size >= mem_limit ) break;
			int tmax = MAX_TABLE_SEQ - (frag_size ? seq->size / frag_size + 1 : 0);
			if( (int)word_table.sequences.size() >= tmax || (int)word_table.frag_count >= tmax ) break;
		}
		m = i;
		if( word_table.size == 0 ) continue;
		float p0 = 0;
		for(int j=m; j<N; j++){
			Sequence *seq = sequences[j];
			if (seq->state & IS_REDUNDANT) continue;
			if ( options.store_disk ) seq->SwapIn();
			CheckOne( seq, word_table, param, buffer, options );
			total_letters -= seq->size;
			if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
			int len_bound = param.len_upper_bound;
			if( word_table.sequences[ word_table.sequences.size()-1 ]->size > len_bound ){
				break;
			}
			float p = (100.0*j)/N;
			if( p > p0+1E-1 ){ // print only if the percentage changed
				if( options.verbose)
				{
					printf( "\r%4.1f%%", p );
					fflush( stdout );
				}
				p0 = p;
			}
		}
		if( word_table.size > tabsize ) tabsize = word_table.size;
		//if( i && i < m ) printf( "\r---------- %6i remaining sequences to the next cycle\n", m-i );
		word_table.Clear();
	}
	if( options.verbose)printf( "\n%9li  finished  %9li  clusters\n", (long)sequences.size(), (long)rep_seqs.size() );
	mem = (mem_need + tabsize*sizeof(IndexCount))/mega;
	if( options.verbose)printf( "\nApprixmated maximum memory consumption: %iM\n", (int)mem );
	CDHIT_CleanUpTempFiles();
	word_table.Clear();

#if 0
	int zeros = 0;
	for(i=0; i<word_table.indexCounts.size(); i++) zeros += word_table.indexCounts[i].Size() ==0;
	printf( "%9i  empty entries out of  %9i\n", zeros, word_table.indexCounts.size() );
#endif
}

void SequenceDB::ClusterTo( SequenceDB & other, const Options & options )
{
	int i, flag;
	int len, len_tmp, len_lower_bound, len_upper_bound;
	int NR2_red_no = 0;
	int aan_no = 0;
	char *seqi;
	int NAA = options.NAA;
	double aa1_cutoff = options.cluster_thd;
	double aas_cutoff = 1 - (1-options.cluster_thd)*4;
	double aan_cutoff = 1 - (1-options.cluster_thd)*options.NAA;
	Vector<int>  word_encodes( MAX_SEQ );
	Vector<INTs> word_encodes_no( MAX_SEQ );

	if( ! options.isEST ){
		cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd,
				options.tolerance, naa_stat_start_percent, naa_stat, NAA);
	}

	WorkingParam param( aa1_cutoff, aas_cutoff, aan_cutoff );
	WorkingBuffer buffer( other.sequences.size(), options );

	size_t mem_limit = options.max_memory / sizeof(IndexCount);
	int N = other.sequences.size();
	int M = sequences.size();

	int T = options.threads;
	valarray<size_t>  counts(T);
	Vector<WorkingParam> params(T);
	Vector<WorkingBuffer> buffers(T);
	for(i=0; i<T; i++){
		params[i].Set( aa1_cutoff, aas_cutoff, aan_cutoff );
		buffers[i].Set( N, options );
	}
	if( T >1 ) omp_set_num_threads(T);


	size_t total_letters = 0;
	for(i=0; i<N; i++) total_letters += other.sequences[i]->size;

	WordTable word_table( options.NAA, NAAN );

	for(i=0; i<N; ){
		size_t sum = 0;
		int m = i;
		if(mem_limit==0)mem_limit=1;
		while( m < N && sum < mem_limit && m < (i + MAX_TABLE_SEQ) ){
			Sequence *seq = other.sequences[m];
			if( ! (seq->state & IS_REDUNDANT) ){
				if ( options.store_disk ) seq->SwapIn();
				sum += seq->size;
			}
			m ++;
		}
		if( m > N ) m = N;
		//printf( "m = %i  %i,  %i\n", i, m, m-i );
		for(int ks=i; ks<m; ks++){
			Sequence *seq = other.sequences[ks];
			len = seq->size;
			seqi = seq->data;
			calc_ann_list(len, seqi, NAA, aan_no, word_encodes, word_encodes_no, options.isEST);
			word_table.AddWordCounts(aan_no, word_encodes, word_encodes_no, ks-i, options.isEST);
			word_table.sequences.Append( seq );
			seq->cluster_id = ks;
			seq->state |= IS_REP;
			if ( (ks+1) % 1000 == 0 ) {
				if( options.verbose)
				{
					printf( "." );
					fflush( stdout );
				}
				if ( (ks+1) % 10000 == 0 ) if( options.verbose)printf( "%9i  finished\n", ks+1 );
			}  
		}
		float p0 = 0;
		if( T > 1 ){
			int JM = M;
			counts = 0;
			#pragma omp parallel for schedule( dynamic, 1 )
			for(int j=0; j<JM; j++){
				Sequence *seq = sequences[j];
				if( seq->state & IS_REDUNDANT ) continue;
				int len = seq->size;
				int len_upper_bound = upper_bound_length_rep(len,options);
				int len_lower_bound = len - options.diff_cutoff_aa2;

				int len_tmp = (int) ( ((double)len) * options.diff_cutoff2);
				if (len_tmp < len_lower_bound) len_lower_bound = len_tmp;

				int tid = omp_get_thread_num();
				params[tid].len_upper_bound = len_upper_bound;
				params[tid].len_lower_bound = len_lower_bound;

				if( word_table.sequences[ word_table.sequences.size()-1 ]->size > len_upper_bound ){
					JM = 0;
					continue;
				}

				int flag = other.CheckOne( seq, word_table, params[tid], buffers[tid], options );
				if ((flag == 1) || (flag == -1)) { // if similar to old one delete it
					if (! options.cluster_best) {
						seq->Clear();
						seq->state |= IS_REDUNDANT ;
						counts[tid] ++;
					}
					if( flag == -1 ) seq->state |= IS_MINUS_STRAND; // for EST only
				}
				float p = (100.0*j)/N;
				if( p > p0+1E-1 ){ // print only if the percentage changed
					if( options.verbose)
					{
						printf( "\r%4.1f%%", p );
						fflush( stdout );
					}
					p0 = p;
				}
			}
			for(int j=0; j<T; j++) NR2_red_no += counts[j];
		}else{
			for(int j=0; j<M; j++){
				Sequence *seq = sequences[j];
				if( seq->state & IS_REDUNDANT ) continue;
				len = seq->size;
				seqi = seq->data;
				len_upper_bound = upper_bound_length_rep(len,options);
				len_lower_bound = len - options.diff_cutoff_aa2;

				len_tmp = (int) ( ((double)len) * options.diff_cutoff2);
				if (len_tmp < len_lower_bound) len_lower_bound = len_tmp;
				param.len_upper_bound = len_upper_bound;
				param.len_lower_bound = len_lower_bound;

				if( word_table.sequences[ word_table.sequences.size()-1 ]->size > len_upper_bound ){
					break;
				}

				flag = other.CheckOne( seq, word_table, param, buffer, options );
				if ((flag == 1) || (flag == -1)) { // if similar to old one delete it
					if (! options.cluster_best) {
						seq->Clear();
						seq->state |= IS_REDUNDANT ;
						NR2_red_no ++;
					}
					if( flag == -1 ) seq->state |= IS_MINUS_STRAND; // for EST only
				}
				float p = (100.0*j)/N;
				if( p > p0+1E-1 ){ // print only if the percentage changed
					if( options.verbose)
					{
						printf( "\r%4.1f%%", p );
						fflush( stdout );
					}
					p0 = p;
				}
			}
		}
		if( options.verbose)printf( "\r..........%9i  compared  %9i  clusters\n", i, NR2_red_no );
		word_table.Clear();
		word_table.size = 0;
		i = m;
	}

	if (options.cluster_best) {//delete redundant sequences in options.cluster_best mode
		for (i=0; i<(int)sequences.size(); i++){
			Sequence *seq = sequences[i];
			if (seq->identity > 0 ){
				seq->state |= IS_REDUNDANT;
				NR2_red_no ++;
			}
		}
	}
	for (i=0; i<(int)sequences.size(); i++){
		Sequence *seq = sequences[i];
		if( seq->identity <0 ) seq->identity *= -1;
		if( !(seq->state & IS_REDUNDANT) ) rep_seqs.Append( i );
	}

	if( options.verbose)cout << endl;
	if( options.verbose)cout << sequences.size() << " compared\t" << NR2_red_no << " clustered" << endl;
	CDHIT_CleanUpTempFiles();
}
