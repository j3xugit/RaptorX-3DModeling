#pragma once
#include "cdhit_util.h"
#include "cdhit_scomat.h"
#include "cdhit_option.h"
#include "cdhit_param.h"
#include "cdhit_seq.h"
#include "cdhit_table.h"
#include "cdhit_buffer.h"
using namespace std;

//====== class: SequenceDB ======//start
class SequenceDB
{
	public:
		//---- constructor and destructor ----//
		SequenceDB(){ }
		~SequenceDB(){ Clear(); }
		void Clear(){
			for(int i=0; i<(int)sequences.size(); i++) delete sequences[i];
			sequences.clear(); rep_seqs.clear();
		}

		//---- basic variables -----//
		int NAAN;
		Vector<Sequence*>  sequences;
		Vector<int>        rep_seqs;

		//---- basic functions ----//
		//-> IO
		void Read( const char *file, const Options & options );
		void WriteClusters( const char *db, const char *newdb, const Options & options );
		void WriteExtra1D( const Options & options );
		void WriteExtra2D( SequenceDB & other, const Options & options );
		void DivideSave( const char *db, const char *newdb, int n, const Options & options );

		//-> simp IO  //- by Sheng Wang -//__140230__//
		void Read_From_String( const Options & options,
			vector <string> & in_head, vector <string> & in_seq );
		void WriteExtra1D_To_String( const Options & options,
			vector <vector <string> > & member, vector <int >& center );

		//-> swap
		void SwapIn( int seg, bool reponly=false );
		void SwapOut( int seg );

		//-> table
		void SortDivide( Options & options, bool sort=true );
		void MakeWordTable( const Options & options );
		size_t MinimalMemory( int frag_no, int bsize, int T, const Options & options );
		void ComputeDistance( const Options & options );

		//-> main 
		void ClusterOne( Sequence *seq, int id, WordTable & table,
				WorkingParam & param, WorkingBuffer & buf, const Options & options );
		void DoClustering( const Options & options );
		void DoClustering( int T, const Options & options );
		void ClusterTo( SequenceDB & other, const Options & Options );
		int  CheckOne( Sequence *seq, WordTable & tab, WorkingParam & par, WorkingBuffer & buf, const Options & opt );
		int  CheckOneEST( Sequence *seq, WordTable & tab, WorkingParam & par, WorkingBuffer & buf, const Options & opt );
		int  CheckOneAA( Sequence *seq, WordTable & tab, WorkingParam & par, WorkingBuffer & buf, const Options & opt );


	//----- other inner functions -----//
	public:
		//--- inner variables ---//
		ScoreMatrix mat;
		Vector<int> Comp_AAN_idx;
		
		//--- inner functions ---//
		//-> main func
		int diag_test_aapn(int NAA1, char iseq2[], int len1, int len2, WorkingBuffer & buffer, const Options & options,
			int &best_sum, int band_width, int &band_left, int &band_center, int &band_right, int required_aa1);
		int diag_test_aapn_est(int NAA1, char iseq2[], int len1, int len2, WorkingBuffer & buffer, const Options & options,
			int &best_sum, int band_width, int &band_left, int &band_center, int &band_right, int required_aa1);
		int local_band_align( char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix &mat, 
			int &best_score, int &iden_no, int &alnln, float &dist, int *alninfo,
			int band_left, int band_center, int band_right, WorkingBuffer & buffer, const Options & options);

		//-> vice func
		void make_comp_iseq(int len, char *iseq_comp, char *iseq);
		void make_comp_short_word_index(int NAA, int *NAAN_array, Vector<int> &Comp_AAN_idx);
		int upper_bound_length_rep(int len, double opt_s, int opt_S, double opt_aL, int opt_AL );
		int upper_bound_length_rep(int len, const Options & options );
		void cal_aax_cutoff(double &aa1_cutoff, double &aa2_cutoff, double &aan_cutoff,
			double cluster_thd, int tolerance, int naa_stat_start_percent,int naa_stat[5][61][4], int NAA);
		void update_aax_cutoff(double &aa1_cutoff, double &aa2_cutoff, double &aan_cutoff,
			int tolerance, int naa_stat_start_percent,int naa_stat[5][61][4], int NAA, double cluster_thd);
		int calc_ann_list(int len, char *seqi, int NAA, int& aan_no, Vector<int> & aan_list, Vector<INTs> & aan_list_no, bool est);

};
//====== class: SequenceDB ======//end
