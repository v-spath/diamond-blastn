#include <iostream>
#include "profile.h"
#include <assert.h>

using std::cout;
using std::endl;

CProfile::LocalAln CProfile::AlignProfProfLocal(CProfile* profile1, CProfile* profile2, vector<int>* column_mapping1, vector<int>* column_mapping2)
{
	size_t prof1_width = profile1->width;
	size_t prof2_width = profile2->width;

	size_t prof1_card = profile1->card;
	size_t prof2_card = profile2->card;

	score_t gap_open = params->gap_open;
	score_t gap_ext = params->gap_ext;
	score_t gap_term_open = params->gap_term_open;
	score_t gap_term_ext = params->gap_term_ext;

	CDPMatrix matrix(prof1_width + 1, prof2_width + 1);
	vector<vector<bool>> zeros(prof1_width + 1, vector<bool>(prof2_width + 1));
	matrix.set_zeros(params->instruction_set);

	dp_row_t curr_row(prof2_width + 1);
	dp_row_t prev_row(prof2_width + 1);

	// Frequency of symbols at all columns of profile1 and profile2
//	vector<pair<int, int>> col1(32);
	array<pair<int, int>, 32> col1;

	CProfileValues<score_t, NO_SYMBOLS>& scores1 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile1->scores);
	CProfileValues<score_t, NO_SYMBOLS>& scores2 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile2->scores);

	// Prepare ranges for computed cells in rows
	bool is_guided = column_mapping1 != nullptr && column_mapping2 != nullptr;
	vector<pair<int, int>> row_ranges;

	if (is_guided)
		FindRowRanges(column_mapping1, column_mapping2, row_ranges);
	else
		row_ranges.assign(prof1_width + 1, make_pair(0, prof2_width));

	// Precompute scores for gaps for profile2
	vector<dp_gap_costs> prof2_gaps(prof2_width + 1);

	//size_t n_gap_open = 0;
	//size_t n_gap_ext = 0;
	//size_t n_gap_term_open = 0;
	//size_t n_gap_term_ext = 0;

	counter_t n_gap_prof1_start_open, n_gap_prof1_start_ext, n_gap_prof1_start_term_open, n_gap_prof1_start_term_ext,
		n_gap_prof1_cont_ext, n_gap_prof1_cont_term_ext;
	vector<dp_gap_corrections> gap_corrections(prof2_width + 1);

	//scoresX.get_value(j, GAP_OPEN) returns information of how the column j of the profile X aligns with a single gap_open
	//scoresX.get_value(j, GAP_EXT) returns information of how the column j of the profile X aligns with a single gap_ext
	//etc.
	//The task is to calculate how the column j of the profile X aligns with a column of gaps of different categories

	for (size_t j = 0; j <= prof2_width; ++j)
	{
		prof2_gaps[j].open = scores2.get_value(j, GAP_OPEN);
		prof2_gaps[j].ext = scores2.get_value(j, GAP_EXT);
		prof2_gaps[j].term_open = scores2.get_value(j, GAP_TERM_OPEN);
		prof2_gaps[j].term_ext = scores2.get_value(j, GAP_TERM_EXT);
	}

	// Boundary conditions
	prev_row[0].D = 0;
	prev_row[0].H = -infty;
	prev_row[0].V = -infty;

	if (prof2_width >= 1)
	{
		prev_row[1].D = -infty;

		prev_row[1].H = prev_row[0].D + prof2_gaps[1].term_open * prof1_card;

		prev_row[1].V = -infty;
		matrix.set_dir_all(0, 1, direction_t::H);
	}

	for (size_t j = 2; j <= prof2_width; ++j)
	{
		prev_row[j].D = -infty;

		prev_row[j].H = prev_row[j - 1].H + prof2_gaps[j].term_ext * prof1_card;

		prev_row[j].V = -infty;
		matrix.set_dir_all(0, j, direction_t::H);
	}
	prev_row[prof2_width].H = -infty;

	// Precomputing for gap correction in the profile2
	vector<score_t> n_gaps_prof2_to_change(prof2_width + 1);
	vector<score_t> n_gaps_prof2_term_to_change(prof2_width + 1);
	vector<score_t> gaps_prof2_change(prof2_width + 1);

	for (size_t j = 1; j <= prof2_width; ++j)
	{
		DP_SolveGapsProblemWhenStarting(j, prof2_width, prof2_card, profile2,
			gap_corrections[j].n_gap_start_open, gap_corrections[j].n_gap_start_ext, gap_corrections[j].n_gap_start_term_open, gap_corrections[j].n_gap_start_term_ext);

		DP_SolveGapsProblemWhenContinuing(j, prof2_width, prof2_card, profile2,
			gap_corrections[j].n_gap_cont_ext, gap_corrections[j].n_gap_cont_term_ext);

#ifndef NO_GAP_CORRECTION
		n_gaps_prof2_to_change[j] = profile2->counters.get_value(j, GAP_OPEN);      //the number of gaps to be changed from gap_open into gap_ext 
		n_gaps_prof2_term_to_change[j] = profile2->counters.get_value(j, GAP_TERM_OPEN); //the number of gaps to be changed from term_open into term_ext

#else
		n_gaps_prof2_to_change[j] = 0;
		n_gaps_prof2_term_to_change[j] = 0;
#endif

		gaps_prof2_change[j] = n_gaps_prof2_to_change[j] * (gap_ext - gap_open) +
			n_gaps_prof2_term_to_change[j] * (gap_term_ext - gap_term_open);
	}

	score_t global_max = std::numeric_limits<score_t>::min();
	size_t global_max_i = -1, global_max_j = -1;
	// Calculate matrix interior
	for (size_t i = 1; i <= prof1_width; ++i)
	{
		// Precompute scores for gaps for current and previous row of profile1
		//score_t prof1_gap_open_prev = scores1.get_value(i - 1, GAP_OPEN);
		score_t prof1_gap_open_curr = scores1.get_value(i, GAP_OPEN);
		//score_t prof1_gap_term_open_prev = scores1.get_value(i - 1, GAP_TERM_OPEN);
		score_t prof1_gap_term_open_curr = scores1.get_value(i, GAP_TERM_OPEN);
		score_t prof1_gap_ext_curr = scores1.get_value(i, GAP_EXT);
		score_t prof1_gap_term_ext_curr = scores1.get_value(i, GAP_TERM_EXT);

		// Boundary conditions
		curr_row[0].D = -infty;
		curr_row[0].H = -infty;
		matrix.set_dir_all(i, 0, direction_t::V);

		if (row_ranges[i].first)
			curr_row[row_ranges[i].first - 1].D = curr_row[row_ranges[i].first - 1].H = curr_row[row_ranges[i].first - 1].V = -infty;

		if (i < prof1_width)
		{
			if (i == 1)
				curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_open_curr * prof2_card;
			else
				curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_ext_curr * prof2_card;

			for (int j = row_ranges[i].second + 1; j <= min(row_ranges[i + 1].second, (int)prof2_width); ++j)
				curr_row[j].D = curr_row[j].H = curr_row[j].V = -infty;
		}
		else
			curr_row[0].V = -infty;

		// Calculate frequency of symbols in (i-1)-th column of profile1
/*		col1.clear();
		size_t col1_n_non_gaps = 0;
		for (size_t k = 0; k < NO_AMINOACIDS_AND_GAPS; ++k)
			if (profile1->counters.get_value(i, k))
			{
				size_t count = profile1->counters.get_value(i, k);
				col1.emplace_back(k, count);
				if (k < NO_AMINOACIDS)
					col1_n_non_gaps += count;
			}
		size_t col1_size = col1.size();
		*/

		size_t col1_size = 0;
		size_t col1_n_non_gaps = 0;
		for (size_t k = 0; k < NO_AMINOACIDS_AND_GAPS; ++k)
			if (profile1->counters.get_value(i, k))
			{
				size_t count = profile1->counters.get_value(i, k);
				col1[col1_size++] = make_pair(k, count);
				if (k < NO_AMINOACIDS)
					col1_n_non_gaps += count;
			}

		// Precomputing for gap correction in the profile1
		n_gap_prof1_start_open = n_gap_prof1_start_ext = n_gap_prof1_start_term_open = n_gap_prof1_start_term_ext = 0;
		DP_SolveGapsProblemWhenStarting(i, prof1_width, prof1_card, profile1,
			n_gap_prof1_start_open, n_gap_prof1_start_ext, n_gap_prof1_start_term_open, n_gap_prof1_start_term_ext);

		DP_SolveGapsProblemWhenContinuing(i, prof1_width, prof1_card, profile1,
			n_gap_prof1_cont_ext, n_gap_prof1_cont_term_ext);

#ifndef NO_GAP_CORRECTION
		size_t n_gaps_prof1_to_change = profile1->counters.get_value(i, GAP_OPEN);   //the number of gaps to be changed from gap_open into gap_ext
		size_t n_gaps_prof1_term_to_change = profile1->counters.get_value(i, GAP_TERM_OPEN); //the number of gaps to be changed from term_open into term_ext
#else
		size_t n_gaps_prof1_to_change = 0;
		size_t n_gaps_prof1_term_to_change = 0;
#endif

		unsigned char* matrix_cell = matrix.get_cell(i, max(1, row_ranges[i].first) - 1);
		size_t max_j = min(row_ranges[i].second, (int)prof2_width);
		score_t row_max = numeric_limits<score_t>::min();
		size_t row_max_j = -1;
		for (size_t j = max(1, row_ranges[i].first); j <= max_j; ++j)
		{
			// Get current cell of matrix for faster access
			++matrix_cell;

			//...........................................................................
			// Calcualte score for D
			//...........................................................................

			//analyzing an alignment of the column j of the profile 2 with the column i of the profile 1
			score_t* scores2_column = scores2.get_column(j);

			//scores2_column contains information on how the column j aligns with a single character that could appear in the profile 1
			//the variable t is to store a cost of the alignment of the column j with a specific character (col1[k].first) 
			//of a concrete number of occurrences (col1[k].second) of the profile 1

			score_t t = 0;
			switch (col1_size)
			{
			case 3:
				t += col1[2].second * scores2_column[col1[2].first];
			case 2:
				t += col1[1].second * scores2_column[col1[1].first];
			case 1:
				t += col1[0].second * scores2_column[col1[0].first];
			case 0:
				break;
			default:
				for (size_t k = 0; k < col1_size; ++k)
					t += col1[k].second * scores2_column[col1[k].first];
			}

			//an alignment of a column with another column, after an an alignment of two columns  
			score_t t_D = prev_row[j - 1].D + t;

			//..........................................................................
			//an alignment of a column with another column, after inserting a column of gaps before, into the profile 1
			score_t t_H = prev_row[j - 1].H;

			//One has to check if there is a need of exchanging gap_opens/gap_term_opens (from the column of the profile 1 that is being just aligned) 
			//into gap_extensions/gap_term_extensions
			//If so, the value of the variable t has to be corrected, because it was computed with:
			//t += col1[k].liczba_gap_open * scores2_column[GAP_OPEN];
			//t += col1[k].liczba_gap_ext * scores2_column[GAP_EXT];
			//t += col1[k].liczba_gap_term_open * scores2_column[GAP_TERM_OPEN];
			//t += col1[k].liczba_gap_term_ext * scores2_column[GAP_TERM_EXT];

			if (n_gaps_prof1_to_change || n_gaps_prof1_term_to_change)
			{
				//correction by delta_t
				score_t delta_t = n_gaps_prof1_to_change * (scores2_column[GAP_EXT] - scores2_column[GAP_OPEN]) +
					n_gaps_prof1_term_to_change * (scores2_column[GAP_TERM_EXT] - scores2_column[GAP_TERM_OPEN]);

				t_H = t_H + t + delta_t;
			}
			else
				t_H = t_H + t;

			//..........................................................................
			//an alignment of a column with another column, after inserting a column of gaps before, into the profile 2
			score_t t_V = prev_row[j - 1].V;

			t_V += t + gaps_prof2_change[j] * col1_n_non_gaps;

			if (t_D > t_H && t_D > t_V)
			{
				curr_row[j].D = t_D;
				matrix.set_dir_D(matrix_cell, direction_t::D);
			}
			else if (t_H > t_V)
			{
				curr_row[j].D = t_H;
				matrix.set_dir_D(matrix_cell, direction_t::H);
			}
			else
			{
				curr_row[j].D = t_V;
				matrix.set_dir_D(matrix_cell, direction_t::V);
			}

			if (curr_row[j].D < 0) {
				curr_row[j].D = 0;
				zeros[i][j] = true;
			}
			//cout << "i=" << i << " j=" << j << " D=" << curr_row[j].D << endl;
			if (curr_row[j].D > row_max) {
				row_max = curr_row[j].D;
				row_max_j = j;
			}
			
			

			//...........................................................................
			// Calcualte score for H
			//...........................................................................
			//inserting a columns of gaps into the profile 1

			//the case when a column (i) was aligned with the column (j-1)
			score_t gap_corr = prof2_gaps[j].open * n_gap_prof1_start_open + prof2_gaps[j].ext * n_gap_prof1_start_ext +
				prof2_gaps[j].term_open * n_gap_prof1_start_term_open + prof2_gaps[j].term_ext * n_gap_prof1_start_term_ext;

			t_D = curr_row[j - 1].D + gap_corr;

			//the case when a sequence of gaps is continued (another columnn of gaps is inserted)
			t_H = curr_row[j - 1].H + prof2_gaps[j].ext * n_gap_prof1_cont_ext +
				prof2_gaps[j].term_ext * n_gap_prof1_cont_term_ext;

#ifdef ALWAYS_3_DIRS
			//the column of gaps is inserted into the profile 1 after a column of gaps in the profile 2
			if ((i > 1) && (j > 1))
			{
				t_V = curr_row[j - 1].V + gap_corr;

				if (t_D > t_H && t_D > t_V)
				{
					curr_row[j].H = t_D;
					matrix.set_dir_H(matrix_cell, direction_t::D);
				}
				else if (t_V > t_H)	//Attention!!! Swapping the checking order 
				{
					curr_row[j].H = t_V;
					matrix.set_dir_H(matrix_cell, direction_t::V);
				}
				else
				{
					curr_row[j].H = t_H;
					matrix.set_dir_H(matrix_cell, direction_t::H);
				}
			}
			else
			{
				if (t_D > t_H)
				{
					curr_row[j].H = t_D;
					matrix.set_dir_H(matrix_cell, direction_t::D);
				}
				else
				{
					curr_row[j].H = t_H;
					matrix.set_dir_H(matrix_cell, direction_t::H);
				}
			}
#else // !ALWAYS_3_DIRS
			if (t_D > t_H)
			{
				curr_row[j].H = t_D;
				matrix.set_dir_H(matrix_cell, direction_t::D);
			}
			else
			{
				curr_row[j].H = t_H;
				matrix.set_dir_H(matrix_cell, direction_t::H);
			}
#endif // !ALWAYS_3_DIRS
			//...........................................................................
			// Calcualte score for V
			//...........................................................................
			//inserting a columns of gaps into the profile 2

			//the case when a column (i-1) was aligned with the column (j)
			gap_corr = prof1_gap_open_curr * gap_corrections[j].n_gap_start_open + prof1_gap_ext_curr * gap_corrections[j].n_gap_start_ext +
				prof1_gap_term_open_curr * gap_corrections[j].n_gap_start_term_open + prof1_gap_term_ext_curr * gap_corrections[j].n_gap_start_term_ext;
			t_D = prev_row[j].D + gap_corr;

			//the case when a sequence of gaps is continued (another columnn of gaps is inserted)
			t_V = prev_row[j].V + prof1_gap_ext_curr * gap_corrections[j].n_gap_cont_ext +
				prof1_gap_term_ext_curr * gap_corrections[j].n_gap_cont_term_ext;

#ifdef ALWAYS_3_DIRS
			if ((i > 1) && (j > 1))
			{
				t_H = prev_row[j].H + gap_corr;

				if (t_D > t_H && t_D > t_V)
				{
					curr_row[j].V = t_D;
					matrix.set_dir_V(matrix_cell, direction_t::D);
				}
				else if (t_H > t_V)
				{
					curr_row[j].V = t_H;
					matrix.set_dir_V(matrix_cell, direction_t::H);
				}
				else
				{
					curr_row[j].V = t_V;
					matrix.set_dir_V(matrix_cell, direction_t::V);
				}
			}
			else
			{
				if (t_D > t_V)
				{
					curr_row[j].V = t_D;
					matrix.set_dir_V(matrix_cell, direction_t::D);
				}
				else
				{
					curr_row[j].V = t_V;
					matrix.set_dir_V(matrix_cell, direction_t::V);
				}
			}
#else
			if (t_D > t_V)
			{
				curr_row[j].V = t_D;
				matrix.set_dir_V(matrix_cell, direction_t::D);
			}
			else
			{
				curr_row[j].V = t_V;
				matrix.set_dir_V(matrix_cell, direction_t::V);
			}
#endif //ALWAYS_3_DIRS

		}
		//cout << "i=" << i << " maxD=" << row_max << " maxj=" << row_max_j << endl;
		if (row_max > global_max) {
			global_max = row_max;
			global_max_i = i;
			global_max_j = row_max_j;
		}

		curr_row.swap(prev_row);
	}
	//cout << "i=" << global_max_i << " global_maxD=" << global_max << " maxj=" << global_max_j << endl;

	// Construct alignment
	const pair<size_t, size_t> bt = Backtrace(profile1, profile2, matrix, zeros, global_max_i, global_max_j);
	return { bt.first, global_max_i, bt.second, global_max_j };
}

std::pair<size_t, size_t> CProfile::Backtrace(CProfile* profile1, CProfile* profile2, CDPMatrix& matrix, const std::vector<std::vector<bool>>& zeros, size_t max_i, size_t max_j, uint32_t no_threads)
{
	size_t i = max_i;
	size_t j = max_j;

	direction_t dir;
	dir = direction_t::D;

	while ((i || j) && !zeros[i][j])
	{
		//cout << "i=" << i << " j=" << j << endl;
		if (dir == direction_t::D)
		{
			dir = matrix.get_dir_D(i--, j--);
		}
		else if (dir == direction_t::H)
		{
			dir = matrix.get_dir_H(i, j--);
		}
		else if (dir == direction_t::V)
		{
			dir = matrix.get_dir_V(i--, j);
		}
		else
		{
			assert(0);
		}
	}
	return { i,j };
}