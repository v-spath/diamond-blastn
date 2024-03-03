#include "../util/simd.h"

#ifdef __SSE2__

#include <numeric>
#include "../basic/config.h"
#include "../util/log_stream.h"
#include "cluster.h"
#include "../lib/famsa/msa.h"
#include "../util/tsv/file.h"
#include "../dp/dp.h"
#include "../util/sequence/sequence.h"

using std::unique_ptr;
using std::endl;
using std::vector;
using std::shared_ptr;
using std::stringstream;

struct Edge {
	Edge(int64_t node1, int64_t node2) :
		node1(node1),
		node2(node2) {
	}
	int64_t node1, node2, cc;
	bool operator<(const Edge& e) const {
		return cc < e.cc;
	}
	struct CC {
		int64_t operator()(const Edge& edge) const {
			return edge.cc;
		}
	};
};

struct Cfg {
	Cfg(bool lazy_titles, const FlatArray<OId>& clusters, const vector<OId>& centroids, SequenceFile& db) :
		lazy_titles(lazy_titles),
		clusters(clusters),
		centroids(centroids),
		db(db)
	{}
	const bool lazy_titles;
	const FlatArray<OId>& clusters;
	const vector<OId>& centroids;
	SequenceFile& db;
	shared_ptr<Block> seqs;
};

namespace Cluster {

static int64_t root(vector<int64_t>& v, int64_t i) {
	const int64_t p = v[i];
	if (p == i)
		return i;
	const int64_t pp = v[p];
	if (pp == p)
		return p;
	const int64_t r = root(v, p);
	v[i] = r;
	return r;
}

std::map<int64_t, CProfile> profs;

static void compute_msa(int64_t cluster, Cfg& cfg) {
	vector<CSequence> seqs;
	seqs.reserve(cfg.clusters.count(cluster));
	int64_t n = 0;
	for (auto it = cfg.clusters.cbegin(cluster); it != cfg.clusters.cend(cluster); ++it) {
		seqs.emplace_back("", cfg.seqs->seqs()[*it].to_string(), n++);
	}

	CParams params;
	params.n_threads = 1;
	CFAMSA famsa(params);
	if (!famsa.ComputeMSA(seqs))
		throw std::runtime_error("Error computing MSA");
	vector<CGappedSequence*> msa;
	famsa.GetAlignment(msa);
	/*for (CGappedSequence* s : msa)
		std::cout << s->Decode() << endl;*/

	CProfile* p = famsa.get_final_profile();
	p->card = p->data.size();
	p->data.clear();
	profs[cluster] = *p;
	//famsa.get_final_profile()->AlignProfProfLocal(p, p, nullptr, nullptr);
}

static void merge_clusters(vector<OId>& clust, OId rep1, OId rep2, vector<int>& indirection, int new_indirection) {
	for (int64_t i = 0; i < clust.size(); ++i)
		if (clust[i] == rep2) {
			clust[i] = rep1;
			indirection[i] += new_indirection;
		}
}

static double mutual_cov(OId id1, OId id2, Block& seqs) {
	const HspValues hsp_values = HspValues::COORDS;
	DP::Targets dp_targets;
	Statistics stats;
	const Sequence seq1 = seqs.seqs()[id1], seq2 = seqs.seqs()[id2];
	const int bin = DP::BandedSwipe::bin(hsp_values, seq1.length(), 0, 0, (int64_t)seq1.length() * (int64_t)seq2.length(), 0, 0);
	dp_targets[bin].emplace_back(seq2, seq2.length(), 0);
	const Bias_correction cbs(seq1);
	DP::Params p{ seq1, "", Frame(0), seq1.length(), config.comp_based_stats == 1 ? cbs.int8.data() : nullptr, DP::Flags::FULL_MATRIX, hsp_values, stats, nullptr};
	list<Hsp> hsps = DP::BandedSwipe::swipe(dp_targets, p);
	if (hsps.empty())
		return 0;
	return std::min(hsps.front().query_cover_percent(seq1.length()), hsps.front().subject_cover_percent(seq2.length()));
}

void profile_recluster() {
	config.database.require();
	config.clustering.require();
	score_matrix.set_db_letters(1000000);
	if (config.query_file.empty())
		throw runtime_error("--query is missing");

	TaskTimer timer("Opening the database");
	unique_ptr<SequenceFile> db(SequenceFile::auto_create({ config.database }, SequenceFile::Flags::NEED_LETTER_COUNT | SequenceFile::Flags::ACC_TO_OID_MAPPING | SequenceFile::Flags::OID_TO_ACC_MAPPING, SequenceFile::Metadata()));
	score_matrix.set_db_letters(config.db_size ? config.db_size : db->letters());
	timer.finish();
	message_stream << "#Database sequences: " << db->sequence_count() << ", #Letters: " << db->letters() << endl;
	if (flag_any(db->format_flags(), SequenceFile::FormatFlags::TITLES_LAZY))
		db->init_random_access(0, 0, false);

	FlatArray<OId> clusters;
	vector<OId> centroids;
	tie(clusters, centroids) = read<OId>(config.clustering, *db, CentroidSorted());
	message_stream << "Found " << centroids.size() << " centroids, " << clusters.data_size() << " mappings in input file." << endl;

	Cfg cfg{ flag_any(db->format_flags(), SequenceFile::FormatFlags::TITLES_LAZY), clusters, centroids, *db };

	SequenceFile::LoadFlags flags = SequenceFile::LoadFlags::SEQS | SequenceFile::LoadFlags::CONVERT_ALPHABET | SequenceFile::LoadFlags::NO_CLOSE_WEAKLY;
	if (!cfg.lazy_titles)
		flags |= SequenceFile::LoadFlags::TITLES;

	timer.go("Building clustering");
	vector<OId> clust = member2centroid_mapping(clusters, centroids);

	timer.go("Loading database sequences");
	cfg.seqs.reset(db->load_seqs(INT64_MAX, nullptr, flags));

	if (false) {
		timer.go("Loading accession mapping");
		const Util::Tsv::Table acc_mapping = db->seqid_file().read(config.threads_);
		timer.go("Computing profiles");
		for (int64_t i = 0; i < centroids.size(); ++i) {
			compute_msa(i, cfg);
		}

		atomic<int64_t> next;
		auto f = [&] {
			for (;;) {
				const int64_t i = next++;
				if (i >= centroids.size())
					break;
				for (int64_t j = i + 1; j < centroids.size(); ++j) {
					const CProfile::LocalAln aln = profs[i].AlignProfProfLocal(&profs[i], &profs[j], nullptr, nullptr);
					const double cov1 = double(aln.i1 - aln.i0) / profs[i].width;
					const double cov2 = double(aln.j1 - aln.j0) / profs[j].width;
					if (cov1 > 0.5 && cov2 > 0.5) {
						stringstream ss;
						ss << acc_mapping[i].template get<string>(0) << '\t' << acc_mapping[j].template get<string>(0) << '\t'
							<< cov1 << '\t' << cov2 << endl;
						cout << ss.str();
					}
				}
			}
			};

		vector<thread> t;
		for (int i = 0; i < config.threads_; ++i)
			t.emplace_back(f);
		for (auto& i : t)
			i.join();
	}

	timer.go("Loading alignments");
	Util::Tsv::File aln_file(Util::Tsv::Schema{ Util::Tsv::Type::STRING, Util::Tsv::Type::STRING }, config.query_file.front());
	const Util::Tsv::Table alns = aln_file.read(config.threads_);
	timer.finish();
	message_stream << "#Alignments: " << alns.size() << endl;

	timer.go("Reading edges");
	atomic<int64_t> next(0);
	mutex mtx;
	vector<Edge> edges;
	auto f = [&]{
		static const int64_t N = 1000000;
		for (;;) {
			const int64_t begin = next.fetch_add(N), end = std::min(begin + N, alns.size());
			if (end <= begin)
				break;
			vector<Edge> v;
			v.reserve(N);
			for (int64_t i = begin; i < end; ++i) {
				const int64_t member1 = db->accession_to_oid(alns[i].get<string>(0)).front();
				const int64_t member2 = db->accession_to_oid(alns[i].get<string>(1)).front();
				const int64_t rep1 = clust[member1], rep2 = clust[member2];
				if (rep1 != rep2) {
					v.emplace_back(member1, member2);
				}
			}
			{
				std::lock_guard<mutex> lock(mtx);
				edges.insert(edges.end(), v.begin(), v.end());
			}
		}
	};
	vector<thread> threads;
	for (int i = 0; i < config.threads_; ++i)
		threads.emplace_back(f);
	for (auto& t : threads)
		t.join();
	timer.finish();
	message_stream << "#Edges: " << edges.size() << endl;

	timer.go("Finding connected components");
	vector<int64_t> parent(db->sequence_count());
	std::iota(parent.begin(), parent.end(), 0);
	for (const Edge& e : edges) {
		const int64_t r1 = root(parent, clust[e.node1]), r2 = root(parent, clust[e.node2]);
		parent[r2] = r1;
	}

	timer.go("Preparing edges");
	for (Edge& e : edges)
		e.cc = root(parent, clust[e.node1]);
	parent.clear();
	parent.shrink_to_fit();

	timer.go("Sorting edges");
	ips4o::parallel::sort(edges.begin(), edges.end(), std::less<Edge>(), config.threads_);

	timer.go("Indexing connected components");
	vector<int64_t> ccs;
	ccs.push_back(0);
	auto it = merge_keys(edges.cbegin(), edges.cend(), Edge::CC());
	while (it.good()) {
		ccs.push_back(std::distance(edges.cbegin(), it.end()));
		++it;
	}
	timer.finish();
	message_stream << "#Connected components: " << ccs.size() - 1 << endl;

	timer.go("Processing alignments");
	vector<int> indirection(db->sequence_count(), 0);

	next = 0;
	atomic<int> merges(0);
	auto f2 = [&] {
		for (;;) {
			const int64_t cc = next++;
			if (cc >= ccs.size() - 1)
				break;
			const auto begin = edges.cbegin() + ccs[cc], end = edges.cbegin() + ccs[cc + 1];
			for (auto it = begin; it != end; ++it) {
				const int64_t rep1 = clust[it->node1], rep2 = clust[it->node2];
				if (rep1 != rep2) {
					//cout << Util::Seq::seqid(cfg.seqs->ids()[rep1], false) << '\t' << Util::Seq::seqid(cfg.seqs->ids()[rep2], false) << '\t'  << mutual_cov(rep1, rep2, *cfg.seqs) << endl;
					/*if (mutual_cov(rep1, rep2, *cfg.seqs) < 50)
						continue;
					if (mutual_cov(rep1, it->node1, *cfg.seqs) < 50)
						continue;
					if (mutual_cov(rep1, it->node2, *cfg.seqs) < 50)
						continue;
					if (mutual_cov(rep2, it->node1, *cfg.seqs) < 50)
						continue;
					if (mutual_cov(rep2, it->node2, *cfg.seqs) < 50)
						continue;
					if (mutual_cov(it->node1, it->node2, *cfg.seqs) < 50)
						continue;*/
					merge_clusters(clust, rep1, rep2, indirection, indirection[it->node1] + 1);
					++merges;
				}
			}
		}
	};

	threads.clear();
	for (int i = 0; i < config.threads_; ++i)
		threads.emplace_back(f2);
	for (auto& t : threads)
		t.join();
	timer.finish();	
	//if (rep1 != rep2 && indirection[member1] + indirection[member2] <= config.max_indirection) {
	//cout << s1.front() << '\t' << s2.front() << endl;
	//cout << alns[i].get<string>(0) << '\t' << alns[i].get<string>(1) << endl;

	message_stream << "Cluster count: " << clusters.size() - merges << endl;
	
	timer.go("Writing output");
	Util::Tsv::File out_file(Util::Tsv::Schema{ Util::Tsv::Type::STRING, Util::Tsv::Type::STRING }, config.output_file, Util::Tsv::Flags::WRITE);
	output_mem(out_file, *db, clust);
}

}

#endif