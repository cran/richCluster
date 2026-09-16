//
//  RichCluster.cpp
//  RichCluster
//
//  Created by Sarah on 6/2/25.
//

#include <stdio.h>
#include <string>
#include <algorithm>
#include <numeric>
#include "RichCluster.h"
#include "StringUtils.h"
#include <Rcpp.h>

void richCluster::computeDistances() {
  if (verbose) Rcpp::Rcout << "Computing distances..." << std::endl;
  int totalGeneCount = StringUtils::countUniqueElements(geneIDs, geneDelim);
  
  for (int i=0; i<n_terms; ++i) {
    // D2: before this line nothing under src/ was interruptible -- measured,
    // four SIGINTs over 8 s did not stop a 500-term run, and the whole call
    // scales about O(n^3.8), so a 2000-term run was a session the user could
    // only leave with kill -9.  Once per OUTER iteration only: the inner loop
    // below runs n_terms times per check, and R_ToplevelExec (what
    // checkUserInterrupt costs) per pair would dominate the distance it guards.
    Rcpp::checkUserInterrupt();

    // the unordered set of term1 genes
    std::unordered_set<std::string> term1_genes = StringUtils::splitStringToUnorderedSet(geneIDs[i], geneDelim);
    
    for (int j=0; j<n_terms; ++j) {
      if (i == j) {
        distMatrix.setDistance(richCluster::SAME_TERM_DISTANCE, i, j);
        continue;
      }
      // unordered set of term2 genes
      std::unordered_set<std::string> term2_genes = StringUtils::splitStringToUnorderedSet(geneIDs[j], geneDelim);
      
      double distanceScore = dm.computeDistance(term1_genes, term2_genes, totalGeneCount);
      distMatrix.setDistance(distanceScore, i, j);
      
      // DS-11: cutoff comparisons are STRICT package-wide.  This site used to
      // read `>=` while the merge predicate below (findBestMergePartner) and
      // all four DavidClustering sites used `>`, so at an exact tie with the
      // cutoff the package disagreed with itself.  Unified on `>`: it is the
      // majority convention (5 sites to 2), it is what NEWS.md already tells
      // users ("a pair is merged only on a strict inequality"), and
      // test-david-fixed-point.R forbids tightening the DAVID path to `>=`.
      // Measured on the bundled data at cutoff 0.5: kappa has 0 exact ties, so
      // the default path is unchanged; jaccard has 6 and dice 26.
      if (distanceScore > dm.getCutoff()) {
        // add to adjacency list bidirectionally
        adjList.addNeighbor(i, j);
        adjList.addNeighbor(j, i);
      }
    }
  }
  if (verbose) Rcpp::Rcout << "Done filling out DistanceMatrix." << std::endl;
}

// go through adjacency list and find the best subset of each seed
void richCluster::filterSeeds() {
  if (verbose) Rcpp::Rcout << "Filtering seeds..." << std::endl;
  
  // SPEC-RC-011 REQ-011-005: visit nodes in ascending index order so results do
  // not depend on the standard library's hash traversal (build portability --
  // tools/baseline/README.md finding 4). Keys 0..n_terms-1 are all pre-inserted
  // by the AdjacencyList constructor (src/AdjacencyList.h:16-20), so .at() cannot
  // throw.
  //
  // APPLIED under the author ruling of 2026-08-23, which resolved the plan.md
  // section F.3 stop rule that fired on 2026-08-23: this swap is NOT
  // output-neutral on this build. libstdc++ traverses this unordered_map in
  // DESCENDING key order (probe: n=580 -> first=579; n=1866 -> first=1865), and
  // mergeClusters() is emission-order dependent, so canonicalising the seed order
  // moves the bundled output. The investigation's premise -- "the map's key set is
  // always 0..n-1, so its traversal order is a fixed function of n, not of the
  // data" -- is true but insufficient: fixed given n does not mean ascending.
  // The author ruled to apply the fix and regenerate the affected baselines
  // (tools/baseline/results/fixed_580_summary.csv and fixed_1866_summary.csv);
  // the pre-R2 figures are recorded in SPEC-RC-011 progress.md section E.2.7 and
  // spec.md section D.7. Order invariance (REQ-011-001) is delivered by the
  // canonical input sort in runRichCluster and does not depend on this loop; it
  // holds identically with and without this change (measured both ways).
  for (int node = 0; node < n_terms; ++node) {
    // D2: filterSeeds dominates the profile (measured: n=1600 runs 13.6 s
    // end-to-end), so this is the loop a user most needs to be able to leave.
    // Per outer node, not inside filterSeed's candidate scan.
    Rcpp::checkUserInterrupt();

    // C9: copy the adjacency hash-set into an ORDERED set before it is scanned,
    // so the candidate order (and every tie in filterSeed) is the same on every
    // standard library.  The AdjacencyList itself is only ever used for
    // membership and this copy, so it can stay hashed.
    const std::unordered_set<int>& neighbors = adjList.getAdjList().at(node);
    std::set<int> neighbors_set(neighbors.begin(), neighbors.end());
    std::set<int> cluster = filterSeed(node, neighbors_set);
    clusList.addCluster(cluster);
  }
  if (verbose) Rcpp::Rcout << "Done filtering." << std::endl;
}

std::set<int> richCluster::filterSeed(
    int node, std::set<int> neighbors
) {
  std::set<int> cluster{node};
  while (true) {
    int bestN = -1;
    double bestLink = -1.0;
    
    for (int n : neighbors) {
      if (cluster.count(n)) continue;
      std::set<int> n_set{n};
      
      try {
        double link = lm.computeLinkage(cluster, n_set);
        if (link > bestLink) {
          bestLink = link;
          bestN = n;
        } 
      } catch (const std::exception& e) {
        Rcpp::Rcerr << "  EXCEPTION during linkage: " << e.what() << std::endl;
        throw; // rethrow to bubble up
      } 
    }
    // DS-11: strict everywhere.  `< cutoff` admitted a candidate sitting
    // exactly ON the cutoff, which is the `>=` convention by negation; `<=`
    // makes the seed filter agree with findBestMergePartner's `> cutoff`.
    if (bestLink <= lm.getCutoff() || bestN == -1)
      break;
    cluster.insert(bestN);
  }
  return cluster;
}


void richCluster::mergeClusters() {
  if (verbose) Rcpp::Rcout << "Starting cluster merging..." << std::endl;

  bool mergingPossible = true;
  int iteration = 0;

  while (mergingPossible) {
    // D2: one check per merge iteration.  Each iteration is O(clusters^2)
    // linkage evaluations, and the loop only ends when an entire pass merges
    // nothing, so without this the last phase of a long run was unstoppable.
    Rcpp::checkUserInterrupt();

    iteration++;
    if (verbose) Rcpp::Rcout << "Merge iteration " << iteration << "..." << std::endl;
    int nMerged = 0;
    auto& clusters = clusList.getList();
    
    for (auto it1 = clusters.begin(); it1 != clusters.end(); ++it1) {
      auto it2 = findBestMergePartner(it1, clusters);
      
      if (it2 != clusters.end() && it2 != it1) {
        // Merge cluster2 into cluster1
        it1->insert(it2->begin(), it2->end());
        clusters.erase(it2);  // immediately erase
        nMerged++;
      }
    } 
    if (verbose) Rcpp::Rcout << "  Number of merges in this iteration: " << nMerged << std::endl;
    if (nMerged == 0) {
      if (verbose) Rcpp::Rcout << "No more merges possible. Merging complete." << std::endl;
      break;
    } 
  }
  clusList.deduplicate();
}

ClusterList::ClusterIt richCluster::findBestMergePartner(
    ClusterList::ClusterIt it1, std::list<std::set<int>>& clusters
) {
  double bestLink = -1.0;
  auto bestIt = clusters.end();
   
  for (auto it2 = clusters.begin(); it2 != clusters.end(); ++it2) {
    if (it1 == it2) continue;
     
    double link = lm.computeLinkage(*it1, *it2);
    if (link > bestLink && link > lm.getCutoff()) {
      bestLink = link;
      bestIt = it2;
    } 
  }
  
  return bestIt;
} 



// the exported function to R
// [[Rcpp::export]]
Rcpp::List runRichCluster(Rcpp::CharacterVector terms,
                          Rcpp::CharacterVector geneIDs,
                          std::string distanceMetric, double distanceCutoff,
                          std::string linkageMethod, double linkageCutoff,
                          std::string geneDelim = ",",
                          bool verbose = false) {
  if (verbose) Rcpp::Rcout << "Starting richCluster..." << std::endl;
  if (verbose) Rcpp::Rcout << "terms.size = " << terms.size() << std::endl;
  if (verbose) Rcpp::Rcout << "geneIDs.size = " << geneIDs.size() << std::endl;
  try {
    if (terms.size() != geneIDs.size())
      throw std::invalid_argument("input vectors (terms, geneIDs) must be the same size");

    const int n = terms.size();

    // D3: runRichCluster is exported, so callers reach it without cluster()'s
    // validate_inputs().  Unvalidated, the cutoffs did not merely misbehave --
    // measured on two terms sharing no genes (kappa = 0), distanceCutoff = -1
    // with linkageCutoff = -1 returned the single cluster "B, A", because
    // 0 > -1 admits every pair; NA_real_, NaN and Inf each returned 2
    // singletons, because every comparison against NaN is false.  !(x > 0.0) is
    // the form the DAVID entry uses (src/DavidClustering.cpp) and is what
    // rejects NA/NaN -- a naive `x <= 0.0 || x > 1.0` passes them.  Messages
    // mirror R/cluster.R validate_inputs() so both entry points read alike.
    if (!(distanceCutoff > 0.0) || distanceCutoff > 1.0)
      throw std::invalid_argument("distance_cutoff must be between 0 and 1.");
    if (!(linkageCutoff > 0.0) || linkageCutoff > 1.0)
      throw std::invalid_argument("linkage_cutoff must be between 0 and 1.");

    // D1: an empty delimiter is an unbreakable hang, not a slow run.  Both
    // tokenisers advance with start = end + delimiter.length(), and
    // std::string::find("", start) returns start while length() is 0, so start
    // never moves (measured directly against the shipped function:
    // StringUtils::splitStringToUnorderedSet("g1,g2", "") had not returned
    // after 8 s, and the R session had to be killed with SIGKILL -- the loop
    // allocates nothing, so it spins silently forever).  Guarding here rather
    // than reshaping the tokenisers keeps the fix at the one entry point every
    // caller passes through, cluster(gene_delim=) included.
    if (geneDelim.empty())
      throw std::invalid_argument(
        "gene_delim must be a non-empty string: an empty delimiter cannot split "
        "a gene list, and the tokeniser would never terminate. Pass the "
        "separator used in your GeneID column (\",\" is the default).");

    // DS-12, second half: the widening in DistanceMatrix removes the undefined
    // behaviour, but a correct 64-bit allocation of n*n doubles still fails on
    // any real machine well before it is reachable -- and the blanket handler
    // below turns that std::length_error / bad_alloc into "Unknown C++
    // exception occurred.", which names nothing the caller can act on.  Check
    // up front and say what is actually wrong.
    {
      const double gib = (static_cast<double>(n) * static_cast<double>(n) *
                          sizeof(double)) / (1024.0 * 1024.0 * 1024.0);
      if (gib > 16.0) {
        throw std::invalid_argument(
          "too many terms: " + std::to_string(n) + " terms need a " +
          std::to_string(n) + " x " + std::to_string(n) +
          " distance matrix of about " + std::to_string(static_cast<long long>(gib)) +
          " GiB. Filter the input harder (raise min_value's stringency) or split "
          "the analysis.");
      }
    }

    // --- SPEC-RC-011 REQ-011-001/003: canonical, locale-independent ordering.
    // Key: (ascii_fold(Term), Term, GeneID), bytewise, stable.
    // ascii_fold maps ONLY bytes 'A'..'Z' to 'a'..'z'; every other byte is
    // unchanged. Do NOT use std::tolower / std::locale (locale-dependent).
    std::vector<std::string> t = Rcpp::as<std::vector<std::string>>(terms);
    std::vector<std::string> g = Rcpp::as<std::vector<std::string>>(geneIDs);
    std::vector<std::string> ft(n);
    for (int i = 0; i < n; ++i) {
      ft[i] = t[i];
      for (char& c : ft[i]) if (c >= 'A' && c <= 'Z') c = char(c + 32);
    }
    std::vector<int> perm(n);                 // perm[k] = caller index of sorted pos k
    std::iota(perm.begin(), perm.end(), 0);
    std::stable_sort(perm.begin(), perm.end(), [&](int a, int b) {
      if (ft[a] != ft[b]) return ft[a] < ft[b];
      if (t[a]  != t[b])  return t[a]  < t[b];
      return g[a] < g[b];
    });
    Rcpp::CharacterVector t_sorted(n), g_sorted(n);
    for (int k = 0; k < n; ++k) {
      t_sorted[k] = terms[perm[k]];
      g_sorted[k] = geneIDs[perm[k]];
    }

    richCluster RC(t_sorted, g_sorted,
                   distanceMetric, distanceCutoff,
                   linkageMethod, linkageCutoff,
                   geneDelim, verbose);
    RC.computeDistances();
    RC.filterSeeds();
    RC.mergeClusters();

    // --- SPEC-RC-011 REQ-011-002: restore the CALLER's index basis on export.
    // distance_matrix: dmc(perm[k], perm[l]) = dms(k, l); caller-order dimnames.
    Rcpp::NumericMatrix dms = RC.export_dm();
    Rcpp::NumericMatrix dmc(n, n);
    for (int k = 0; k < n; ++k)
      for (int l = 0; l < n; ++l)
        dmc(perm[k], perm[l]) = dms(k, l);
    dmc.attr("dimnames") = Rcpp::List::create(terms, terms);

    // all_clusters: remap TermIndices tokens k -> perm[k] ELEMENT-WISE,
    // preserving emission order.  Since C9 (converge ledger, 2026-09-04) that
    // order is canonical: clusters are std::set<int> over the sorted-basis
    // index, so members come out in canonical term order, and export_r()
    // numbers clusters by content (lexicographic over members).  This
    // element-wise remap carries that order into the caller's basis unchanged.
    // (Before C9 the emission order was the standard library's hash-set
    // iteration and "never re-sort" kept the bundled output byte-identical;
    // that invariant is superseded -- the order differed between macOS and
    // Linux, and the bundled output has been regenerated in canonical order.)
    Rcpp::DataFrame cls = RC.export_cl();
    Rcpp::CharacterVector ti = cls["TermIndices"];
    Rcpp::CharacterVector ti2(ti.size());
    for (int r = 0; r < ti.size(); ++r) {
      std::string s = Rcpp::as<std::string>(ti[r]);
      std::string out;
      size_t pos = 0;
      while (pos < s.size()) {
        size_t comma = s.find(", ", pos);
        std::string tok = (comma == std::string::npos)
                            ? s.substr(pos) : s.substr(pos, comma - pos);
        if (!out.empty()) out += ", ";
        out += std::to_string(perm[std::stoi(tok)]);
        if (comma == std::string::npos) break;
        pos = comma + 2;
      }
      ti2[r] = out;
    }
    Rcpp::DataFrame cls2 = Rcpp::DataFrame::create(
      Rcpp::Named("Cluster")     = cls["Cluster"],
      Rcpp::Named("TermNames")   = cls["TermNames"],
      Rcpp::Named("TermIndices") = ti2);

    return Rcpp::List::create(
      Rcpp::_["distance_matrix"] = dmc,
      Rcpp::_["all_clusters"]    = cls2);
  } catch (Rcpp::internal::InterruptedException&) {
    // D2: Rcpp::checkUserInterrupt() signals by throwing this, and it does NOT
    // derive from std::exception -- the catch (...) below would swallow every
    // Ctrl-C into "Unknown C++ exception occurred." and leave the run
    // uninterruptible after all.  END_RCPP in RcppExports.cpp is what must see
    // it (it answers with Rf_onintr()), so pass it through untouched.
    throw;
  } catch (const std::exception& e) {
    Rcpp::stop("C++ exception: %s", e.what());
  } catch (...) { 
    Rcpp::stop("Unknown C++ exception occurred.");
  } 
}
