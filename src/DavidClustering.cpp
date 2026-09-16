#include "DavidClustering.h"
#include "StringUtils.h"
#include <Rcpp.h>

DavidClustering::DavidClustering(
    const Rcpp::CharacterVector& terms,
    const Rcpp::CharacterVector& geneIDs,
    double similarityThreshold,
    int initialGroupMembership,
    int finalGroupMembership,
    double multipleLinkageThreshold,
    bool verbose
) : terms(terms),
    geneIDs(geneIDs),
    similarityThreshold(similarityThreshold),
    initialGroupMembership(initialGroupMembership),
    finalGroupMembership(finalGroupMembership),
    multipleLinkageThreshold(multipleLinkageThreshold),
    verbose(verbose) {

    // SPEC-RC-009 / DS-04 -- reject out-of-domain parameters loudly.  The
    // size_t cast at findInitialSeeds() (initialGroupMembership - 1) wraps
    // silently for non-positive values, and no kappa can exceed 1, so these
    // previously produced zero clusters with no diagnostic.  Messages mirror
    // the R-side validate_david_inputs() exactly.  !(x > 0.0) also rejects NaN.
    if (!(similarityThreshold > 0.0) || similarityThreshold > 1.0) {
        Rcpp::stop("similarity_threshold must be between 0 and 1.");
    }
    if (!(multipleLinkageThreshold > 0.0) || multipleLinkageThreshold > 1.0) {
        Rcpp::stop("multiple_linkage_threshold must be between 0 and 1.");
    }
    if (initialGroupMembership < 1) {
        Rcpp::stop("initial_group_membership must be a whole number >= 1.");
    }
    if (finalGroupMembership < 1) {
        Rcpp::stop("final_group_membership must be a whole number >= 1.");
    }

    n_terms = terms.size();
    // Convert the incoming CharacterVector of comma-separated gene IDs to a
    // standard vector so StringUtils utilities can operate on it.
    std::vector<std::string> geneIDsVector = Rcpp::as<std::vector<std::string>>(geneIDs);
    totalGeneCount = StringUtils::countUniqueElements(geneIDsVector);
    kappaMatrix.resize(n_terms, std::vector<double>(n_terms, 0.0));
}

Rcpp::List DavidClustering::run() {
    if (verbose) Rcpp::Rcout << "Calculating kappa scores..." << std::endl;
    calculateKappaScores();

    if (verbose) Rcpp::Rcout << "Finding initial seeds..." << std::endl;
    findInitialSeeds();

    if (verbose) Rcpp::Rcout << "Merging seeds..." << std::endl;
    mergeSeeds();

    // Format final clusters for output
    std::vector<int> clusterColumn;
    std::vector<std::string> termNamesColumn;
    std::vector<std::string> termIndicesColumn;

    int clusterNum = 1;
    for (const auto& cluster : finalClusters) {
        if (cluster.size() >= static_cast<size_t>(finalGroupMembership)) {
            std::vector<std::string> clusterTermNames;
            std::string termIndicesStr;

            for (int term_index : cluster) {
                clusterTermNames.push_back(Rcpp::as<std::string>(terms[term_index]));
                if (!termIndicesStr.empty()) {
                    termIndicesStr += ", ";
                }
                termIndicesStr += std::to_string(term_index);
            }

            clusterColumn.push_back(clusterNum++);
            termNamesColumn.push_back(StringUtils::vectorToString(clusterTermNames, ", "));
            termIndicesColumn.push_back(termIndicesStr);
        }
    }

    Rcpp::DataFrame clusters_df = Rcpp::DataFrame::create(
        Rcpp::Named("Cluster") = clusterColumn,
        Rcpp::Named("TermNames") = termNamesColumn,
        Rcpp::Named("TermIndices") = termIndicesColumn
    );

    return Rcpp::List::create(
        Rcpp::Named("clusters") = clusters_df
    );
}

void DavidClustering::calculateKappaScores() {
    for (int i = 0; i < n_terms; ++i) {
        // SPEC-RC-007 / T1-16, OQ-2 Reading B (author ruling, 2026-08-25):
        // term1_genes is invariant in j, so it is split ONCE per outer
        // iteration rather than once per pair -- O(n) splits instead of
        // O(n^2).  This is exactly the hoisting src/RichCluster.cpp already
        // does in computeDistances(); this makes DavidClustering match it.
        //
        // A pure refactor: the value bound on each pass is byte-for-byte the
        // value the inner-loop split produced, so the proof obligation is a
        // before/after bit-identity across all six artifact shapes, not a new
        // assertion.  There is no RED step for a hoist.
        std::unordered_set<std::string> term1_genes = StringUtils::splitStringToUnorderedSet(Rcpp::as<std::string>(geneIDs[i]), ",");

        for (int j = i + 1; j < n_terms; ++j) {
            std::unordered_set<std::string> term2_genes = StringUtils::splitStringToUnorderedSet(Rcpp::as<std::string>(geneIDs[j]), ",");

            // DS-01: a term with no genes shares nothing with anything.
            // Without this guard two empty gene sets fall into the aab == 1
            // branch below and score kappa = 1.0 (and, when every term is
            // empty, totalGeneCount is 0 and the divisions below are 0/0).
            if (term1_genes.empty() && term2_genes.empty()) {
                kappaMatrix[i][j] = 0.0;
                kappaMatrix[j][i] = 0.0;
                continue;
            }

            int term1term2 = 0;
            for (const auto& gene : term1_genes) {
                if (term2_genes.count(gene)) {
                    term1term2++;
                }
            }

            int posTerm1Total = term1_genes.size();
            int posTerm2Total = term2_genes.size();
            int term1only = posTerm1Total - term1term2;
            int term2only = posTerm2Total - term1term2;
            int term1term2Non = totalGeneCount - term1term2 - term1only - term2only;

            double oab = static_cast<double>(term1term2 + term1term2Non) / totalGeneCount;
            double aab = (static_cast<double>(posTerm1Total) * posTerm2Total + static_cast<double>(totalGeneCount - posTerm1Total) * (totalGeneCount - posTerm2Total)) / (static_cast<double>(totalGeneCount) * totalGeneCount);

            // N3: THE TWO KAPPA IMPLEMENTATIONS DIVERGE HERE, DELIBERATELY.
            //
            // This is the second of two kappa implementations. The other is
            // DistanceMetric::getKappa (src/DistanceMetric.cpp), which backs
            // cluster(); this one backs david_cluster() and is exported to R
            // through runDavidClusteringWithKappa(). The core algebra is
            // identical term for term (verified bit-identical), and since N2
            // the degenerate branch agrees too: aab == 1 is perfect agreement,
            // and BOTH now return 1.0.
            //
            // ONE divergence remains, and it is intentional: getKappa FLOORS a
            // negative kappa at 0 (its OD-5 ruling, so cluster() exports a
            // [0, 1] scale), while this path carries the raw value. Do NOT add
            // the floor here -- DAVID's algorithm compares raw kappa against
            // similarityThreshold and multipleLinkageThreshold, and flooring
            // would move those decisions. Measured on 8 random terms over 30
            // genes, this matrix spans [-0.4884, 0.3636] with 17 of 28 pairs
            // negative; cluster()'s spans [0, 1]. Both are user-visible.
            double kappa = (aab == 1) ? 1.0 : (oab - aab) / (1 - aab);

            kappaMatrix[i][j] = kappa;
            kappaMatrix[j][i] = kappa;
        }
    }
}

void DavidClustering::findInitialSeeds() {
    for (int i = 0; i < n_terms; ++i) {
        TermSet neighbors;
        for (int j = 0; j < n_terms; ++j) {
            if (i == j) continue;
            if (kappaMatrix[i][j] > similarityThreshold) {
                neighbors.insert(j);
            }
        }

        if (neighbors.size() >= static_cast<size_t>(initialGroupMembership - 1)) {
            TermSet current_seed = neighbors;
            current_seed.insert(i);

            int totalPairs = 0;
            int passedPair = 0;
            std::vector<int> seed_vec(current_seed.begin(), current_seed.end());
            for (size_t k = 0; k < seed_vec.size(); ++k) {
                for (size_t l = k + 1; l < seed_vec.size(); ++l) {
                    totalPairs++;
                    if (kappaMatrix[seed_vec[k]][seed_vec[l]] > similarityThreshold) {
                        passedPair++;
                    }
                }
            }

            if (totalPairs > 0 && (static_cast<double>(passedPair) / totalPairs) > multipleLinkageThreshold) {
                initialSeeds.push_back(current_seed);
            }
        }
    }
}

void DavidClustering::mergeSeeds() {
    std::list<TermSet> working(initialSeeds.begin(), initialSeeds.end());

    while (true) {
        std::list<TermSet> remainingSeeds = std::move(working);
        const size_t sizeBefore = remainingSeeds.size();
        std::list<TermSet> merged;

        while (!remainingSeeds.empty()) {
            TermSet currentSeed = remainingSeeds.front();
            remainingSeeds.pop_front();

            while (true) {
                double bestScore = 0.0;
                auto bestIt = remainingSeeds.end();

                for (auto it = remainingSeeds.begin(); it != remainingSeeds.end(); ++it) {
                    double score = calculateDiceCoefficient(currentSeed, *it);
                    if (score > multipleLinkageThreshold && score > bestScore) {
                        bestScore = score;
                        bestIt = it;
                    }
                }

                if (bestIt != remainingSeeds.end()) {
                    currentSeed.insert(bestIt->begin(), bestIt->end());
                    remainingSeeds.erase(bestIt);
                } else {
                    break;
                }
            }
            merged.push_back(currentSeed);
        }

        const bool converged = (merged.size() == sizeBefore);
        working = std::move(merged);
        if (converged) {
            break;
        }
    }

    finalClusters.assign(working.begin(), working.end());
}

double DavidClustering::calculateDiceCoefficient(const TermSet& seed1, const TermSet& seed2) {
    int commonCount = 0;
    for (int term : seed1) {
        if (seed2.count(term)) {
            commonCount++;
        }
    }
    return 2.0 * commonCount / (seed1.size() + seed2.size());
}

// Exported function to be called from R
// [[Rcpp::export]]
Rcpp::List runDavidClustering(
    Rcpp::CharacterVector terms,
    Rcpp::CharacterVector geneIDs,
    double similarityThreshold,
    int initialGroupMembership,
    int finalGroupMembership,
    double multipleLinkageThreshold,
    bool verbose = false) {

    DavidClustering david(
        terms,
        geneIDs,
        similarityThreshold,
        initialGroupMembership,
        finalGroupMembership,
        multipleLinkageThreshold,
        verbose
    );

    return david.run();
}

// ---------------------------------------------------------------------------
// T3-04a -- ADDITIVE, OUTPUT-NEUTRAL EXPOSURE OF THE DAVID KAPPA MATRIX
//
// runDavidClustering() above is UNTOUCHED, and so are run() and
// calculateKappaScores().  This is a SECOND, additional entry point: it runs
// exactly the same pipeline and additionally returns the kappa matrix that
// pipeline used, appended AFTER every element run() already returns.  Nothing
// existing is removed, renamed or reordered, and no numeric code moved --
// which is what lets T3-04a V1 prove output neutrality by bit-identity against
// the v102-david and v102-dc tags T3-03 froze on the stock build.
//
// The matrix is read back through DavidClustering::getKappaMatrix() AFTER
// run() has completed, so it is the very matrix the clustering consumed, not a
// separately computed copy (T3-04a V3).
// ---------------------------------------------------------------------------

namespace {

// vector<vector<double>> -> NumericMatrix.  No dimnames: the exposed object is
// the raw n_terms x n_terms matrix, indexed by the same term order as the
// `terms` argument.
Rcpp::NumericMatrix kappaMatrixToNumericMatrix(const std::vector<std::vector<double>>& m) {
    int n = static_cast<int>(m.size());
    Rcpp::NumericMatrix out(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            out(i, j) = m[i][j];
        }
    }
    return out;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List runDavidClusteringWithKappa(
    Rcpp::CharacterVector terms,
    Rcpp::CharacterVector geneIDs,
    double similarityThreshold,
    int initialGroupMembership,
    int finalGroupMembership,
    double multipleLinkageThreshold,
    bool verbose = false) {

    DavidClustering david(
        terms,
        geneIDs,
        similarityThreshold,
        initialGroupMembership,
        finalGroupMembership,
        multipleLinkageThreshold,
        verbose
    );

    Rcpp::List result = david.run();
    result.push_back(kappaMatrixToNumericMatrix(david.getKappaMatrix()), "kappa_matrix");
    return result;
}
