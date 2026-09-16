#ifndef DavidClustering_h
#define DavidClustering_h

#include <Rcpp.h>
#include <string>
#include <vector>
#include <unordered_set>

class DavidClustering {
public:
    DavidClustering(
        const Rcpp::CharacterVector& terms,
        const Rcpp::CharacterVector& geneIDs,
        double similarityThreshold,
        int initialGroupMembership,
        int finalGroupMembership,
        double multipleLinkageThreshold,
        bool verbose = false
    );

    Rcpp::List run();

    // --- T3-04a: additive, output-neutral observability -------------------
    // Read-only access to the kappa matrix this instance computed.
    //
    // WHY THIS EXISTS.  calculateKappaScores() already counts the gene
    // intersection correctly with a manual membership loop, so it is the
    // in-repo reference implementation the Wave-0 oracle gate validates
    // against.  As shipped, `kappaMatrix` and `calculateKappaScores()` were
    // both private and run() returned only `clusters`, so that matrix was not
    // observable from R by any means and the gate was not executable at all.
    //
    // WHAT THIS IS NOT.  This accessor computes nothing and changes nothing:
    // it returns a const reference to the member run() already filled.  The
    // computation is untouched; only its observability is added.  Call it
    // after run() (or after calculateKappaScores()); before either, the matrix
    // is still the constructor's 0.0 fill.
    //
    // NOTE the diagonal is 0, not 1: calculateKappaScores() loops
    // `for (j = i + 1; ...)` and never writes it.  That is characterised, not
    // fixed here -- the diagonal is never read by the clustering, and this
    // item is additive and output-neutral by construction.
    const std::vector<std::vector<double>>& getKappaMatrix() const { return kappaMatrix; }

private:
    using TermSet = std::unordered_set<int>;

    void calculateKappaScores();
    void findInitialSeeds();
    void mergeSeeds();
    double calculateDiceCoefficient(const TermSet& seed1, const TermSet& seed2);

    // Input data
    Rcpp::CharacterVector terms;
    Rcpp::CharacterVector geneIDs;
    int n_terms;
    int totalGeneCount;

    // Parameters
    double similarityThreshold;
    int initialGroupMembership;
    int finalGroupMembership;
    double multipleLinkageThreshold;
    bool verbose;   // false => run() narrates nothing (CRAN default)

    // Internal data structures
    std::vector<std::vector<double>> kappaMatrix;
    std::vector<TermSet> initialSeeds;
    std::vector<TermSet> finalClusters;
};

#endif // DavidClustering_h
