//
//  DistanceMetric.cpp
//  richCluster
//
//  Created by Sarah on 6/2/25.
//  edited by Junguk Hur on 8/22/2025
//

#include <stdio.h>
#include "DistanceMetric.h"
#include <unordered_set>
#include <string>
#include <stdexcept>

double DistanceMetric::computeDistance(const std::unordered_set<std::string>& t1_genes,
                                       const std::unordered_set<std::string>& t2_genes,
                                       int totalGeneCount) {
  if (metric=="kappa")
    return getKappa(t1_genes, t2_genes, totalGeneCount);
  else if (metric=="jaccard")
    return getJaccard(t1_genes, t2_genes);
  else if (metric=="dice")
    return getDice(t1_genes, t2_genes);
  else
    throw std::invalid_argument("unsupported distance metric: " + metric);
}


// the various distance metric computations
// kappa is the standard
double DistanceMetric::getKappa(const std::unordered_set<std::string>& t1_genes,
                                const std::unordered_set<std::string>& t2_genes,
                                int totalGeneCount) {
  
  // Count the genes shared by t1_genes and t2_genes.
  // NOTE: std::set_intersection requires SORTED ranges; an unordered_set is a hash
  // table and is not sorted, so it must not be used here. Count by membership test
  // instead, as getJaccard below already does.
  double common = 0;
  for (const auto& gene : t1_genes) {
    if (t2_genes.count(gene)) {
      ++common;
    }
  }

  // No shared genes -> no agreement to score. This branch also swallows the
  // genuinely-empty case (both terms gene-less), which is why the
  // chance_agree == 1 guard below can only ever be reached by the OTHER
  // degenerate solution, perfect agreement -- see the N2 note there.
  if (common == 0) {
    return 0.0; // return 0 if no overlapping genes
  } 
  
  double t1_only = t1_genes.size() - common; // Genes unique to t1_genes
  double t2_only = t2_genes.size() - common; // Genes unique to t2_genes
  
  double unique = totalGeneCount - common - t1_only - t2_only; // Count of all genes not found in either term
  
  double relative_observed_agree = (common + unique) / totalGeneCount;
  double chance_yes = ((common + t1_only) / totalGeneCount) * ((common + t2_only) / totalGeneCount);
  double chance_no = ((unique + t1_only) / totalGeneCount) * ((unique + t2_only) / totalGeneCount);
  double chance_agree = chance_yes + chance_no;
  
  // N2: chance_agree == 1 is PERFECT AGREEMENT, not a nothing-in-common case.
  // With p1 = |t1|/N and p2 = |t2|/N, chance_agree = p1*p2 + (1-p1)*(1-p2)
  // reaches 1 in exactly two ways:
  //
  //   p1 = p2 = 0  both terms empty.  UNREACHABLE from here -- `common == 0`
  //                returns above, so this branch never sees it.
  //   p1 = p2 = 1  both terms hold every gene in the universe, i.e. the two
  //                terms are IDENTICAL.
  //
  // So the only branch that arrives here is the identical one, where
  // relative_observed_agree is 1 as well: kappa is 0/0, and the limit along
  // perfect agreement is 1, the MAXIMUM of the scale. Returning 0 -- the
  // minimum -- put two identical terms in separate singleton clusters while
  // jaccard and dice both scored them 1.0. Verified: two identical 12-gene
  // terms in a 12-gene universe (tools/cpp-units case K2).
  //
  // The same value is correct for a near-degenerate chance_agree that rounds to
  // exactly 1: relative_observed_agree is then ~1 too, so 1 remains the limit.
  if (chance_agree == 1)
    return 1.0;

  double kappa = (relative_observed_agree - chance_agree) / (1 - chance_agree);

  // OD-5: floor negative kappa at 0, so the exported scale is [0, 1].
  // DAVID documents kappa as [0, 1] and disqualifies low-overlap pairs rather than
  // carrying them as negatives; a negative kappa has no fixed lower bound and no
  // meaningful interpretation as term similarity. This also makes the `common == 0`
  // early return above consistent with the general case rather than a special case.
  //
  // N3: this floor is DELIBERATE and SCOPED TO cluster(). The second kappa
  // implementation -- DavidClustering::calculateKappaScores(), which backs
  // david_cluster() and exports its matrix through runDavidClusteringWithKappa
  // -- does NOT floor, and must not be changed to: DAVID's algorithm consumes
  // raw kappa against its own thresholds. The two therefore disagree on sign
  // by design, and only here: the algebra is otherwise identical term for term,
  // and the chance_agree == 1 branch above now agrees with it at 1.0. Measured
  // on 8 random terms over 30 genes, DAVID's exported matrix spans
  // [-0.4884, 0.3636] with 17 of 28 pairs negative; cluster()'s spans [0, 1].
  // See the mirror note at src/DavidClustering.cpp (calculateKappaScores).
  return kappa < 0 ? 0.0 : kappa;
}

double DistanceMetric::getJaccard(const std::unordered_set<std::string>& t1_genes,
                                  const std::unordered_set<std::string>& t2_genes) {
  double common = 0;
  for (const auto& gene : t1_genes) {
    if (t2_genes.count(gene)) {
      ++common;
    }
  }
  // Jaccard index = |A n B| / |A u B|, and |A u B| = |A| + |B| - |A n B|.
  // Dividing by |A| + |B| instead yields Dice/2, which is capped at 0.5 and can
  // therefore never clear a distance_cutoff of 0.5.
  double union_size = static_cast<double>(t1_genes.size()) + static_cast<double>(t2_genes.size()) - common;

  if (union_size == 0)
    return 0.0; // both terms carry no genes

  return common / union_size;
}

double DistanceMetric::getDice(const std::unordered_set<std::string>& t1_genes,
                               const std::unordered_set<std::string>& t2_genes) {
  // Count shared genes by membership test.
  // NOTE: std::set_intersection requires SORTED ranges; an unordered_set is a
  // hash table and is not sorted, so it must not be used here. This mirrors
  // getKappa and getJaccard above.
  double common = 0;
  for (const auto& gene : t1_genes) {
    if (t2_genes.count(gene)) {
      ++common;
    }
  }

  // Dice coefficient = 2|A n B| / (|A| + |B|).
  // The factor of 2 in the numerator is load-bearing: common / (|A| + |B|)
  // yields Dice/2, which was the pre-existing jaccard defect.
  double denom = static_cast<double>(t1_genes.size()) + static_cast<double>(t2_genes.size());

  if (denom == 0)
    return 0.0; // both terms carry no genes

  return (2.0 * common) / denom;
}
