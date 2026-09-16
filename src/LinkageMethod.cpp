//
//  LinkageMethod.cpp
//  richCluster
//
//  Created by Sarah on 6/2/25.
//

#include <stdio.h>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include "LinkageMethod.h"

using Cluster = LinkageMethod::Cluster;

// Returned when two clusters share every term, so no cross-pair remains to measure.
// -Inf can never exceed a cutoff, so such a pair is never merged.
static const double NO_LINK = -std::numeric_limits<double>::infinity();

// NOTE ON DIRECTION: distFct returns a SIMILARITY (higher = more alike), not a
// distance. So single linkage takes the MAXIMUM similarity (the closest pair) and
// complete linkage the MINIMUM (the furthest pair) -- the opposite of what the
// distance-based names suggest. This yields single >= average >= complete.

double LinkageMethod::computeLinkage(const Cluster& cluster1, const Cluster& cluster2) {
  if (method=="single")
    return LinkageMethod::single(cluster1, cluster2);
  else if (method=="complete")
    return LinkageMethod::complete(cluster1, cluster2);
  else if (method=="average")
    return LinkageMethod::average(cluster1, cluster2);
  else if (method=="ward")
    return LinkageMethod::ward(cluster1, cluster2);
  else
    // DS-07: mirror DistanceMetric::computeDistance -- an unknown method is an
    // error, never a silent linkage of 0 (which returns wrong clusters).
    throw std::invalid_argument("unsupported linkage method: " + method);
}

double LinkageMethod::single(const Cluster& cluster1, const Cluster& cluster2) {
  double maxSim = NO_LINK; // seed below every attainable similarity
  int nCompared = 0;
  for (auto i = cluster1.begin(); i!= cluster1.end(); ++i) {
    for (auto j = cluster2.begin(); j!= cluster2.end(); ++j) {
      // Skip a term the two clusters share: its self-distance is a sentinel, not a
      // similarity. Compare the TERM INDICES -- iterators into two different
      // containers never compare equal, so `i==j` here never fires.
      if (*i==*j)
        continue;
      // de-reference the two term iterators
      int t1 = *i;
      int t2 = *j;
      // get similarity between them
      double sim = distFct(t1, t2);
      if (sim > maxSim)
        maxSim = sim;
      nCompared++;
    }
  }
  if (nCompared == 0)
    return NO_LINK;
  return maxSim;
}

double LinkageMethod::complete(const Cluster& cluster1, const Cluster& cluster2) {
  // seed above every attainable similarity, so a negative similarity is not
  // silently floored at 0 the way a seed of 0 would floor it
  double minSim = std::numeric_limits<double>::infinity();
  int nCompared = 0;
  for (auto i = cluster1.begin(); i!= cluster1.end(); ++i) {
    for (auto j = cluster2.begin(); j!= cluster2.end(); ++j) {
      // see the note in single(): compare term indices, not iterators
      if (*i==*j)
        continue;
      // de-reference the two term iterators
      int t1 = *i;
      int t2 = *j;
      // get similarity between them
      double sim = distFct(t1, t2);
      if (sim < minSim)
        minSim = sim;
      nCompared++;
    }
  }
  if (nCompared == 0)
    return NO_LINK;
  return minSim;
}

double LinkageMethod::average(const Cluster& cluster1, const Cluster& cluster2) {
  double totalDist = 0;
  int n_terms = 0;
  for (auto i = cluster1.begin(); i!= cluster1.end(); ++i) {
    for (auto j = cluster2.begin(); j!= cluster2.end(); ++j) {
      // see the note in single(): compare term indices, not iterators
      if (*i==*j)
        continue;
      // de-reference the two term iterators
      int t1 = *i;
      int t2 = *j;
      // get distance between them
      double dist = distFct(t1, t2);
      totalDist += dist;
      n_terms ++;
    }
  }
  if (n_terms == 0)
    return NO_LINK; // would otherwise divide by zero and yield NaN
  return totalDist/n_terms;
}


// Ward's minimum-variance criterion, expressed on the selection scale.
//
// DIRECTION.  Ward yields a merge COST (lower = better merge), but
// findBestMergePartner (RichCluster.cpp) selects the HIGHEST value and gates it
// on a cutoff in (0, 1].  The cost is mapped through
//
//     link = 1 - 2 * dESS
//
// (the pre-F5 pair was d = 1 - s with link = 1 - sqrt(2 * dESS); F5 moved both
// together, see the tail of this function.)  It is strictly decreasing in the
// cost.  The transform is DERIVED, not chosen: for two singleton clusters
// dESS = d^2/2 = (1 - s)/2, so the expression collapses to s, the raw
// similarity.  That anchor is what makes linkage_cutoff mean approximately the
// same thing here as for single/complete/average.  Do NOT substitute
// 1/(1+dESS): it merges every singleton pair at the default cutoff, including
// dissimilar and anti-similar ones.  See SPEC-RC-003 s3.
//
// The value may be negative for costly merges; that is safe, because the gate
// is link > cutoff with cutoff > 0, so such a pair is rejected on the cutoff
// condition before the bestLink seed is ever consulted.  It may also EXCEED 1,
// for the nested-cluster case the N1 block below describes.
//
// The matrix stores a SIMILARITY, and the SQUARED distance is d^2 = 1 - s.  A
// negative similarity is deliberately NOT floored: it yields d^2 > 1, which is
// correct.  (wardDistFrom does clamp d^2 at 0, which only bites for s > 1 --
// unreachable from the three metrics, all of which cap at 1.)
// F5: Ward's criterion is defined on SQUARED EUCLIDEAN distances -- the Huygens
// identity used for sepSq below is exact only when the points embed in L2.
//
// The similarity-to-distance transform decides whether that holds.  Measured on
// the bundled data (376 terms, Gower test: D is Euclidean iff -0.5 * J D^2 J is
// positive semi-definite):
//
//   d = 1 - s        min eigenvalue  kappa -3.69e-01, jaccard -7.19e-02,
//                    dice -2.69e-01                        -> NOT Euclidean
//   d = sqrt(1 - s)  min eigenvalue  all three ~ -1e-15     -> Euclidean
//
// So (1 - s) is the SQUARED distance, not the distance.  That is what makes
// dESS below a genuine minimum-variance increment rather than a Ward-shaped
// number.  It is the standard transform for these coefficients (Gower &
// Legendre 1986, on which similarity coefficients yield Euclidean distances);
// the measurement above confirms it for kappa on this data too.
//
// THE SQUARE ROOT ITSELF IS NOT LOAD-BEARING -- the earlier claim that it was
// is withdrawn.  Every caller squares the return value straight back
// (`const double d = wardDist(...); ... d * d`), so this function computes
// max(1 - s, 0) the long way round.  It is kept ONLY because removing it is
// not bit-identical: sqrt() is correctly rounded, but squaring the result is
// not exact, so (sqrt(1 - s))^2 differs from (1 - s) by up to one ulp.
// Measured over 2,000,000 draws on [0, 1.5], 49.8% of values differ, worst
// case 2.22e-16; the singleton anchor at s = 0.5 returns 0.49999999999999989
// with the round trip and exactly 0.5 without it.  Dropping it is therefore a
// strictly-more-accurate but RESULTS-CHANGING edit, deliberately not bundled
// with the N1 correction below so that fix's numeric delta stays attributable.
//
// Callers never reach the same-term sentinel: every loop skips *i == *j.
static inline double wardDistFrom(double similarity) {
  const double dSq = 1.0 - similarity;
  return dSq > 0.0 ? std::sqrt(dSq) : 0.0;
}

double LinkageMethod::ward(const Cluster& cluster1, const Cluster& cluster2) {
  auto wardDist = [this](int a, int b) { return wardDistFrom(distFct(a, b)); };

  const double nA = static_cast<double>(cluster1.size());
  const double nB = static_cast<double>(cluster2.size());

  // Cross term: sum of squared distances over A x B.
  // A term shared by both clusters is skipped -- its true self-distance is 0,
  // so a skipped pair contributes exactly what the identity says it should.
  // The denominators below keep the FULL cardinalities regardless.
  double crossSumSq = 0.0;
  int nCompared = 0;
  for (auto i = cluster1.begin(); i != cluster1.end(); ++i) {
    for (auto j = cluster2.begin(); j != cluster2.end(); ++j) {
      if (*i == *j)
        continue;
      const double d = wardDist(*i, *j);
      crossSumSq += d * d;
      nCompared++;
    }
  }
  if (nCompared == 0)
    return NO_LINK; // same guard as single/complete/average

  // Within-cluster terms, over ORDERED pairs: each unordered pair is counted
  // twice, which the 2*n^2 denominators account for.  Do not halve these.
  // The diagonal contributes d(x,x)^2 = 0, so skipping it changes nothing.
  double withinASumSq = 0.0;
  for (auto i = cluster1.begin(); i != cluster1.end(); ++i) {
    for (auto j = cluster1.begin(); j != cluster1.end(); ++j) {
      if (*i == *j)
        continue;
      const double d = wardDist(*i, *j);
      withinASumSq += d * d;
    }
  }

  double withinBSumSq = 0.0;
  for (auto i = cluster2.begin(); i != cluster2.end(); ++i) {
    for (auto j = cluster2.begin(); j != cluster2.end(); ++j) {
      if (*i == *j)
        continue;
      const double d = wardDist(*i, *j);
      withinBSumSq += d * d;
    }
  }

  // Squared centroid separation, from pairwise squared distances alone.
  double sepSq = crossSumSq   / (nA * nB)
               - withinASumSq / (2.0 * nA * nA)
               - withinBSumSq / (2.0 * nB * nB);

  // With the Euclidean transform above, sepSq is non-negative in exact
  // arithmetic, so this guard catches floating-point rounding only -- not a
  // geometry violation.  It stays because it is load-bearing either way: an
  // unclamped negative would drag the disjoint closed form below negative and
  // reverse the merge ordering for that pair.  Measured on ~9,000 real cluster
  // pairs across all three metrics, it fired 0 times.
  if (sepSq < 0.0)
    sepSq = 0.0;

  // N1 -- THE MERGED CLUSTER IS THE *SET* UNION, NOT THE MULTISET UNION.
  //
  // richCluster's clusters are NOT disjoint.  mergeClusters() (RichCluster.cpp)
  // merges with `it1->insert(it2->begin(), it2->end())` -- set insertion -- and
  // filterSeeds() emits one seed PER NODE, so mutual neighbours produce seeds
  // that share terms.  Measured on the bundled data at min_value = 1e-4, 28 of
  // 4465 final ward cluster pairs share at least one term.
  //
  // The Lance-Williams closed form (nA*nB/(nA+nB)) * sepSq is the ESS increment
  // for the MULTISET union, in which a shared term is counted TWICE, and is
  // exact ONLY for disjoint A and B.  Applied to an overlapping pair it reads
  // nA + nB for a union that actually holds nA + nB - shared terms (nA = 6,
  // nB = 13 sharing 4 gives a 15-term union scored as 19).  On the 14
  // overlapping final-cluster pairs of a real run it flipped 10 of 14 merge
  // decisions; a synthetic probe showed 0% error at shared = 0, 75% at
  // shared = 1, and 180% -- with a sign change -- at shared = 2.  The control:
  // on 200 disjoint pairs from the same run, closed form and definition agreed
  // to 3.55e-15, so the divergence is exactly the overlap term.
  //
  // The increment is therefore taken from its DEFINITION whenever the clusters
  // overlap.  With the same Huygens form used above, for any set S
  //
  //     ESS(S) = ( sum over ORDERED pairs i != j in S of d^2(i,j) ) / (2|S|)
  //     dESS   = ESS(A u B) - ESS(A) - ESS(B)
  //
  // which is exact for overlapping and disjoint pairs alike.
  //
  // WHY THE DISJOINT BRANCH SURVIVES.  The two agree exactly when A and B are
  // disjoint: then |A u B| = nA + nB, and the ordered-pair sum over the union
  // is withinA + withinB + 2*cross, so
  //
  //   dESS = (withinA + withinB + 2*cross) / (2*(nA + nB))
  //          - withinA/(2*nA) - withinB/(2*nB)
  //        = cross/(nA + nB)
  //          - withinA*nB / (2*nA*(nA + nB))          [ 1/(nA+nB) - 1/nA ]
  //          - withinB*nA / (2*nB*(nA + nB))          [ 1/(nA+nB) - 1/nB ]
  //        = (nA*nB/(nA + nB)) * sepSq
  //
  // term for term.  Keeping that branch is not cosmetic: it is the path ~99%
  // of real pairs take, and it avoids a second O(|A u B|^2) pass.  The two
  // branches are asserted to agree numerically by tools/cpp-units case W9,
  // over 500 random pairs spanning both.
  std::size_t nShared = 0;
  for (int t : cluster2) {
    if (cluster1.count(t))
      ++nShared;
  }

  double dESS;
  if (nShared == 0) {
    // Disjoint -- the closed form is exact here, and cheaper.
    dESS = (nA * nB / (nA + nB)) * sepSq;
  } else {
    // Overlapping -- build the union the merge would actually produce and take
    // the increment from the definition.
    //
    // dESS may come out NEGATIVE on this branch, and is deliberately NOT
    // clamped.  For nested clusters (A a subset of B) the union IS B, so
    // dESS = ESS(B) - ESS(A) - ESS(B) = -ESS(A) < 0 and the link exceeds 1,
    // which always clears the cutoff.  That is the correct reading: merging a
    // cluster into a superset of itself costs nothing and deduplicates.
    Cluster merged = cluster1;
    merged.insert(cluster2.begin(), cluster2.end());
    const double nU = static_cast<double>(merged.size());

    // Same ordered-pair convention as withinASumSq / withinBSumSq above.
    double withinUSumSq = 0.0;
    for (auto i = merged.begin(); i != merged.end(); ++i) {
      for (auto j = merged.begin(); j != merged.end(); ++j) {
        if (*i == *j)
          continue;
        const double d = wardDist(*i, *j);
        withinUSumSq += d * d;
      }
    }

    dESS = withinUSumSq / (2.0 * nU)
         - withinASumSq / (2.0 * nA)
         - withinBSumSq / (2.0 * nB);
  }

  // The rescale is the MATCHED PARTNER of the distance transform above; changing
  // one without the other decalibrates the cutoff.  RC-003 chose the rescale so
  // that a singleton pair returns the similarity itself, which is what lets
  // linkage_cutoff mean the same thing for ward as for the other three methods.
  //
  // Under d^2 = (1 - s), dESS is LINEAR in (1 - s), not quadratic, so the square
  // root that the old d = (1 - s) form needed is gone:
  //   singletons -> sepSq = 1 - s, dESS = (1 - s)/2, 1 - 2*dESS = s.
  // (The old pair was d = 1 - s with 1 - sqrt(2*dESS), which gave the same s.)
  return 1.0 - 2.0 * dESS;
}
