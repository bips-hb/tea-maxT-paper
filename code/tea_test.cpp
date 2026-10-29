// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(BH)]]

// Closed permutation tests for binary endpoints.
#include <Rcpp.h>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace tea_detail {

// Round statistics to a common grid when identifying ties.
using Key = std::int64_t;
using ScoreMatrix = boost::multi_array<Key, 2>;
using ScoreCube = boost::multi_array<Key, 3>;
using PairCounts = boost::numeric::ublas::triangular_matrix<
    int, boost::numeric::ublas::strict_upper>;
constexpr double scale = 1e10;
constexpr Key missing = std::numeric_limits<Key>::lowest();

// Generates reproducible pseudo-random numbers from a seed.
// Provides random indices and uniform values for the test.
class Generator {
    std::mt19937 engine;
public:
    Generator(int seed, std::uint32_t domain) {
        std::seed_seq state{static_cast<std::uint32_t>(seed), domain};
        engine.seed(state);
    }

    std::uint32_t below(std::uint32_t bound) {
        const std::uint32_t threshold = (std::uint32_t(0) - bound) % bound;
        std::uint32_t value;
        do { value = engine(); } while (value < threshold);
        return value % bound;
    }

    double uniform() {
        double value;
        do {
            const std::uint32_t a = engine() >> 5;
            const std::uint32_t b = engine() >> 6;
            value = (a * 67108864.0 + b) / 9007199254740992.0;
        } while (value == 0.0);
        return value;
    }
};

// Use the observed assignment as row 0. Then sample with replacement
// or enumerate every other valid assignment exactly once.
class Permutations {
    int n, n1;
    bool exact, exhausted = false;
    std::vector<int> observed, enumeration;
    Generator rng;
public:
    int total;
    std::vector<int> labels;

    Permutations(const Rcpp::IntegerVector& group, int draws, bool enumerate,
                 int seed) : n(group.size()), n1(Rcpp::sum(group)),
        exact(enumerate), observed(group.begin(), group.end()),
        rng(seed, 0x5045524dU), total(draws + 1), labels(observed) {
        if (!exact) return;
        long double number = 1;
        const int k = std::min(n1, n - n1);
        for (int i = 1; i <= k; ++i) {
            number *= static_cast<long double>(n - k + i) / i;
            if (number > std::numeric_limits<int>::max() - 1.0L)
                Rcpp::stop("Too many assignments for exact enumeration.");
        }
        total = static_cast<int>(std::round(number));
        enumeration.assign(n, 0);
        std::fill_n(enumeration.begin(), n1, 1);
    }

    void next(int row) {
        if (row == 0) return;
        if (!exact) {
            labels = observed;
            for (int i = n - 1; i > 0; --i)
                std::swap(labels[i], labels[rng.below(i + 1)]);
            return;
        }
        do {
            if (exhausted) Rcpp::stop("Assignment enumeration failed.");
            labels = enumeration;

            exhausted = !std::prev_permutation(enumeration.begin(), enumeration.end());
        } while (labels == observed);
    }
};

// Pooled standardized differences in proportions. Endpoint totals remain
// fixed under permutation, so all attainable statistics can be computed first.
class Scores {
    int m;
    std::vector<int> lower;
    std::vector<std::vector<int>> events;
    std::vector<std::vector<Key>> lookup;
public:
    Scores(const Rcpp::NumericMatrix& x, int n1, int alternative) :
        m(x.ncol()), lower(m), events(m), lookup(m) {
        const int n = x.nrow(), n0 = n - n1;
        for (int j = 0; j < m; ++j) {
            if (j % 128 == 0) Rcpp::checkUserInterrupt();
            for (int i = 0; i < n; ++i) {
                const double value = x(i, j);
                if (value != 0.0 && value != 1.0)
                    Rcpp::stop("'x' must contain only 0 and 1, without missing values.");
                if (value == 1.0) events[j].push_back(i);
            }
            const int successes = events[j].size();
            lower[j] = std::max(0, successes - n0);
            const int upper = std::min(n1, successes);
            lookup[j].resize(upper - lower[j] + 1);
            const double pooled = static_cast<double>(successes) / n;
            const double se = std::sqrt(pooled * (1.0 - pooled) *
                                       (1.0 / n1 + 1.0 / n0));
            for (int k = lower[j]; k <= upper; ++k) {
                double z = se == 0.0 ? 0.0 :
                    (static_cast<double>(k) / n1 -
                     static_cast<double>(successes - k) / n0) / se;
                if (alternative == 1) z = -z;
                if (alternative == 2) z = std::fabs(z);
                lookup[j][k - lower[j]] = static_cast<Key>(std::llround(z * scale));
            }
        }
    }

    void evaluate(const std::vector<int>& labels, std::vector<Key>& values) const {
        for (int j = 0; j < m; ++j) {
            const int successes = std::accumulate(events[j].begin(), events[j].end(), 0,
                [&](int sum, int subject) { return sum + labels[subject]; });
            values[j] = lookup[j][successes - lower[j]];
        }
    }
};

// Classical tails include all ties; shared-U counts the same random fraction
// of ties in every local test.
struct Count {
    int greater = 0, equal = 0;
    void add(Key value, Key observed) {
        greater += value > observed;
        equal += value == observed;
    }
    double p(int total, double u = 1.0) const {
        return (greater + u * equal) / total;
    }
};

struct Capture {
    bool keep;
    Rcpp::IntegerMatrix assignments;
    Rcpp::NumericMatrix statistics;
    Capture(bool save, int rows, int n, int m) : keep(save),
        assignments(save ? rows : 0, save ? n : 0),
        statistics(save ? rows : 0, save ? m : 0) {}
    void add(int row, const std::vector<int>& labels, const std::vector<Key>& values) {
        if (!keep) return;
        for (std::size_t i = 0; i < labels.size(); ++i) assignments(row, i) = labels[i];
        for (std::size_t j = 0; j < values.size(); ++j) statistics(row, j) = values[j] / scale;
    }
    void append(Rcpp::List& result) const {
        if (!keep) return;
        result["permutations"] = assignments;
        result["permutation_statistics"] = statistics;
    }
};

// Fix the observed endpoint order, resolving ties by column order.
inline std::vector<int> ranking(const std::vector<Key>& values) {
    std::vector<int> order(values.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
        [&](int a, int b) { return values[a] > values[b]; });
    return order;
}

// Use the same permutation list for every intersection.
class PermutationProblem {
    Permutations permutations;
    Scores scores;
public:
    int m, total;
    std::vector<Key> observed, current;
    std::vector<int> order;
    Capture capture;

    PermutationProblem(const Rcpp::NumericMatrix& x, const Rcpp::IntegerVector& group,
                       int n_perm, bool exact, int perm_seed, int alternative, bool keep) :
        permutations(group, n_perm, exact, perm_seed),
        scores(x, Rcpp::sum(group), alternative), m(x.ncol()), total(permutations.total),
        observed(m), current(m), capture(keep, total, x.nrow(), m) {
        scores.evaluate(permutations.labels, observed);
        order = ranking(observed);
    }

    template<class Function>
    void for_each_assignment(Function function) {
        for (int row = 0; row < total; ++row) {
            if (row % 32 == 0) Rcpp::checkUserInterrupt();
            permutations.next(row);
            scores.evaluate(permutations.labels, current);
            capture.add(row, permutations.labels, current);
            function(row, current);
        }
    }

    Rcpp::List finish(Rcpp::List result) const {
        Rcpp::NumericVector statistic(m);
        std::transform(observed.begin(), observed.end(), statistic.begin(),
            [](Key value) { return value / scale; });
        result["statistic"] = statistic;
        result["n_assignments"] = total;
        capture.append(result);
        return result;
    }
};

// Evaluate suffix maxima, then take cumulative maxima of their local p-values.
Rcpp::List run_maxT(PermutationProblem& problem, double u, bool randomized, double alpha) {
    const int m = problem.m, total = problem.total;
    const auto& observed = problem.observed;
    const auto& order = problem.order;
    std::vector<Count> singles(m), maxima(m);

    problem.for_each_assignment([&](int, const std::vector<Key>& current) {

        Key maximum = missing;
        for (int a = m - 1; a >= 0; --a) {
            maximum = std::max(maximum, current[order[a]]);
            maxima[a].add(maximum, observed[order[a]]);
            singles[a].add(current[order[a]], observed[order[a]]);
        }
    });

    Rcpp::NumericVector marginal(m), adjusted(m), classic(m);

    double running = 0.0, running_classic = 0.0;
    for (int a = 0; a < m; ++a) {
        running = std::max(running, maxima[a].p(total, u));
        running_classic = std::max(running_classic, maxima[a].p(total));
        marginal[order[a]] = singles[a].p(total, u);
        adjusted[order[a]] = running;
        classic[order[a]] = running_classic;
    }
    Rcpp::List result = Rcpp::List::create(
        Rcpp::_["p_value"] = marginal,
        Rcpp::_["p_global"] = maxima[0].p(total, u),
        Rcpp::_["p_classic"] = classic);
    if (randomized) result["u"] = u;
    if (alpha < 0) {
        result["p_adjusted"] = adjusted;
    } else {

        Rcpp::LogicalVector reject(m);
        for (int j = 0; j < m; ++j) reject[j] = adjusted[j] <= alpha;
        result["reject"] = reject;
        result["p_adjusted_lower"] = adjusted;
        result["p_adjusted_upper"] = Rcpp::clone(adjusted);
    }
    return problem.finish(result);
}

// Use the second-largest statistic only to resolve ties in the maximum.
struct TopTwo {
    Key first = missing, second = missing;
    void insert(Key value) {
        if (value > first) { second = first; first = value; }
        else if (value > second) second = value;
    }
};

// Test singletons and completions formed by one endpoint and a later suffix.
Rcpp::List run_tea(PermutationProblem& problem) {
    const int m = problem.m, total = problem.total;
    const auto& observed = problem.observed;
    const auto& order = problem.order;
    std::vector<Count> singles(m), maxima(m);

    PairCounts pairs(m, m);
    pairs.clear();
    std::vector<TopTwo> suffix(m + 1);

    problem.for_each_assignment([&](int, const std::vector<Key>& current) {

        suffix[m] = TopTwo{};
        for (int a = m - 1; a >= 0; --a) {
            suffix[a] = suffix[a + 1];
            suffix[a].insert(current[order[a]]);
        }
        for (int a = 0; a < m; ++a) {
            if (a % 256 == 0) Rcpp::checkUserInterrupt();
            const Key threshold = observed[order[a]];
            singles[a].add(current[order[a]], threshold);
            maxima[a].add(suffix[a].first, threshold);
            for (int b = a + 1; b < m; ++b) {

                TopTwo value = suffix[b];
                value.insert(current[order[a]]);

                pairs(a, b) +=
                    value.first > threshold ||
                    (value.first == threshold && value.second >= observed[order[b]]);
            }
        }
    });

    Rcpp::NumericVector marginal(m), adjusted(m), classic(m);
    // Each completion contributes to its selected endpoint and its entire suffix.
    std::vector<double> anchor(m), starts(m, 0.0);
    for (int a = 0; a < m; ++a) {
        marginal[order[a]] = singles[a].p(total);
        anchor[a] = marginal[order[a]];
        for (int b = a + 1; b < m; ++b) {
            const double p = static_cast<double>(pairs(a, b)) / total;
            anchor[a] = std::max(anchor[a], p);
            starts[b] = std::max(starts[b], p);
        }
    }
    double running = 0.0, running_classic = 0.0;
    for (int a = 0; a < m; ++a) {
        running = std::max(running, starts[a]);
        running_classic = std::max(running_classic, maxima[a].p(total));
        adjusted[order[a]] = std::max(anchor[a], running);
        classic[order[a]] = running_classic;
    }
    const double global = m == 1 ? marginal[0] : static_cast<double>(pairs(0, 1)) / total;
    return problem.finish(Rcpp::List::create(
        Rcpp::_["p_value"] = marginal,
        Rcpp::_["p_adjusted"] = adjusted,
        Rcpp::_["p_global"] = global,
        Rcpp::_["p_classic"] = classic));
}

// Evaluate top-k completion tests, using bounds to limit the closure search.
class KeaClosure {
    int endpoint_count, assignment_count, order_k, cutoff;
    bool decision_mode;
    const std::vector<Key>& observed;
    const ScoreMatrix& assignment_scores;
    const ScoreCube& suffix_top_scores;
    std::vector<int> lower_counts, upper_counts, prefix_ranks;
    double nodes_evaluated = 0, branches_pruned = 0;
    int global_count = 0;

    struct NodeCounts {
        int completion = 0;
        int bound = 0;
        int prefix_only = 0;
    };

    bool needs_improvement(int rank) const {
        return decision_mode
            ? lower_counts[rank] <= cutoff && upper_counts[rank] > cutoff
            : lower_counts[rank] < upper_counts[rank];
    }

    int pruning_target(int last_rank) const {
        int target = std::numeric_limits<int>::max();
        for (int rank : prefix_ranks)
            if (needs_improvement(rank))
                target = std::min(target, decision_mode ? cutoff : lower_counts[rank]);
        for (int rank = last_rank + 1; rank < endpoint_count; ++rank)
            if (needs_improvement(rank))
                target = std::min(target, decision_mode ? cutoff : lower_counts[rank]);
        return target;
    }

    void update_lower_counts(int last_rank, int count, bool include_tail) {
        for (int rank : prefix_ranks)
            lower_counts[rank] = std::max(lower_counts[rank], count);
        if (include_tail)
            for (int rank = last_rank + 1; rank < endpoint_count; ++rank)
                lower_counts[rank] = std::max(lower_counts[rank], count);
    }

    std::vector<Key> observed_threshold(int last_rank) const {
        const int prefix_size = prefix_ranks.size();
        const int length = std::min(order_k,
            prefix_size + endpoint_count - last_rank - 1);
        std::vector<Key> threshold(length);
        for (int position = 0; position < prefix_size; ++position)
            threshold[position] = observed[prefix_ranks[position]];
        for (int position = prefix_size; position < length; ++position)
            threshold[position] = observed[last_rank + 1 + position - prefix_size];
        return threshold;
    }

    // Compare the ordered statistics lexicographically, including remaining ties.
    NodeCounts count_node(int last_rank, const ScoreMatrix& prefix_scores,
                         const std::vector<Key>& threshold) const {
        const int prefix_size = prefix_ranks.size();
        const int tail_size = endpoint_count - last_rank - 1;
        NodeCounts counts;
        for (int row = 0; row < assignment_count; ++row) {
            const auto tail = suffix_top_scores[last_rank + 1][row];
            int prefix_position = 0, tail_position = 0;
            int comparison = 0, comparison_at_prefix = 0;
            for (int position = 0; position < static_cast<int>(threshold.size()); ++position) {
                Key value;
                if (prefix_position < prefix_size &&
                    (tail_position >= tail_size ||
                     prefix_scores[row][prefix_position] >= tail[tail_position]))
                    value = prefix_scores[row][prefix_position++];
                else
                    value = tail[tail_position++];
                if (comparison == 0)
                    comparison = (value > threshold[position]) - (value < threshold[position]);
                if (position + 1 == prefix_size) comparison_at_prefix = comparison;
                if (comparison != 0 && position + 1 >= prefix_size) break;
            }
            counts.completion += comparison >= 0;
            counts.bound += prefix_size == 0 || comparison_at_prefix >= 0;

            // Intersections smaller than k use all their endpoint statistics.
            if (prefix_size > 1 && prefix_size < order_k) {
                const auto first = prefix_scores[row].begin();
                counts.prefix_only += !std::lexicographical_compare(
                    first, first + prefix_size,
                    threshold.begin(), threshold.begin() + prefix_size);
            }
        }
        return counts;
    }

    void make_child_scores(int rank, const ScoreMatrix& prefix_scores,
                           int prefix_size, ScoreMatrix& child_scores) const {
        for (int row = 0; row < assignment_count; ++row) {
            const Key value = assignment_scores[rank][row];
            if (prefix_size == 0) {
                child_scores[row][0] = value;
                continue;
            }
            const auto source = prefix_scores[row].begin();
            const auto target = child_scores[row].begin();
            const auto insertion = std::upper_bound(
                source, source + prefix_size, value, std::greater<Key>());
            const auto position = insertion - source;
            std::copy(source, insertion, target);
            target[position] = value;
            std::copy(insertion, source + prefix_size, target + position + 1);
        }
    }

    void visit(int last_rank, const ScoreMatrix& prefix_scores, int inherited_bound) {
        const int prefix_size = prefix_ranks.size();
        if (prefix_size != 0 && inherited_bound <= pruning_target(last_rank)) {
            ++branches_pruned;
            return;
        }
        if (static_cast<std::uint64_t>(++nodes_evaluated) % 1024 == 0)
            Rcpp::checkUserInterrupt();

        const std::vector<Key> threshold = observed_threshold(last_rank);
        NodeCounts counts = count_node(last_rank, prefix_scores, threshold);
        counts.bound = std::min(counts.bound, inherited_bound);
        if (prefix_size == 0) global_count = counts.completion;
        update_lower_counts(last_rank, counts.completion, true);
        if (prefix_size > 1 && prefix_size < order_k)
            update_lower_counts(last_rank, counts.prefix_only, false);

        // Stop when further intersections cannot change the p-values or decisions.
        if (prefix_size == order_k || counts.bound <= pruning_target(last_rank)) {
            ++branches_pruned;
            return;
        }

        ScoreMatrix child_scores(boost::extents[assignment_count][prefix_size + 1]);
        for (int rank = last_rank + 1; rank < endpoint_count; ++rank) {
            prefix_ranks.push_back(rank);
            const int child_bound = std::min(counts.bound, upper_counts[rank]);
            if (child_bound > pruning_target(rank)) {
                make_child_scores(rank, prefix_scores, prefix_size, child_scores);
                visit(rank, child_scores, child_bound);
            } else {
                ++branches_pruned;
            }
            prefix_ranks.pop_back();
        }
    }

public:
    KeaClosure(int endpoints, int assignments, int k, double alpha,
               const std::vector<Key>& observed_scores,
               const ScoreMatrix& all_scores,
               const ScoreCube& suffix_scores,
               const std::vector<int>& singleton_counts,
               const std::vector<int>& certified_upper_counts)
        : endpoint_count(endpoints), assignment_count(assignments), order_k(k),
          cutoff(-1), decision_mode(alpha >= 0),
          observed(observed_scores), assignment_scores(all_scores),
          suffix_top_scores(suffix_scores), lower_counts(singleton_counts),
          upper_counts(certified_upper_counts) {
        if (decision_mode) {
            cutoff = static_cast<int>(std::floor(alpha * assignment_count));

            while (cutoff < assignment_count &&
                   static_cast<double>(cutoff + 1) / assignment_count <= alpha) ++cutoff;
            while (cutoff >= 0 &&
                   static_cast<double>(cutoff) / assignment_count > alpha) --cutoff;
        }
    }

    Rcpp::List run(const std::vector<int>& order) {
        ScoreMatrix empty_prefix(boost::extents[assignment_count][0]);
        visit(-1, empty_prefix, assignment_count);
        Rcpp::NumericVector lower(endpoint_count), upper(endpoint_count);
        Rcpp::LogicalVector reject(endpoint_count);
        for (int rank = 0; rank < endpoint_count; ++rank) {
            lower[order[rank]] = static_cast<double>(lower_counts[rank]) / assignment_count;
            reject[order[rank]] = lower_counts[rank] <= cutoff;
            const int upper_count = decision_mode && reject[order[rank]]
                ? std::min(upper_counts[rank], cutoff) : upper_counts[rank];
            upper[order[rank]] = static_cast<double>(upper_count) / assignment_count;
        }
        Rcpp::List result = Rcpp::List::create(
            Rcpp::Named("p_global") = static_cast<double>(global_count) / assignment_count,
            Rcpp::Named("nodes_evaluated") = nodes_evaluated,
            Rcpp::Named("branches_pruned") = branches_pruned);
        if (decision_mode) {
            result["reject"] = reject;
            result["p_adjusted_lower"] = lower;
            result["p_adjusted_upper"] = upper;
        } else {
            result["p_adjusted"] = lower;
        }
        return result;
    }
};

Rcpp::List run_kea(PermutationProblem& problem, int k, double alpha) {
    const int full_m = problem.m;
    const int total = problem.total;
    const auto& observed = problem.observed;
    const auto& endpoint_order = problem.order;

    ScoreMatrix all_T(boost::extents[full_m][total]);
    std::vector<int> singleton_counts(full_m, 0);
    problem.for_each_assignment([&](int row, const std::vector<Key>& current) {
        for (int rank = 0; rank < full_m; ++rank) {
            const Key value = current[endpoint_order[rank]];
            all_T[rank][row] = value;
            singleton_counts[rank] += value >= observed[endpoint_order[rank]];
        }
    });

    // Keep up to k copies of identical score vectors to preserve top-k multiplicity.
    std::map<std::vector<Key>, int> groups;
    std::vector<int> kept, representative(full_m), multiplicity;
    for (int rank = 0; rank < full_m; ++rank) {
        if (rank % 128 == 0) Rcpp::checkUserInterrupt();
        const auto begin = all_T[rank].begin();
        std::vector<Key> key(begin, begin + total);
        const auto inserted = groups.emplace(std::move(key), static_cast<int>(kept.size()));
        const int first = inserted.first->second;
        representative[rank] = first;
        if (inserted.second) {
            kept.push_back(rank);
            multiplicity.push_back(1);
        } else if (multiplicity[first] < k) {
            ++multiplicity[first];
            kept.push_back(rank);
            multiplicity.push_back(0);
        }
    }

    const int m = kept.size();
    std::vector<int> order(m), singles(m), classic(m, 0);
    std::iota(order.begin(), order.end(), 0);
    std::vector<Key> observed_ranked(m);
    ScoreMatrix T(boost::extents[m][total]);
    // Retain the k largest statistics of each suffix in every permutation.
    ScoreCube suffix_top(boost::extents[m + 1][total][k]);
    std::fill(suffix_top.data(), suffix_top.data() + suffix_top.num_elements(), missing);
    for (int rank = 0; rank < m; ++rank) {
        observed_ranked[rank] = observed[endpoint_order[kept[rank]]];
        singles[rank] = singleton_counts[kept[rank]];
        std::copy_n(all_T[kept[rank]].begin(), total, T[rank].begin());
    }
    for (int row = 0; row < total; ++row) {
        if (row % 32 == 0) Rcpp::checkUserInterrupt();
        for (int rank = m - 1; rank >= 0; --rank) {
            const Key value = T[rank][row];
            const auto here = suffix_top[rank][row].begin();
            const auto next = suffix_top[rank + 1][row].begin();
            const auto insertion = std::upper_bound(next, next + k, value, std::greater<Key>());
            const int position = insertion - next;
            std::copy(next, insertion, here);
            if (position < k) {
                here[position] = value;
                std::copy_n(insertion, k - position - 1, here + position + 1);
            }
            classic[rank] += here[0] >= observed_ranked[rank];
        }
    }

    std::partial_sum(classic.begin(), classic.end(), classic.begin(),
                     [](int left, int right) { return std::max(left, right); });

    // Lower augmentation orders provide upper bounds for higher orders.
    std::vector<int> bounds = classic;
    Rcpp::List result;
    Rcpp::NumericVector nodes_by_k(k);
    double nodes_evaluated = 0, branches_pruned = 0;
    const int first_order = k == 1 ? 1 : 2;
    for (int current_k = first_order; current_k <= k; ++current_k) {
        KeaClosure closure(m, total, current_k, alpha, observed_ranked, T, suffix_top, singles, bounds);
        result = closure.run(order);
        nodes_by_k[current_k - 1] = Rcpp::as<double>(result["nodes_evaluated"]);
        nodes_evaluated += nodes_by_k[current_k - 1];
        branches_pruned += Rcpp::as<double>(result["branches_pruned"]);
        Rcpp::NumericVector upper = alpha < 0
            ? Rcpp::as<Rcpp::NumericVector>(result["p_adjusted"])
            : Rcpp::as<Rcpp::NumericVector>(result["p_adjusted_upper"]);
        for (int rank = 0; rank < m; ++rank)
            bounds[rank] = static_cast<int>(std::llround(upper[order[rank]] * total));
    }
    result["nodes_evaluated"] = nodes_evaluated;
    result["branches_pruned"] = branches_pruned;
    result["nodes_by_k"] = nodes_by_k;

    for (const char* name : {"p_adjusted", "p_adjusted_lower", "p_adjusted_upper"}) {
        if (!result.containsElementNamed(name)) continue;
        Rcpp::NumericVector reduced = result[name], expanded(full_m);
        for (int rank = 0; rank < full_m; ++rank)
            expanded[endpoint_order[rank]] = reduced[representative[rank]];
        result[name] = expanded;
    }
    if (result.containsElementNamed("reject")) {
        Rcpp::LogicalVector reduced = result["reject"], expanded(full_m);
        for (int rank = 0; rank < full_m; ++rank)
            expanded[endpoint_order[rank]] = reduced[representative[rank]];
        result["reject"] = expanded;
    }
    Rcpp::NumericVector marginal(full_m), classic_expanded(full_m);
    for (int rank = 0; rank < full_m; ++rank) {
        marginal[endpoint_order[rank]] = static_cast<double>(singleton_counts[rank]) / total;
        classic_expanded[endpoint_order[rank]] = static_cast<double>(classic[representative[rank]]) / total;
    }
    result["p_value"] = marginal;
    result["p_classic"] = classic_expanded;
    result["n_score_groups"] = static_cast<int>(groups.size());
    result["n_endpoints_retained"] = m;
    return problem.finish(result);
}

}

// [[Rcpp::export(rng = false)]]
Rcpp::List tea_test_cpp(Rcpp::NumericMatrix x, Rcpp::IntegerVector group,
                       int n_perm, bool exact, int perm_seed, int u_seed,
                       int alternative, bool keep_permutations,
                       std::string method, int k, double alpha) {
    tea_detail::PermutationProblem problem(x, group, n_perm, exact, perm_seed,
                                          alternative, keep_permutations);
    if (method == "classic")
        return tea_detail::run_maxT(problem, 1.0, false, alpha);
    if (method == "shared_u") {
        tea_detail::Generator uniform_rng(u_seed, 0x554e4946U);
        return tea_detail::run_maxT(problem, uniform_rng.uniform(), true, alpha);
    }
    if (method == "tea") {
        if (alpha < 0) return tea_detail::run_tea(problem);
        return tea_detail::run_kea(problem, std::min(2, problem.m), alpha);
    }
    if (method == "kea") return tea_detail::run_kea(problem, k, alpha);
    Rcpp::stop("Unknown method.");
    return Rcpp::List();
}

