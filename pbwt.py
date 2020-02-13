import os

import allel
import networkx as nx
import numpy as np
from scipy.spatial.distance import pdist
from intervaltree import Interval, IntervalTree

DIR = os.path.dirname(os.path.realpath(__file__))
DEFAULT_MACS_DIR = os.path.join(DIR, 'resources/panels/macs')
DEFAULT_SCRM_DIR = os.path.join(DIR, 'resources/panels/scrm')


def create_macs_panel_subset(in_file, m=None, n=None, out_file=None):
    with open(in_file, 'r') as f:
        lines = f.readlines()
    m0, n0 = int(lines[-5].split('\t')[-1]), int(lines[-4].split('\t')[-1])
    m, n = m if m else m0, n if n else n0
    assert 0 <= n < n0 and 0 < m <= m0
    out_file = out_file if out_file else os.path.join(DEFAULT_MACS_DIR, f'm{m}_n{n}.macs')
    lines = lines[:2] + [l[:m - m0 - 1] + '\n' for l in lines[2:2 + n]] + [f'TOTAL_SAMPLES:\t{m}\n',
                                                                           f'TOTAL_SITES:\t{n}\n',
                                                                           'BEGIN_SELECTED_SITES\n',
                                                                           lines[-2][:2 * n - 1] + '\n',
                                                                           'END_SELECTED_SITES\n']
    with open(out_file, 'w') as f:
        f.writelines(lines)


def read_vcf_panel(in_file):
    return np.transpose(allel.read_vcf(in_file)['calldata/GT'].reshape(110, 200))


def read_macs_panel(in_file):
    with open(in_file, 'r') as f:
        lines = f.readlines()
    m = int(lines[-5].split('\t')[-1])
    return np.array([list(l[-m - 1:-1])[::-1] for l in lines[2:-5]], dtype=int).transpose()


# X[i][j] => j'th site of i'th individual.
def build_prefix_array(X, k, a):
    M, N = X.shape
    u, v, a0, b = 0, 0, np.zeros((X[:, k] == 0).sum()), np.zeros((X[:, k] == 1).sum())
    for i in range(M):
        # y_i^k[k] => x_{a_k[i]}
        if X[a[i, k], k] == 0:
            a0[u], u = a[i, k], u + 1
        else:
            b[v], v = a[i, k], v + 1
    a[:, k + 1] = np.concatenate([a0, b])


def build_prefix_and_divergence_arrays(X, k, a, d):
    M, N = X.shape
    u, v, p, q = 0, 0, k + 1, k + 1
    num_zeros, num_ones = (X[:, k] == 0).sum(), (X[:, k] == 1).sum()
    a0, b, d0, e = np.zeros(num_zeros), np.zeros(num_ones), np.zeros(num_zeros), np.zeros(num_ones)
    for i in range(M):
        if d[i, k] > p:
            p = d[i, k]
        if d[i, k] > q:
            q = d[i, k]
        if X[a[i, k], k] == 0:
            a0[u], d0[u], u, p = a[i, k], p, u + 1, 0
        else:
            b[v], e[v], v, q = a[i, k], q, v + 1, 0
    a[:, k + 1] = np.concatenate([a0, b])
    d[:, k + 1] = np.concatenate([d0, e])


def build_pbwt(X):
    M, N = X.shape
    a, d = np.zeros((M, N + 1), dtype=int), np.zeros((M, N + 1), dtype=int)
    a[:, 0] = np.array(list(range(M)), dtype=int)
    for k in range(N):
        build_prefix_and_divergence_arrays(X, k, a, d)
    return a, d


# def report_matches(X, a, d, L=1):
#     M, N = X.shape
#     for k in range(N):
#         for i in range(1, M):
#             if k + 1 - d[i, k+1] >= L:
#                 if k == N-1 or X[a[i, k+1], k+1] != X[a[i-1, k+1], k+1]:
#                     print(f'Match from X[{a[i, k+1]}] to X[{a[i-1, k+1]}] from sites {d[i, k+1]} to {k}.')
#                     print(f'Length: {k + 1 - d[i, k+1]}')
#                     print(f'X[{a[i-1, k+1]}, {d[i, k+1]}:{k+1}]: {X[a[i-1, k+1], d[i, k+1]:k+1]}')
#                     print(f'X[{a[i, k+1]}, {d[i, k+1]}:{k+1}]: {X[a[i, k+1], d[i, k+1]:k+1]}')


def report_matches(X, a, d, L=1):
    M, N = X.shape
    total = 0
    for k in range(N):
        i = 1
        while i < M:
            if k + 1 - d[i, k+1] >= L and (k == N-1 or X[a[i, k+1], k+1] != X[a[i-1, k+1], k+1]):
                j = i + 1
                while j < M and d[j, k+1] == d[i, k+1] and (k == N-1 or X[a[j, k+1], k+1] != X[a[j-1, k+1], k+1]):
                    j += 1
                total += 1
                print(f'Block of matches between {a[i-1:j, k+1]} from sites {d[i, k+1]} to {k}')
                print(f'Number of individuals: {len(a[i-1:j, k+1])}')
                print(f'Length: {k + 1 - d[i, k+1]}')
                for x in range(i-1, j):
                    print(f'X[{a[x, k+1]}, {d[i, k+1]}:{k+1}]: {X[a[x, k+1], d[i, k+1]:k+1]}')
                print()
                i = j
            else:
                i += 1
    print(f'Total blocks: {total}')


def report_match_blocks(X, a, d, m=2, L=1):
    M, N = X.shape
    total, it = 0, IntervalTree()
    for k in range(N):
        i = 1
        while i < M:
            if k + 1 - d[i, k+1] >= L:
                j, max = i + 1, bool(k == N-1 or X[a[i, k+1], k+1] != X[a[i-1, k+1], k+1])
                while j < M and d[j, k+1] == d[i, k+1]:
                    if max is False and (k == N-1 or X[a[j, k+1], k+1] != X[a[j-1, k+1], k+1]):
                        max = True
                    j += 1
                if max and j - i + 1 >= m:
                    total += 1
                    it.add(Interval(
                        begin=d[i, k+1],
                        end=k+1,
                        data={
                            'a_i': i-1,
                            'a_j': j,
                            'len': k + 1 - d[i, k+1]
                        }
                    ))
                    print(f'Block of matches between {a[i-1:j, k+1]} from sites {d[i, k+1]} to {k}')
                    print(f'Number of individuals: {len(a[i-1:j, k+1])}')
                    print(f'Length: {k + 1 - d[i, k+1]}')
                    for x in range(i-1, j):
                        print(f'X[{a[x, k+1]}, {d[i, k+1]}:{k+1}]: {X[a[x, k+1], d[i, k+1]:k+1]}')
                    print()
                i = j
            else:
                i += 1
    print(f'Total blocks: {total}')

    it1 = it[1]
    for i in it[1]:
        print(f"{a[i.data['a_i']:i.data['a_j'], i.end]} -> {i.data['len']}")

    return it


# def report_long_matches(X, a, d, L=1):
#     M, N = X.shape
#     for k in range(1, N):
#         u, v, a0, b0 = 0, 0, [0]*M, [0]*M
#         for i in range(M):
#             # if d[i, k] > k - L:
#             if k - d[i, k] < L:
#             # if k - d[i, k] >= L:
#                 if u > 0 and v > 0:
#                     for i_u in range(u):
#                         for i_v in range(v):
#                             # report match from a[i_u] to b[i_v] ending at k
#                             print(f'match from a0[i_u = {i_u}] = {a0[i_u]} to b0[i_v = {i_v}] = {b0[i_v]} ending at k = {k}')
#                             print(f'u: {u}, v: {v}')
#                             print(f'd[i = {i}, k = {k}]: {d[i,k]}')
#                             # print(f'{X[a0[i_u], d[i, k]:k]}')
#                             # print(f'{X[b0[i_v], d[i, k]:k]}')
#                             # print('Match:')
#                             # print(f'k: {k}, d[i = {i}, k = {k}]: {d[i, k]}')
#                             # print(f'u: {u}, a0: {a0}, v: {v}, b0: {b0}')
#                             # print(f'Length: {k - d[i, k]}')
#                             # print(f'a[i_u = {i_u}]: {a0[i_u]}')
#                             # print(f'b[i_v = {i_v}]: {b0[i_v]}')
#                             # print(f':{X[a[a0[i_u], k], d[i, k]:k]}')
#                             # print(f'X[a[a0[i_u], k], d[i, k]:k]: {X[a[a0[i_u], k], d[i, k]:k]}')
#                             # print(f'X[a[b0[i_v], k], d[i, k]:k]: {X[a[b0[i_v], k], d[i, k]:k]}\n')
#                 u, v = 0, 0
#
#             if X[a[i, k], k] == 0:
#                 a0[u], u = a[i, k], u+1
#             else:
#                 b0[v], v = a[i, k], v+1
#             # print(f'u: {u}, v: {v} (after)')


def build_tcpbwt(X):
    M, N = X.shape
    a = np.zeros((M, N + 1), dtype=int)
    a[:, 0] = np.array(list(range(M)), dtype=int)
    for k in range(N):
        # Let p be a minimal index such that y[p] = 1.
        p = X[a[:, k], k].tolist().index(1)
        u, v, a0, b = 0, 0, np.zeros((X[:, k] == 0).sum()), np.zeros((X[:, k] == 1).sum())
        for i in range(M):
            if X[a[i, k], k] == 0:
                a0[u], u = a[i, k], u + 1
            else:
                b[v], v = a[i, k], v + 1
        # We put the block of ones starting from the position of the first appearance of 1 in the column.
        a[:, k + 1] = np.concatenate([a0[:p], b, a0[p:]])
    return a


def naive_distance(X):
    M, N = X.shape
    ret = np.array([np.sum(np.abs(X[i, :] - X[j, :])) for i in range(M) for j in range(i + 1, M) if i is not j])
    return ret


def naive_similarity(X):
    ret = naive_distance(X)
    return 1 - (ret / ret.max())


def compute_distances(X, metric='naive_distance'):
    assert metric in {'jaccard', 'naive_distance', 'naive_distance'}
    if metric == 'jaccard':
        return pdist(X, metric='jaccard')
    elif metric == 'naive_distance':
        return naive_distance(X)
    elif metric == 'naive_similarity':
        return naive_similarity(X)


def dist_index(i, j, n):
    assert i != j
    if i < j:
        i, j = j, i
    return int(n * j - j * (j + 1) / 2 + i - 1 - j)


def local_tree(X, a, k, metric='naive_distance'):
    M, N = X.shape
    assert 0 <= k < N
    sig = a[:, k + 1]
    distances = compute_distances(X, metric)
    d = np.array([distances[dist_index(sig[i], sig[i + 1], M)] for i in range(M - 1)], dtype=float)
    return sig, d


def local_trees(X, a, metric='naive_distance'):
    M, N = X.shape
    distances = compute_distances(X, metric)
    return [(a[:, k + 1], np.array([distances[dist_index(a[:, k + 1][i], a[:, k + 1][i + 1], M)] for i in range(M - 1)],
                                   dtype=float)) for k in range(N)]


def subtree_root(G, n):
    while G.in_edges(n):
        n = list(G.in_edges(n))[0][0]
    return n


def planar_encoding_to_networkx(sig, d):
    M = len(sig)
    G = nx.DiGraph()
    G.add_nodes_from(list(range(M)))
    d_argsort = np.argsort(d)
    for i in d_argsort:
        a, b, c = subtree_root(G, sig[i]), subtree_root(G, sig[i + 1]), list(G.nodes())[-1] + 1
        G.add_node(c)
        G.add_edges_from([(c, a), (c, b)])
    return G


def y(X, a, k=0):
    assert 0 <= k < len(X)
    return X[a[:, k], :]


def y_k(X, a, k=0):
    assert 0 <= k < len(X)
    return X[a[:, k], k]


def parse_interval(L, R, h, y, d):
    M, maximal_constant_subtrees = len(y), []
    hBr, HBr, allele = -1, -1, y[L]
    if L == R - 1:
        maximal_constant_subtrees.append({'L': L, 'R': R, 'allele': allele, 'h': -1, 'H': d[L]})
    else:
        i = L
        while i < h:
            j, D_L, hBr = i + 1, d[h], -1
            while j < M and d[j] > D_L:
                if hBr == -1 or HBr > d[j]:
                    hBr, HBr = j, d[j]
                j += 1
            maximal_constant_subtrees.append({'L': i, 'R': j, 'allele': allele, 'h': hBr, 'H': D_L})
            i = j
        i, rPack_tmp = R, []
        while i > h:
            j, D_R, hBr = i - 1, d[i], -1
            while d[j] < D_R:
                if hBr == -1 or HBr > d[j]:
                    hBr, HBr = j, d[j]
                j -= 1
            rPack_tmp.append((j, i, hBr))
            i = j
        for i in range(len(rPack_tmp) - 1, -1, -1):
            L, R, hBr = rPack_tmp[i]
            maximal_constant_subtrees.append({'L': i, 'R': j, 'allele': allele, 'h': hBr, 'H': d[L]})
    return maximal_constant_subtrees


def reduce_tree(y, d):
    d_list = d.tolist()
    d = np.array(d_list + [d_list[-1]] * 2)
    M, maximal_constant_subtrees = len(y), []
    L, H = 0, -1.0
    for i in range(M):
        if d[i] < H or H == -1:
            H, h = d[i], i
        if i == M - 1 or y[i] != y[i + 1]:
            if i == M - 1 or d[i + 1] < H:
                h += 1
            maximal_constant_subtrees += parse_interval(L, i + 1, h, y, d)
            H, L = -1, i + 1
    return maximal_constant_subtrees


def maximal_constant_subtrees_rec(L, R, y, sig, d):
    if R - L == 1:
        H = min(d[L - 1], d[L]) if L > 0 else d[L]
        return [{'L': L, 'R': R, 'allele': y[L], 'h': -1, 'H': H}]
    if np.all(y[L:R] == y[L]):
        h = np.argmax(d[L:R - 1]) + L
        return [{'L': L, 'R': R, 'allele': y[L], 'h': h, 'H': d[h]}]

    subtree_intervals = [L] + (np.argwhere(d[L:R - 1] == np.amax(d[L:R - 1])).flatten() + L + 1).tolist() + [R]
    subtree_intervals = [(subtree_intervals[i], subtree_intervals[i + 1]) for i in range(len(subtree_intervals) - 1)]

    return [a for b in [maximal_constant_subtrees_rec(L, R, y, sig, d) for L, R in subtree_intervals] for a in b]


def maximal_constant_subtrees(y, sig, d):
    return maximal_constant_subtrees_rec(0, len(y), y, sig, d)


# def tree_refinement_rec(L, R, y, sig, d_red, MCS, l):
#     print('\t' * l + f'L={L}, R={R}')
#     M = len(y)
#     if R - L == 1:
#         return [MCS[L]]
#     if np.all(d_red[L:R - 1] == d_red[L]):
#         num_ones = sum([subtree['allele'] for subtree in MCS[0:3]])
#         if num_ones > 1:
#
#             print('\t' * l + f'Subtree [{L}:{R}) has {R - L} total subtrees nodes & {num_ones} 1-subtrees.')
#
#             zero_subtree_indices = [i for i in range(L, R) if MCS[i]['allele'] == 0]
#             mcs_intervals_ref = [(L, zero_subtree_indices[0], 1)] if zero_subtree_indices[0] > L else []
#             for i in range(len(zero_subtree_indices) - 1):
#                 mcs_intervals_ref.append((zero_subtree_indices[i], zero_subtree_indices[i]+1, 0))
#                 if zero_subtree_indices[i+1] - zero_subtree_indices[i] > 1:
#                     mcs_intervals_ref.append((zero_subtree_indices[i] + 1, zero_subtree_indices[i + 1], 1))
#             mcs_intervals_ref.append((zero_subtree_indices[-1], zero_subtree_indices[-1]+1, 0))
#             if R - zero_subtree_indices[-1] > 1:
#                 mcs_intervals_ref.append((zero_subtree_indices[-1] + 1, R, 1))
#
#             mcs_ref = [MCS[L] if allele == 0 else {'L': L, 'R': R, 'allele': y[L], 'h': h, 'H': d[h]} for L, R, allele in mcs_intervals_ref]
#
#             print('\t' * l + f'mcs_intervals_ref={mcs_intervals_ref}')
#
#             # one_subtree_intervals = [L - 1] + [i for i in range(L, R) if MCS[i]['allele'] == 0] + [R]
#             # one_subtree_intervals = [(one_subtree_intervals[i] + 1, one_subtree_intervals[i + 1])
#             #                          for i in range(len(one_subtree_intervals) - 1)
#             #                          if one_subtree_intervals[i + 1] - one_subtree_intervals[i] - 1 > 0]
#             # print('\t' * l + f'one_subtree_intervals={one_subtree_intervals}')
#             return
#         else:
#             return [MCS[L:R]]
#
#     subtree_intervals = [L] + (np.argwhere(d_red[L:R - 1] == np.amax(d_red[L:R - 1])).flatten() + L + 1).tolist() + [R]
#     subtree_intervals = [(subtree_intervals[i], subtree_intervals[i + 1]) for i in range(len(subtree_intervals) - 1)]
#
#     for L, R in subtree_intervals:
#         tree_refinement_rec(L, R, y, sig, d_red, MCS, l + 1)


# def tree_refinement(y, sig, d, MCS=None):
#     if MCS is None:
#         MCS = maximal_constant_subtrees(y, sig, d)
#     M = len(y)
#     d_red = [d[subtree['R'] - 1] for subtree in MCS if subtree['R'] - 1 < M - 1]
#     return tree_refinement_rec(0, len(MCS), y, sig, d_red, MCS, 0)


def reduce_tree(MCS, sig, d):
    M, sig_red, d_red = len(sig), [], []
    for subtree in MCS:
        L, R, allele, h, H = subtree['L'], subtree['R'], subtree['allele'], subtree['h'], subtree['H']
        sig_red.append(sig[L:R])
        if R - 1 < M - 1:
            d_red.append(d[R - 1])
    return np.array(sig_red), np.array(d_red)


def sorted_sites_to_intervals(l):
    if l is None or len(l) == 0:
        return [{'L': -1, 'R': -1, 'len': 0}]
    if not isinstance(l, list):
        l = list(l)
    if len(l) == 1:
        return [{'L': l[0], 'R': l[0] + 1, 'len': 1}]
    L, R, ret = 0, 0, []
    for i in range(1, len(l)):
        if l[i] - l[i - 1] == 1:
            R = i
        else:
            ret.append({'L': l[L], 'R': l[R] + 1, 'len': l[R] + 1 - l[L]})
            L, R = i, i
    ret.append({'L': l[L], 'R': l[R] + 1, 'len': l[R] + 1 - l[L]})
    return ret


def extend_interval_compatible_sites(X, lt, interval):
    M, N = X.shape
    L, R, = interval['L'], interval['R']

    if L > 0:
        one_indices = np.argwhere(X[lt, L - 1] == 1).flatten()
        while len(one_indices) == one_indices[-1] - one_indices[0] + 1 and L > 0:
            L -= 1
            if L > 0:
                one_indices = np.argwhere(X[lt, L - 1] == 1).flatten()

    if R < N:
        one_indices = np.argwhere(X[lt, R] == 1).flatten()
        while len(one_indices) == one_indices[-1] - one_indices[0] + 1 and R < N:
            R += 1
            if R < N:
                one_indices = np.argwhere(X[lt, R] == 1).flatten()

    return {'L': L, 'R': R, 'len': R - L}


def extend_intervals_compatible_sites(X, lt, intervals):
    return [extend_interval_compatible_sites(X, lt, i) for i in intervals]


def _site_interval_to_planar_order_map(X, max_interval_at_site):
    M, N = X.shape
    L, R, ret = 0, 0, {}
    for i in range(1, N):
        if max_interval_at_site[i]['sig'] == max_interval_at_site[i - 1]['sig']:
            R = i
        else:
            ret[(L, R + 1)] = {'len': max_interval_at_site[L]['len'], 'sig': max_interval_at_site[L]['sig']}
            L, R = i, i
    ret[(L, R + 1)] = {'len': max_interval_at_site[L]['len'], 'sig': max_interval_at_site[L]['sig']}
    return ret


def site_interval_to_planar_order_map(X):
    M, N = X.shape

    # a[:,k] -- Planar ordering at site k-1.
    a, a_rev = build_tcpbwt(X), np.flip(build_tcpbwt(np.flip(X, axis=1)), axis=1)
    a_t, a_t_rev = np.transpose(a), np.transpose(a_rev)

    # lt_site_map[i] = j, where i is a tuple of a planar ordering and j is list of sites for which i is the planar
    # ordering from either the forward or backwards run of tcPBWT.
    lt_site_map = {}
    for i in range(1, N + 1):
        a_t_i, a_t_rev_i = tuple(a_t[i]), tuple(a_t_rev[i])
        if a_t_i not in lt_site_map:
            lt_site_map[a_t_i] = {i - 1}
        else:
            lt_site_map[a_t_i].add(i - 1)
        if a_t_rev_i not in lt_site_map:
            lt_site_map[a_t_rev_i] = {i - 1}
        else:
            lt_site_map[a_t_rev_i].add(i - 1)

    for i in lt_site_map.keys():
        lt_site_map[i] = list(lt_site_map[i])
        lt_site_map[i].sort()
        lt_site_map[i] = extend_intervals_compatible_sites(X=X, lt=i,
                                                           intervals=sorted_sites_to_intervals(lt_site_map[i]))

    # A local tree topology is uniquely defined by a planar ordering & distance vector.
    # Given an idempotent distance function d(i,j) == d(j,i), a planar ordering implies the distance vector.
    # lt_set -- Set of unique local trees (planar orderings).
    lt_set = set(lt_site_map.keys())

    max_interval_at_site = {}
    for sig in lt_site_map.keys():
        for interval in lt_site_map[sig]:
            for site in range(interval['L'], interval['R']):
                if site not in max_interval_at_site or interval['len'] > max_interval_at_site[site]['len']:
                    max_interval_at_site[site] = interval
                    max_interval_at_site[site]['sig'] = sig

    site_interval_to_planar_order = _site_interval_to_planar_order_map(X, max_interval_at_site)

    return lt_set, site_interval_to_planar_order


def is_comptible_ordering(a, b):
    assert len(a) == len(b)
    for site in range(len(a)):
        if not (a[site] == b[site] or a[site] + b[site] == 4 * (a[site] // 2) + 1):
            return False
    return True


def find_compatible_orderings_lt(lt, sites):
    for i in range(len(lt) - 1):
        if np.any(lt[i][0] != lt[i + 1][0]) and is_comptible_ordering(lt[i][0], lt[i + 1][0]):
            print(sites[i + 1])


def find_compatible_orderings_min_lt(site_interval_to_planar_order, sites):
    intervals = list(site_interval_to_planar_order.keys())
    intervals = sorted(intervals, key=lambda k: k[0])

    for i in range(len(intervals) - 1):
        if site_interval_to_planar_order[intervals[i]]['sig'] != site_interval_to_planar_order[intervals[i + 1]]['sig'] \
                and is_comptible_ordering(site_interval_to_planar_order[intervals[i]]['sig'],
                                          site_interval_to_planar_order[intervals[i + 1]]['sig']):
            print(sites[intervals[i + 1][0]])


if __name__ == '__main__':
    pass
