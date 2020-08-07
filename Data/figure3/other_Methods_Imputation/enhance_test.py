# -*- coding:utf-8 _*-

import argparse
import hashlib
import os
import sys
import time
from math import ceil
from os import makedirs
from os.path import exists, basename, dirname, join
from time import time
from typing import Tuple, Iterable, Optional, Dict, List

import numpy as np
import scprep
from scipy.stats import poisson
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances


def enhance(
        X: np.ndarray,
        target_transcript_count: Optional[int] = 200000,
        max_neighbor_frac: Optional[float] = 0.02,
        pc_var_fold_thresh: Optional[float] = 2.0,
        max_components: Optional[int] = 50,
        k: Optional[int] = None,
        use_double_precision: Optional[bool] = False,
        seed: Optional[int] = 0) -> np.ndarray:
    """Remove technical noise from a scRNA-Seq expression matrix."""


    if use_double_precision:
        X = np.array(X, dtype=np.float64, order='C', copy=False)
    else:
        X = np.array(X, dtype=np.float32, order='C', copy=False)

    transcript_count = np.median(X.sum(axis=1))

    if k is None:
        k = int(ceil(target_transcript_count / transcript_count))

        k_max = int(max_neighbor_frac * X.shape[0])
        if k <= k_max:
            print(
                'Will perform denoising with k=%d '
                '(value was determined automatically '
                'based on a target transcript count of %d).',
                k, target_transcript_count)
        else:
            print(
                'Performing denoising with k=%d, to not exceed %.1f %% of '
                'the total number of cells. However, based on a target '
                'transcript count of %d, we should use k=%d. As a result, '
                'denoising results may be biased towards highly expressed '
                'genes.',
                k_max, 100 * max_neighbor_frac,
                target_transcript_count, k)
            k = k_max
    else:
        print('Will perform denoising with k=%d'
                     '(value was pre-specified).', k)


    num_components = determine_num_components(
        X, pc_var_fold_thresh, max_components,
        seed=seed)

    sys.stdout.flush()
    X_agg, cell_sizes = knn_aggregate(X, k, num_components, seed=seed)

    # denoise using PCA
    sys.stdout.flush()
    D, scores, components, mean = denoise_pca(
        X_agg, num_components, cell_sizes,
        seed=seed)

    return D, scores, components, mean, cell_sizes


def get_hash(X: np.ndarray) -> str:
    return str(hashlib.md5(X.data.tobytes()).hexdigest())


def determine_num_components(
        X: np.ndarray,
        var_fold_thresh: float = 2.0,
        max_components: int = 50,
        seed: int = 0) -> int:
    """Determine the number of significant principal components."""

    transcript_count = np.median(X.sum(axis=1))

    # apply PCA to real matrix
    _, real_pca_model = apply_pca(X, max_components, transcript_count, seed)

    # simulate pure noise matrix
    np.random.seed(seed)
    mean = normalize(X, transcript_count).mean(axis=0)
    X_noise = np.empty(X.shape, dtype=X.dtype)
    for i in range(X.shape[0]):
        X_noise[i, :] = poisson.rvs(mean)

    # apply PCA to pure noise matrix
    _, random_pca_model = apply_pca(X_noise, max_components, transcript_count, seed)
    var_thresh = var_fold_thresh * random_pca_model.explained_variance_[0]

    # determine number of components
    num_components = np.sum(real_pca_model.explained_variance_ >= var_thresh)

    return num_components


def normalize(
        X: np.ndarray,
        transcript_count: float = None) -> np.ndarray:
    """Perform median-normalization."""

    num_transcripts = X.sum(axis=1)

    if transcript_count is None:
        transcript_count = np.median(num_transcripts)

    N = ((transcript_count / num_transcripts) * X.T).T
    return N


def ft_transform(X: np.ndarray) -> np.ndarray:
    """Apply the Freeman-Tukey transformation."""

    # work around a bug where np.sqrt() says input is invalid for arrays
    # of type np.float32 that contain zeros
    invalid_errstate = 'warn'
    if np.issubdtype(X.dtype, np.float32):
        if np.amin(X) >= 0:
            invalid_errstate = 'ignore'
    with np.errstate(invalid=invalid_errstate):
        T = np.sqrt(X) + np.sqrt(X + 1)

    return T


def apply_pca(
        X, num_components: int = 50,
        transcript_count=None,
        seed: int = 0) -> Tuple[np.ndarray, PCA]:
    """Apply principal component analysis."""

    pca_model = PCA(
        n_components=num_components,
        svd_solver='randomized',
        random_state=seed)

    X_trans = ft_transform(normalize(X, transcript_count))
    scores = pca_model.fit_transform(X_trans)

    return scores, pca_model


def knn_aggregate(
        X: np.ndarray, k: int, num_components: int,
        seed: int = 0) -> np.ndarray:
    """Aggregate measurements from nearest neighbors."""

    transcript_count = np.median(X.sum(axis=1))

    scores, _ = apply_pca(X, num_components, transcript_count, seed=seed)
    X_agg, _ = aggregate_neighbors(X, scores, k)

    _, pca_model = apply_pca(X_agg, num_components, transcript_count, seed=seed)
    input_matrix = ft_transform(normalize(X, transcript_count))
    scores = pca_model.transform(input_matrix)
    X_agg, cell_sizes = aggregate_neighbors(X, scores, k)

    return X_agg, cell_sizes


def aggregate_neighbors(
        X: np.ndarray, scores: np.ndarray, k: int) \
        -> Tuple[np.ndarray, np.ndarray]:
    """Sub-routine for nearest neighbor aggregation."""

    num_transcripts = X.sum(axis=1)
    dtype = X.dtype

    # make sure score matrix is C-contiguous
    scores = np.array(scores, dtype=dtype, order='C', copy=False)

    # work around a bug where np.sqrt() says input is invalid for arrays
    # of type np.float32 that contain zeros
    invalid_errstate = 'warn'
    if np.issubdtype(scores.dtype, np.float32):
        invalid_errstate = 'ignore'
    with np.errstate(invalid=invalid_errstate):
        D = pairwise_distances(scores, n_jobs=1, metric='euclidean')

    S = np.argsort(D, axis=1, kind='mergesort')
    X_agg = np.empty(X.shape, dtype=dtype)
    cell_sizes = np.empty(X.shape[0], dtype=dtype)
    for i in range(X.shape[0]):
        ind = S[i, :k]
        X_agg[i, :] = np.sum(X[ind, :], axis=0, dtype=dtype)
        cell_sizes[i] = np.median(num_transcripts[ind])

    return X_agg, cell_sizes


def restore_matrix(
        scores: np.ndarray, components: np.ndarray,
        mean: np.ndarray, cell_sizes: np.ndarray) -> np.ndarray:
    """Restore the expression matrix from PCA results and cell sizes."""

    # transform from PC space to original space
    D = scores.dot(components)

    # add gene means
    D = D + mean

    # invert the Freeman-Tukey transform
    D[D < 1] = 1
    D = np.power(D, 2)
    D = np.power(D - 1, 2) / (4 * D)

    D = ((cell_sizes / D.sum(axis=1)) * D.T).T
    return D


def denoise_pca(
        X: np.ndarray, num_components: int,
        cell_sizes: np.ndarray,
        seed: int = 0) -> np.ndarray:
    """Denoise data using PCA."""

    scores, pca_model = apply_pca(X, num_components, seed=seed)
    components = pca_model.components_
    mean = pca_model.mean_

    D = restore_matrix(scores, components, mean, cell_sizes)

    return D, scores, components, mean


def write_factorized(
        file_path: str,
        scores: np.ndarray, components: np.ndarray,
        mean: np.ndarray, cell_sizes: np.ndarray,
        cells: Iterable[str], genes: Iterable[str],
        compressed: bool = True) -> None:
    """Write ENHANCE results in factorized form."""

    file_path = os.path.expanduser(file_path)

    data = {}
    data['scores'] = np.array(scores, copy=False)
    data['components'] = np.array(components, copy=False)
    data['mean'] = np.array(mean, copy=False)
    data['cell_sizes'] = np.array(cell_sizes, copy=False)
    data['cells'] = np.array(list(cells))
    data['genes'] = np.array(list(genes))

    if compressed:
        np.savez_compressed(file_path, **data)
    else:
        np.savez(file_path, **data)


def read_factorized(file_path: str) \
        -> Tuple[np.ndarray, List[str], List[str], Dict[str, np.ndarray]]:
    """Read ENHANCE output in factorized form."""
    file_path = os.path.expanduser(file_path)
    data = np.load(file_path)
    scores = data['scores']
    components = data['components']
    mean = data['mean']
    cell_sizes = data['cell_sizes']
    cells = data['cells'].tolist()
    genes = data['genes'].tolist()

    D = restore_matrix(scores, components, mean, cell_sizes)
    return D, cells, genes, data


def get_passed_time(func):
    def wrapper(*args, **kwargs):
        start = time()
        func(*args, **kwargs)
        end = time()
        return end - start

    return wrapper


def entry():
    parser = argparse.ArgumentParser()
    path2 = '~/10x/'
    parser.add_argument('--path', type=str, default=path2, help='data directory to run')
    parser.add_argument('--times', type=int, default=1, help='run how many times')
    args = parser.parse_args()
    log_dir = join(dirname(__file__), 'enhance_log')
    log_file = join(log_dir, 'log.txt')
    if not exists(log_dir):
        makedirs(log_dir)
    return log_file, args.path, args.times


@get_passed_time
def enhence_test(matrix, transcript_count=200000, max_neighbor_frac=0.02, pc_var_fold_thresh=2.0,
                 max_components=50, num_neighbors=50,
                 use_double_precision=True, seed=0):
    """
    matrix==> raw-> cells
    :return:
    """
    enhance(
        matrix,
        target_transcript_count=transcript_count,
        max_neighbor_frac=max_neighbor_frac,
        pc_var_fold_thresh=pc_var_fold_thresh,
        max_components=max_components,
        k=num_neighbors,
        use_double_precision=use_double_precision,
        seed=seed)


if __name__ == '__main__':
    log, directory, times = entry()
    # log, directory, times = './log', '../data/amount_100', 2
    print("reading data")
    print(directory)
    raw_matrix = scprep.io.load_10X(directory)
    cell_nums = [50000]
    print("runing enhance")
    #cell_nums = [1000, 5000, 10000, 50000, 100000, 500000, 1000000]
    with open(log, 'a') as f:
        for cell_num in cell_nums:
            f.write('amount_{0}:\n'.format(str(cell_num)))
            interval_list = []
            for i in range(times):
                index = np.random.randint(low=0, high=raw_matrix.shape[0], size=cell_num)
                emt_data_tmp = raw_matrix.values[index, :]
                intreval = enhence_test(emt_data_tmp)
                print(intreval)
                interval_list.append(intreval)
            interval_list_1 = [str(item) for item in interval_list]
            f.write(str(sum(interval_list) / times) + '\t' + '\t'.join(interval_list_1) + '\n')
            f.flush()
