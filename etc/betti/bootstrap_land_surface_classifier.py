from __future__ import print_function

import numpy as np
from argparse import ArgumentParser
import os
import sys
import cPickle as pickle
import matplotlib.pyplot as plt

from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import classification_report

plt.style.use('seaborn-talk-ii')
DATA_DIR = '/home/licon/uni-koeln/tr32/stats/vertical_profiles'


def score_classifier(data, clf):

    X_train, X_test, y_train, y_test = data
    # clf = KNeighborsClassifier(k)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    rep = classification_report(y_test, y_pred, output_dict=True)

    return rep['weighted avg']['f1-score']


def boot(data, stat, clf, **kwargs):
    """
    stat : function that computes a statistic (or statistics) to compute the
    bootstrap for.
    """
    X_train, X_test, y_train, y_test = data
    n_samples = 1000
    sample_size = 1.0
    results = []

    for s in range(n_samples):
        # Sample from the data with replacement
        sample_idx = np.random.choice(
            X_train.shape[0],
            size=int(sample_size * X_train.shape[0]))
        X_sample = X_train[sample_idx]
        y_sample = y_train[sample_idx]

        statistic = stat(
            (X_sample, X_test, y_sample, y_test), clf, **kwargs)
        results.append(statistic)

    print('avg. f1 score: %0.3f +/- %0.2f' %
          (np.mean(results), np.std(results)))
    return (np.mean(results), np.std(results))


def pre_processor(dataset, varname, z_lims=None):
    """
    Take random subsets of both the timesteps (observations) and height
    levels (features), to circumvent the issue of autocorrelation
    idx_features = np.random.choice(range(1, 100), size=5, replace=False)
    Include only the timesteps between 10h and 17h, inclusive endpoints
    Also remove the first row, corresponding to the level next to the surface
    since this is all 0's

    Time range 1: timesteps 120 - 540 (10:00 - 17:00); take every 5th
    observation.

    z_lims: range of height levels to extract.

    """
    if z_lims is None:
        idx_features = range(1, 100) if varname != 'T' else range(0, 100)
    else:
        idx_features = range(z_lims[0], z_lims[1] + 1)
    idx_obs = np.array(range(120, 541, 5))
    sample_size = idx_obs.shape[0]

    if varname == 'b_0_1_pos':
        b0 = {k: dataset[k]['betti_0_pos'][
            idx_obs[:, None], idx_features] for k in dataset.keys()}
        b1 = {k: dataset[k]['betti_1_pos'][
            idx_obs[:, None], idx_features] for k in dataset.keys()}
        var_data = {k: np.concatenate([b0[k], b1[k]], axis=1)
                    for k in b0.keys()}
    elif varname == 'b_0_1_neg':
        b0 = {k: dataset[k]['betti_0_neg'][
            idx_obs[:, None], idx_features] for k in dataset.keys()}
        b1 = {k: dataset[k]['betti_1_neg'][
            idx_obs[:, None], idx_features] for k in dataset.keys()}
        var_data = {k: np.concatenate([b0[k], b1[k]], axis=1)
                    for k in b0.keys()}
    else:
        var_data = {k: dataset[k][varname][
            idx_obs[:, None], idx_features] for k in dataset.keys()}

    # Build data array for current variable
    X = np.concatenate([v for v in var_data.values()])
    y = np.array([0] * sample_size + [1] * sample_size +
                 [2] * sample_size + [3] * sample_size)
    # Standardize
    X_scaled = (X - np.mean(X, axis=0)[None, :]) / np.std(X, axis=0)[None, :]
    # Remove NaN columns
    X_scaled = X_scaled[:, ~np.isnan(X_scaled).any(axis=0)]

    # Train-test split
    train_size = int(0.7 * X_scaled.shape[0])
    train = np.random.choice(
        range(X_scaled.shape[0]), size=train_size, replace=False)
    test = [x for x in range(X_scaled.shape[0]) if x not in train]

    return (X_scaled[train], X_scaled[test], y[train], y[test])


def parse_args():

    ap = ArgumentParser(
        description='Train and evaluate land pattern classifiers.')
    ap.add_argument(
        '--clf', help='Classifier to use: kNN or logistic regression.',
        choices=['knn', 'logistic'], required=True)
    ap.add_argument('-C', help='Regularization value for logistic regression.',
                    default=1.0, type=float)
    ap.add_argument(
        '--all-c', help='Use all values of C in a pre-defined range',
        action='store_true')
    ap.add_argument('-k', help='Number of neighbors to use.', type=int,
                    default=5)
    ap.add_argument('--all-k', help='Use all values of k between 1 and 20.',
                    action='store_true')
    ap.add_argument(
        '--z-lims', help='Range of height levels to use as features',
        nargs=2, default=None, type=int)
    ap.add_argument('--out', help='Name of output file.', default=None)
    args = ap.parse_args()
    return args


if __name__ == '__main__':

    args = parse_args()
    # print(args)
    # Load data
    # Models for LES-ALM simulations
    simulations = ['SP1', 'SP2', 'SP3', 'SP4']
    names = ['betti_0_pos', 'betti_0_neg', 'betti_1_pos',
             'betti_1_neg', 'b_0_1_pos', 'b_0_1_neg', 'W', 'T']
    betti = {k: {} for k in simulations}
    # Load data
    for s in simulations:
        for n in names:
            if n != 'b_0_1_pos' and n != 'b_0_1_neg':
                data = np.load('%s/%s_%s.npy' % (DATA_DIR, s, n))
                betti[s][n] = data

    if args.clf == 'knn':

        if args.all_k:
            k_vals = range(1, 21)
        else:
            k_vals = [args.k]

        scores = {k: {} for k in k_vals}

        for k in k_vals:
            print('####    k = %d' % k)
            for n in names:
                clf = KNeighborsClassifier(k)
                data = pre_processor(betti, n, args.z_lims)
                print('bootstrapping for %-11s --- ' % n, end='')
                score = boot(data, score_classifier, clf)
                scores[k][n] = score
    elif args.clf == 'logistic':
        if args.all_c:
            c_vals = np.linspace(0.1, 2.0, 20)
        else:
            c_vals = [args.C]

        scores = {c: {} for c in c_vals}

        for c in c_vals:
            print('####    c = %0.2f' % c)
            for n in names:
                clf = LogisticRegression(
                    C=c, penalty='l2', tol=1e-6, max_iter=1000, solver='lbfgs',
                    multi_class='multinomial')
                data = pre_processor(betti, n, args.z_lims)
                print('bootstrapping for %-11s --- ' % n, end='')
                score = boot(data, score_classifier, clf)
                scores[c][n] = score

    if args.out is not None:
        with open(args.out, 'w') as f:
            pickle.dump(scores, f)
