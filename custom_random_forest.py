import numpy as np
from sklearn.tree import DecisionTreeClassifier
from sklearn.base import BaseEstimator
from concurrent.futures import ThreadPoolExecutor

class RandomForestClassifierCustom(BaseEstimator):
    def __init__(
        self, n_estimators = 10, max_depth = None, max_features = None, random_state = 42
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state
        self.trees = []
        self.feat_ids_by_tree = []

    def fit(self, X, y, n_jobs=1):
        self.classes_ = sorted(np.unique(y))
        if n_jobs == 1:
            for i in range(self.n_estimators):
                np.random.seed(self.random_state + i)
                if self.max_features is None:
                    feat_ids = np.arange(X.shape[1])
                else:
                    feat_ids = np.random.choice(range(X.shape[1]), size=self.max_features, replace=False)
                self.feat_ids_by_tree.append(feat_ids)
                indices = np.random.choice(range(X.shape[0]), X.shape[0], replace=True)
                X_bootstrap = X[indices]
                y_bootstrap = y[indices]
                tree = DecisionTreeClassifier(max_depth=self.max_depth, random_state=self.random_state)
                tree.fit(X_bootstrap[:, feat_ids], y_bootstrap)
                self.trees.append(tree)
        else:
            with ThreadPoolExecutor(max_workers=n_jobs) as executor:
                futures = []
                for i in range(self.n_estimators):
                    np.random.seed(self.random_state + i)
                    if self.max_features is None:
                        feat_ids = np.arange(X.shape[1])
                    else:
                        feat_ids = np.random.choice(range(X.shape[1]), size=self.max_features, replace=False)
                    self.feat_ids_by_tree.append(feat_ids)
                    indices = np.random.choice(range(X.shape[0]), X.shape[0], replace=True)
                    X_bootstrap = X[indices]
                    y_bootstrap = y[indices]
                    tree = DecisionTreeClassifier(max_depth=self.max_depth, random_state=self.random_state)
                    futures.append(executor.submit(tree.fit, X_bootstrap[:, feat_ids], y_bootstrap))
                for future in futures:
                    self.trees.append(future.result())
        return self

    def predict_proba(self, X, n_jobs=1):
        probas = np.zeros((X.shape[0], len(self.classes_)))
        if n_jobs == 1:
            for i, tree in enumerate(self.trees):
                feat_ids = self.feat_ids_by_tree[i]
                probas += tree.predict_proba(X[:, feat_ids])
        else:
            with ThreadPoolExecutor(max_workers=n_jobs) as executor:
                futures = []
                for i, tree in enumerate(self.trees):
                    feat_ids = self.feat_ids_by_tree[i]
                    futures.append(executor.submit(tree.predict_proba, X[:, feat_ids]))
                for future, tree in zip(futures, self.trees):
                    probas += future.result()
        probas /= len(self.trees)
        return probas
    
    def predict(self, X, n_jobs=1):
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)
        return predictions