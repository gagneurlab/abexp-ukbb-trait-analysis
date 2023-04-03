# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python [conda env:anaconda-florian4]
#     language: python
#     name: conda-env-anaconda-florian4-py
# ---

# %%
from IPython.display import display

# %%
import os
import numpy as np
import pandas as pd

import json
import yaml

import pyspark
import pyspark.sql.types as t
import pyspark.sql.functions as f

# import glow

# %%
import polars as pl

# %%
import re
import patsy

# %%
import plotnine as pn
# import seaborn as sns

import matplotlib
import matplotlib.pyplot as plt

# %%
# %matplotlib inline
# %config InlineBackend.figure_format='retina'

# %%
# import os
# # os.environ["RAY_ADDRESS"] = os.environ.get("RAY_ADDRESS", 'ray://192.168.16.30:10001')
# os.environ["RAY_ADDRESS"] = 'ray://192.168.16.28:10001'
# os.environ["RAY_ADDRESS"]

# %% [raw]
#

# %%
# we need ray to efficiently share the covariates df in memory
# use_ray = True
use_ray = False

# %%
ray = None
ray_context = None
if use_ray:
    import ray
    from rep.notebook_init import init_ray
    ray_context = init_ray(
        plasma_store_memory_fraction=0.3
    )

    from rep.notebook_init import init_spark_on_ray
    spark = init_spark_on_ray(
        # executor_cores=128,
        executor_memory_overhead=0.9,
        # configs={
        #     "spark.default.parallelism": int(ray.cluster_resources()["CPU"]),
        #     "spark.sql.shuffle.partitions": 2048,
        # }
        enable_glow=False,
    )
else:
    from rep.notebook_init import init_spark
    spark = init_spark(enable_glow=False)

# %%
ray_context

# %%
spark

# %%
snakefile_path = os.getcwd() + "/../../Snakefile"
snakefile_path

# %%
# del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'associate__regression',
        default_wildcards={
            # "phenotype_col": "CAD_SOFT",
            # "phenotype_col": "c_reactive_protein",
            # "phenotype_col": "severe_LDL",
            # "phenotype_col": "Asthma",
            # "phenotype_col": "Triglycerides",
            # "phenotype_col": "IGF1",
            # "phenotype_col": "Direct_bilirubin",
            # "phenotype_col": "triglycerides_f30870_0_0",
            # "phenotype_col": "standing_height_f50_0_0",
            # "phenotype_col": "body_mass_index_bmi_f21001_0_0",
            # "phenotype_col": "systolic_blood_pressure_automated_reading_f4080_0_0",
            "phenotype_col": "HDL_cholesterol",
            # "feature_set": "LOFTEE_pLoF",
            # "feature_set": "max_AbExp",
            # "feature_set": "AbExp_all_tissues",
            "feature_set": "AbExp_best_tissue",
            # "feature_set": "LOFTEE_pLoF",
            # "covariates": "randomized_sex_age_genPC_CLMP_PRS",
            "covariates": "sex_age_genPC_BMI_smoking_CLMP_PRS",
            # "covariates": "sex_age_genPC_CLMP_PRS",
            # "covariates": "sex_age_genPC_CLMP",
            # "covariates": "sex_age_genPC",
        }
    )

# %%
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %% [markdown]
# # Load configuration

# %%
with open(snakemake.input["featureset_config"], "r") as fd:
    config = yaml.safe_load(fd)

# %%
print(json.dumps(config, indent=2, default=str))

# %%
phenotype_col = snakemake.wildcards["phenotype_col"]
phenotype_col

# %%
phenocode = config["covariates"]["phenocode"]
phenocode

# %%
restricted_formula = config["covariates"]["restricted_formula"]
print(restricted_formula)

# %%
# Should we only keep the most associating variable?
aggregate_variables = config.get("aggregate_variables", False)
aggregate_variables

# %% [markdown]
# # Broadcast/deref function definitions

# %%
import tempfile

# %%
tmpdir = tempfile.TemporaryDirectory()
tmpdir

# %%
from typing import Union
import cloudpickle

def spark_broadcast_feather(df: Union[pd.DataFrame, pl.DataFrame, pl.LazyFrame], tmpdir: tempfile.TemporaryDirectory=tmpdir) -> str:
    tmpfile = tempfile.mktemp(dir=tmpdir.name, suffix=".feather")
    if isinstance(df, pd.DataFrame):
        df.to_feather(tmpfile)
    elif isinstance(df, pl.LazyFrame):
        df.collect().write_ipc(tmpfile)
    elif isinstance(df, pl.DataFrame):
        df.write_ipc(tmpfile)
    else:
        raise ValueError(f"Unknown dataframe type: '{type(df)}'")
        
    
    # add file to all spark nodes
    spark.sparkContext.addFile(tmpfile)
    
    # file will be accessible on all hosts by its name
    sparkfilename = os.path.basename(tmpfile)
    
    return "spark://" + sparkfilename

def spark_broadcast_pickle(obj, tmpdir: tempfile.TemporaryDirectory=tmpdir) -> str:
    tmpfile = tempfile.mktemp(dir=tmpdir.name, suffix=".pkl")
    with open(tmpfile, "wb") as fd:
        cloudpickle.dump(obj, fd)
    
    # add file to all spark nodes
    spark.sparkContext.addFile(tmpfile)
    
    sparkfilename = os.path.basename(tmpfile)
    
    return "spark://" + sparkfilename

def get_sparkfile_path(obj: str):
    if obj.startswith("spark://"):
        # remove 'spark://' prefix
        sparkfilename = obj[8:]
    
    return pyspark.SparkFiles.get(sparkfilename)


# %%
def broadcast(obj, idempotent=True):
    if idempotent:
        # check if is already a reference and return
        if isinstance(obj, pyspark.Broadcast):
            return obj
        elif ray is not None and isinstance(obj, ray.ObjectRef):
            return obj
        elif isinstance(obj, str):
            return obj
    
    if use_ray:
        return ray.put(obj)
    else:
        return spark.sparkContext.broadcast(obj)


# %%
import pickle

def deref(obj):
    if isinstance(obj, pyspark.Broadcast):
        return obj.value
    elif ray is not None and isinstance(obj, ray.ObjectRef):
        with ray_context:
            return ray.get(obj)
    elif isinstance(obj, str) and obj.startswith("spark://"):
        # remove 'spark://' prefix
        sparkfilename = obj[8:]
        
        if sparkfilename.endswith(".feather"):
            return pd.read_feather(pyspark.SparkFiles.get(sparkfilename))
        elif sparkfilename.endswith(".pkl"):
            with open(pyspark.SparkFiles.get(sparkfilename), "rb") as fd:
                return pickle.load(fd)
    else:
        return obj


# %% [markdown]
# # Read covariates

# %%
snakemake.input["covariates_pq"]

# %%
snakemake.input["train_test_split_pq"]

# %%
snakemake.input["train_test_split_pq"]

# %%
train_test_split_df = pl.scan_parquet(snakemake.input["train_test_split_pq"]).select("individual", "fold")
train_test_split_df.schema

# %%
covariates_df = pl.scan_parquet(snakemake.input["covariates_pq"])
covariates_df = (
    covariates_df
    .join(train_test_split_df, on="individual")
    .filter(pl.col("fold") == pl.lit("association_testing"))
)
# covariates_df.schema

# %%
covariates_df_columns = covariates_df.columns

# %% [raw]
# display(
#     "Size of 'covariates_df': %.3fGb" % (covariates_df.to_arrow().nbytes / 1024**3)
# )

# %%
import sklearn

# perform data normalization
normalized_covariates_df = covariates_df.collect().to_pandas()
cols_to_normalize = list(set(covariates_df.columns).difference(phenotype_col))
# cols_to_normalize = ["age_when_attended_assessment_centre_f21003_0_0"]
cols_to_normalize = [col for col, dtype in normalized_covariates_df.dtypes[cols_to_normalize].to_dict().items() if (
    pd.api.types.is_float_dtype(dtype)
    or (pd.api.types.is_integer_dtype(dtype) and (np.abs(np.max(normalized_covariates_df[col])) > 2))
)]
display(cols_to_normalize)

normalized_covariates_df = normalized_covariates_df.assign(**{
    col: sklearn.preprocessing.StandardScaler().fit_transform(normalized_covariates_df[[col]])[:, 0] for col in cols_to_normalize
})

# %%
normalized_covariates_df

# %% [markdown]
# ## clumping

# %%
if config["covariates"]["add_clumping"]:
    clumping_variants_df = pd.read_parquet(snakemake.input["clumping_variants_pq"])
else:
    clumping_variants_df = None
clumping_variants_df


# %%
def get_variants_by_gene(clumping_variants_df, gene_id):
    included_vars = clumping_variants_df.query(f"gene == '{gene_id}'")["variant"].values
    return included_vars


# %%
# get_variants_by_gene(clumping_variants_df, gene_id="ENSG00000084674")

# %%
def format_formula(formula, keys, add_clumping=True, clumping_variants_df=clumping_variants_df) -> patsy.ModelDesc:
    if not isinstance(formula, patsy.ModelDesc):
        model_desc = patsy.ModelDesc.from_formula(formula)
    else:
        model_desc = formula
    
    if add_clumping:
        if not "gene" in keys:
            raise ValueError(f"missing gene in keys: '{keys}'!")

        gene_id = keys["gene"]

        variants = get_variants_by_gene(clumping_variants_df=clumping_variants_df, gene_id=gene_id)

        if len(variants) > 0:
            model_desc.rhs_termlist += [
                patsy.Term([patsy.EvalFactor(f"Q('{c}')")]) for c in variants
            ]
    
    return model_desc


# %% [raw]
# test_formula = format_formula(
#     formula=restricted_formula,
#     clumping_variants_df=clumping_variants_df,
#     add_clumping=True,
#     keys={
#         "gene": "ENSG00000160584",
#     }
# )

# %% [raw]
# print(test_formula)

# %% [raw]
# test_formula_2 = format_formula(
#     formula=restricted_formula,
#     clumping_variants_df=clumping_variants_df,
#     add_clumping=True,
#     keys={
#         "gene": "",
#     }
# )

# %%
# print(test_formula_2)

# %%
import patsy

def get_variables_from_formula(formula, lhs=True, rhs=True):
    if not isinstance(formula, patsy.ModelDesc):
        model_desc = patsy.ModelDesc.from_formula(formula)
    else:
        model_desc = formula

    covariate_cols = [
        *([factor.name() for term in model_desc.lhs_termlist for factor in term.factors] if lhs else []),
        *([factor.name() for term in model_desc.rhs_termlist for factor in term.factors] if rhs else []),
    ]
    
    # find Q($i) -> $i
    q_regex = re.compile(r"""Q\(['"](.*)['"]\)""")
    
    parsed_covariate_cols = []
    for cov in covariate_cols:
        q_matches = q_regex.findall(cov)
        
        if len(q_matches) > 0:
            parsed_covariate_cols += q_matches
        else:
            parsed_covariate_cols.append(cov)

    # deduplicate
    parsed_covariate_cols = list(dict.fromkeys(parsed_covariate_cols))
    
    return parsed_covariate_cols


# %% [raw]
# get_variables_from_formula(test_formula)

# %%
# broadcast_covariates_df = broadcast(covariates_df.to_arrow())
broadcast_clumping_variants_df = broadcast(clumping_variants_df)

# %%
# covariates_pq_path = snakemake.input["covariates_pq"]
# covariates_df_columns = covariates_df.columns

# %%
from pyarrow import feather
from pyarrow import parquet as pq
import pyarrow as pa

# %%
data_df_schema = pq.read_table(
    snakemake.input["covariates_pq"],
    use_threads=False,
    memory_map=True,
    columns=[phenotype_col],
).schema
data_df_schema

# %%
data_type = data_df_schema[data_df_schema.get_field_index(phenotype_col)].type
data_type

# %%
is_boolean_dtype = pa.types.is_boolean(data_type)
if is_boolean_dtype:
    print("Regression type is boolean!")
else:
    print("Regression type is continuous!")

# %%
# from sklearn.decomposition import PCA
# pca = PCA(n_components = X_train_std.shape[1])
# pca_data = pca.fit_transform(X_train_std)

# %%
# restricted_model = smf.logit(
#     restricted_formula,
#     data = data_df.to_pandas().astype({snakemake.wildcards["phenotype_col"]: "float32"})
# ).fit()
# # broadcast_restricted_model = spark.sparkContext.broadcast(restricted_model)

# %%
from scipy.stats.distributions import chi2
from statsmodels.discrete.discrete_model import BinaryResultsWrapper
import statsmodels.api as sm

def likelihood_ratio(ll0, ll1):
    return -2 * (ll0-ll1)

def lr_test(
    restricted_model: BinaryResultsWrapper, 
    full_model: BinaryResultsWrapper
) -> (float, float, float):
    """
    Adopted from:
    - https://stackoverflow.com/a/70488612/2219819
    - statsmodels.regression.linear_model.RegressionResults class
    """
    df_diff = full_model.df_model - restricted_model.df_model
    
    chi2_stat = likelihood_ratio(
        restricted_model.llf,
        full_model.llf,
    )
    p = chi2.sf(chi2_stat, df_diff)

    return (
        chi2_stat,
        p,
        df_diff,
    )

def prsquared_adj(binary_model: BinaryResultsWrapper):
    """
    Adopted from statsmodels.regression.linear_model.RegressionResults class:
    Adjusted pseudo R-squared for binary results.
    This is defined here as 1 - (`nobs`-1)/`df_resid` * (1-`prsquared`)
    if a constant is included and 1 - `nobs`/`df_resid` * (1-`prsquared`) if
    no constant is included.
    """
    return 1 - (np.divide(binary_model.nobs - binary_model.k_constant, binary_model.df_resid)
                * (1 - binary_model.prsquared))


# %%
import warnings
from copy import deepcopy
from importlib.resources import open_text
from math import sqrt

import numpy as np
import scipy.linalg
import scipy.special
from scipy.stats import chi2, norm
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.exceptions import ConvergenceWarning
from sklearn.preprocessing import LabelEncoder
from sklearn.utils.extmath import fast_logdet
from sklearn.utils.multiclass import check_classification_targets
from sklearn.utils.validation import check_is_fitted
from tabulate import tabulate


class FirthLogisticRegression(BaseEstimator, ClassifierMixin):
    """
    Logistic regression with Firth's bias reduction method.

    This is based on the implementation in the `logistf` R package. Please see the
    `logistf` reference and Heinze & Schemper (2002) for details about the procedure.

    Parameters
    ----------
    max_iter
        The maximum number of Newton-Raphson iterations.
    max_halfstep
        The maximum number of step-halvings in one Newton-Raphson iteration.
    max_stepsize
        The maximum step size - for each coefficient, the step size is forced to
        be less than max_stepsize.
    pl_max_iter
        The maximum number of Newton-Raphson iterations for finding profile likelihood
        confidence intervals.
    pl_max_halfstep
        The maximum number of step-halvings in one iteration for finding profile
        likelihood confidence intervals.
    pl_max_stepsize
        The maximum step size while finding PL confidence intervals.
    tol
        Convergence tolerance for stopping.
    fit_intercept
        Specifies if intercept should be added.
    skip_pvals
        If True, p-values will not be calculated. Calculating the p-values can be
        time-consuming if `wald=False` since the fitting procedure is repeated for each
        coefficient.
    skip_ci
        If True, confidence intervals will not be calculated. Calculating the confidence
        intervals via profile likelihoood is time-consuming.
    alpha
        Significance level (confidence interval = 1-alpha). 0.05 as default for 95% CI.
    wald
        If True, uses Wald method to calculate p-values and confidence intervals.
    test_vars
        Index or list of indices of the variables for which to calculate confidence
        intervals and p-values. If None, calculate for all variables. This option has
        no effect if wald=True.

    Attributes
    ----------
    bse_
        Standard errors of the coefficients.
    classes_
        A list of the class labels.
    ci_
        The fitted profile likelihood confidence intervals.
    coef_
        The coefficients of the features.
    intercept_
        Fitted intercept. If `fit_intercept = False`, the intercept is set to zero.
    loglik_
        Fitted penalized log-likelihood.
    n_iter_
        Number of Newton-Raphson iterations performed.
    pvals_
        p-values calculated by penalized likelihood ratio tests.
    has_converged_
        information whether the Newton-Raphson iterations converged

    References
    ----------
    Firth D (1993). Bias reduction of maximum likelihood estimates.
    Biometrika 80, 27–38.

    Heinze G, Schemper M (2002). A solution to the problem of separation in logistic
    regression. Statistics in Medicine 21: 2409-2419.
    """

    def __init__(
            self,
            max_iter=25,
            max_halfstep=0,
            max_stepsize=5,
            pl_max_iter=100,
            pl_max_halfstep=0,
            pl_max_stepsize=5,
            tol=0.0001,
            fit_intercept=True,
            skip_pvals=False,
            skip_ci=False,
            alpha=0.05,
            wald=False,
            test_vars=None,
    ):
        self.max_iter = max_iter
        self.max_stepsize = max_stepsize
        self.max_halfstep = max_halfstep
        self.pl_max_iter = pl_max_iter
        self.pl_max_halfstep = pl_max_halfstep
        self.pl_max_stepsize = pl_max_stepsize
        self.tol = tol
        self.fit_intercept = fit_intercept
        self.skip_pvals = skip_pvals
        self.skip_ci = skip_ci
        self.alpha = alpha
        self.wald = wald
        self.test_vars = test_vars

    def _more_tags(self):
        return {"binary_only": True}

    def _validate_input(self, X, y):
        if self.max_iter < 0:
            raise ValueError(
                f"Maximum number of iterations must be positive; "
                f"got max_iter={self.max_iter}"
            )
        if self.max_halfstep < 0:
            raise ValueError(
                f"Maximum number of step-halvings must >= 0; "
                f"got max_halfstep={self.max_iter}"
            )
        if self.tol < 0:
            raise ValueError(
                f"Tolerance for stopping criteria must be positive; got tol={self.tol}"
            )
        X, y = self._validate_data(X, y, dtype=np.float64, ensure_min_samples=2)
        check_classification_targets(y)

        self.classes_ = np.unique(y)
        if len(self.classes_) != 2:
            raise ValueError(f"Got {len(self.classes_)} - only 2 classes supported.")
        y = LabelEncoder().fit_transform(y).astype(X.dtype, copy=False)

        return X, y

    def fit(self, X, y):
        X, y = self._validate_input(X, y)
        if self.fit_intercept:
            X = np.hstack((X, np.ones((X.shape[0], 1))))

        self.coef_, self.loglik_, self.n_iter_, self.has_converged_ = _firth_newton_raphson(
            X, y, self.max_iter, self.max_stepsize, self.max_halfstep, self.tol
        )

        self.bse_ = _bse(X, self.coef_)

        if not self.skip_ci:
            if not self.wald:
                self.ci_ = _profile_likelihood_ci(
                    X=X,
                    y=y,
                    fitted_coef=self.coef_,
                    full_loglik=self.loglik_,
                    max_iter=self.pl_max_iter,
                    max_stepsize=self.pl_max_stepsize,
                    max_halfstep=self.pl_max_halfstep,
                    tol=self.tol,
                    alpha=self.alpha,
                    test_vars=self.test_vars,
                )
            else:
                self.ci_ = _wald_ci(self.coef_, self.bse_, self.alpha)
        else:
            self.ci_ = np.full((self.coef_.shape[0], 2), np.nan)

        # penalized likelihood ratio tests
        if not self.skip_pvals:
            if not self.wald:
                self.pvals_ = _penalized_lrt(
                    self.loglik_,
                    X,
                    y,
                    self.max_iter,
                    self.max_stepsize,
                    self.max_halfstep,
                    self.tol,
                    self.test_vars,
                )

            else:
                self.pvals_ = _wald_test(self.coef_, self.bse_)
        else:
            self.pvals_ = np.full(self.coef_.shape[0], np.nan)

        if self.fit_intercept:
            self.intercept_ = self.coef_[-1]
            self.coef_ = self.coef_[:-1]
        else:
            self.intercept_ = 0

        return self

    def summary(self, xname=None, tablefmt="simple"):
        """
        Prints a summary table.

        Parameters
        ----------
        xname
            Names for the X variables. Default is x1, x2, ... Must match the number of
            parameters in the model.
        tablefmt
            `tabulate` table format for output. Please see the documentation for
            `tabulate` for options.
        """
        check_is_fitted(self)
        if xname and len(xname) != len(self.coef_):
            raise ValueError(
                f"Length of xname ({len(xname)}) does not match the number of "
                f"parameters in the model ({len(self.coef_)})"
            )

        if not xname:
            xname = [f"x{i}" for i in range(1, len(self.coef_) + 1)]

        coef = list(self.coef_)
        if self.fit_intercept:
            xname.append("Intercept")
            coef.append(self.intercept_)

        headers = [
            "",
            "coef",
            "std err",
            f"[{self.alpha / 2}",
            f"{1 - self.alpha / 2}]",
            "p-value",
        ]
        table = zip(xname, coef, self.bse_, self.ci_[:, 0], self.ci_[:, 1], self.pvals_)
        table = tabulate(table, headers, tablefmt=tablefmt)
        table += "\n\n"
        table += f"Log-Likelihood: {round(self.loglik_, 4)}\n"
        table += f"Newton-Raphson iterations: {self.n_iter_}\n"
        print(table)
        if self.fit_intercept:
            xname.pop()
        return

    def decision_function(self, X):
        check_is_fitted(self)
        X = self._validate_data(X, reset=False)
        scores = X @ self.coef_ + self.intercept_
        return scores

    def predict(self, X):
        decision = self.decision_function(X)
        if len(decision.shape) == 1:
            indices = (decision > 0).astype(int)
        else:
            indices = decision.argmax(axis=1)
        return self.classes_[indices]

    def predict_proba(self, X):
        decision = self.decision_function(X)
        if decision.ndim == 1:
            decision = np.c_[-decision, decision]
        proba = scipy.special.expit(decision)
        return proba


def _firth_newton_raphson(X, y, max_iter, max_stepsize, max_halfstep, tol, mask=None):
    # see logistf reference manual for explanation of procedure
    coef = np.zeros(X.shape[1])
    for iter in range(1, max_iter + 1):
        preds = scipy.special.expit(X @ coef)
        XW = _get_XW(X, preds, mask)

        fisher_info_mtx = XW.T @ XW
        hat = _hat_diag(XW)
        U_star = np.matmul(X.T, y - preds + np.multiply(hat, 0.5 - preds))
        step_size = np.linalg.lstsq(fisher_info_mtx, U_star, rcond=None)[0]
        # if mask:
        #     step_size[mask] = 0

        # step-halving
        mx = np.max(np.abs(step_size)) / max_stepsize
        if mx > 1:
            step_size = step_size / mx  # restrict to max_stepsize
        coef_new = coef + step_size
        preds_new = scipy.special.expit(X @ coef_new)
        loglike = _loglikelihood(X, y, preds)
        loglike_new = _loglikelihood(X, y, preds_new)
        steps = 0
        while loglike < loglike_new:
            step_size *= 0.5
            coef_new = coef + step_size
            preds_new = scipy.special.expit(X @ coef_new)
            loglike_new = _loglikelihood(X, y, preds_new)
            steps += 1
            if steps == max_halfstep:
                warning_msg = "Step-halving failed to converge."
                warnings.warn(warning_msg, ConvergenceWarning, stacklevel=2)
                return coef_new, -loglike_new, iter, False

        if iter > 1 and np.linalg.norm(coef_new - coef) < tol:
            return coef_new, -loglike_new, iter, True

        coef += step_size
    warning_msg = (
        "Firth logistic regression failed to converge. Try increasing max_iter."
    )
    warnings.warn(warning_msg, ConvergenceWarning, stacklevel=2)
    return coef, -loglike_new, max_iter, False


def _loglikelihood(X, y, preds):
    # penalized log-likelihood
    XW = _get_XW(X, preds)
    fisher_info_mtx = XW.T @ XW
    penalty = 0.5 * fast_logdet(fisher_info_mtx)
    if not np.isfinite(penalty):
        penalty = 0
    return -1 * (np.sum(y * np.log(preds) + (1 - y) * np.log(1 - preds)) + penalty)


def _get_XW(X, preds, mask=None):
    # mask is 1-indexed because 0 == None
    rootW = np.sqrt(preds * (1 - preds))
    XW = rootW[:, np.newaxis] * X

    # is this equivalent??
    # https://github.com/georgheinze/logistf/blob/master/src/logistf.c#L150-L159
    if mask is not None:
        XW[:, mask] = 0
    return XW


def _get_aug_XW(X, preds, hats):
    rootW = np.sqrt(preds * (1 - preds) * (1 + hats))
    XW = rootW[:, np.newaxis] * X
    return XW


def _hat_diag(XW):
    # Get diagonal elements of the hat matrix
    qr, tau, _, _ = scipy.linalg.lapack.dgeqrf(XW, overwrite_a=True)
    Q, _, _ = scipy.linalg.lapack.dorgqr(qr, tau, overwrite_a=True)
    hat = np.einsum("ij,ij->i", Q, Q)
    return hat


def _bse(X, coefs):
    # se in logistf is diag(object$var) ^ 0.5, where var is the covariance matrix,
    # which is the inverse of the observed fisher information matrix
    # https://stats.stackexchange.com/q/68080/343314
    preds = scipy.special.expit(X @ coefs)
    XW = _get_XW(X, preds)
    fisher_info_mtx = XW.T @ XW
    return np.sqrt(np.diag(np.linalg.pinv(fisher_info_mtx)))


def _penalized_lrt(
        full_loglik, X, y, max_iter, max_stepsize, max_halfstep, tol, test_vars
):
    if test_vars is None:
        test_var_indices = range(X.shape[1])
    elif isinstance(test_vars, int):  # single index
        test_var_indices = [test_vars]
    else:  # list, tuple, or set of indices
        test_var_indices = sorted(test_vars)

    pvals = []
    for mask in test_var_indices:
        _, null_loglik, _, _ = _firth_newton_raphson(
            X,
            y,
            max_iter,
            max_stepsize,
            max_halfstep,
            tol,
            mask,
        )
        pvals.append(_lrt(full_loglik, null_loglik))
    if len(pvals) < X.shape[1]:
        pval_array = np.full(X.shape[1], np.nan)
        for idx, test_var_idx in enumerate(test_var_indices):
            pval_array[test_var_idx] = pvals[idx]
        return pval_array
    return np.array(pvals)


def _lrt(full_loglik, null_loglik):
    # in logistf: 1-pchisq(2*(fit.full$loglik-fit.i$loglik),1)
    lr_stat = 2 * (full_loglik - null_loglik)
    p_value = chi2.sf(lr_stat, df=1)
    return p_value


def _predict(X, coef):
    preds = scipy.special.expit(X @ coef)
    np.clip(preds, a_min=1e-15, a_max=1 - 1e-15, out=preds)
    return preds


def _profile_likelihood_ci(
        X,
        y,
        fitted_coef,
        full_loglik,
        max_iter,
        max_stepsize,
        max_halfstep,
        tol,
        alpha,
        test_vars,
):
    LL0 = full_loglik - chi2.ppf(1 - alpha, 1) / 2
    lower_bound = []
    upper_bound = []
    if test_vars is None:
        test_var_indices = range(fitted_coef.shape[0])
    elif isinstance(test_vars, int):  # single index
        test_var_indices = [test_vars]
    else:  # list, tuple, or set of indices
        test_var_indices = sorted(test_vars)
    for side in [-1, 1]:
        # for coef_idx in range(fitted_coef.shape[0]):
        for coef_idx in test_var_indices:
            coef = deepcopy(fitted_coef)
            for iter in range(1, max_iter + 1):
                preds = _predict(X, coef)
                loglike = -_loglikelihood(X, y, preds)
                XW = _get_XW(X, preds)
                hat = _hat_diag(XW)
                XW = _get_aug_XW(X, preds, hat)  # augmented data using hat diag
                fisher_info_mtx = XW.T @ XW
                U_star = np.matmul(X.T, y - preds + np.multiply(hat, 0.5 - preds))
                # https://github.com/georgheinze/logistf/blob/master/src/logistf.c#L780-L781
                inv_fisher = np.linalg.pinv(fisher_info_mtx)
                tmp1x1 = U_star @ np.negative(inv_fisher) @ U_star
                underRoot = (
                        -2
                        * ((LL0 - loglike) + 0.5 * tmp1x1)
                        / (inv_fisher[coef_idx, coef_idx])
                )
                lambda_ = 0 if underRoot < 0 else side * sqrt(underRoot)
                U_star[coef_idx] += lambda_

                step_size = np.linalg.lstsq(fisher_info_mtx, U_star, rcond=None)[0]
                mx = np.max(np.abs(step_size)) / max_stepsize
                if mx > 1:
                    step_size = step_size / mx  # restrict to max_stepsize
                coef += step_size
                loglike_old = deepcopy(loglike)

                for halfs in range(1, max_halfstep + 1):
                    preds = _predict(X, coef)
                    loglike = -_loglikelihood(X, y, preds)
                    if (abs(loglike - LL0) < abs(loglike_old - LL0)) and loglike > LL0:
                        break
                    step_size *= 0.5
                    coef -= step_size
                if abs(loglike - LL0) <= tol:
                    if side == -1:
                        lower_bound.append(coef[coef_idx])
                    else:
                        upper_bound.append(coef[coef_idx])
                    break
            if abs(loglike - LL0) > tol:
                if side == -1:
                    lower_bound.append(np.nan)
                else:
                    upper_bound.append(np.nan)
                warning_msg = (
                    f"Non-converged PL confidence limits - max number of "
                    f"iterations exceeded for variable x{coef_idx}. Try "
                    f"increasing pl_max_iter."
                )
                warnings.warn(warning_msg, ConvergenceWarning, stacklevel=2)
    bounds = np.column_stack([lower_bound, upper_bound])
    if len(lower_bound) < fitted_coef.shape[0]:
        ci = np.full([fitted_coef.shape[0], 2], np.nan)
        for idx, test_var_idx in enumerate(test_var_indices):
            ci[test_var_idx] = bounds[idx]
        return ci

    return bounds


def _wald_ci(coef, bse, alpha):
    lower_ci = coef + norm.ppf(alpha / 2) * bse
    upper_ci = coef + norm.ppf(1 - alpha / 2) * bse
    return np.column_stack([lower_ci, upper_ci])


def _wald_test(coef, bse):
    # 1 - pchisq((beta^2/vars), 1), in our case bse = vars^0.5
    return chi2.sf(np.square(coef) / np.square(bse), 1)


def load_sex2():
    """
    Load the sex2 dataset from `logistf`.

    Returns
    -------
    X
        sex2 data as numpy array
    y
        sex2 `case` target column
    feature_names
        List of feature names

    References
    ----------
    Cytel Inc., (2010) LogXact 9 user manual, Cambridge, MA:Cytel Inc
    """
    with open_text("firthlogist.datasets", "sex2.csv") as sex2:
        X = np.loadtxt(sex2, skiprows=1, delimiter=",")
    y = X[:, 0]
    X = X[:, 1:]
    feature_names = ["age", "oc", "vic", "vicl", "vis", "dia"]
    return X, y, feature_names


def load_endometrial():
    """
    Load the endometrial cancer dataset analyzed in Heinze and Schemper (2002). The data
    was originally provided by Dr E. Asseryanis from the Vienna University Medical
    School

    Returns
    -------
    X
        endometrial data as numpy array
    y
        endometrial `HG` target column
    feature_names
        List of feature names

    References
    ----------
    Agresti, A (2015). Foundations of Linear and Generalized Linear Models.
    Wiley Series in Probability and Statistics.

    Heinze G, Schemper M (2002). A solution to the problem of separation in logistic
    regression. Statistics in Medicine 21: 2409-2419.
    """
    with open_text("firthlogist.datasets", "endometrial.csv") as sex2:
        X = np.loadtxt(sex2, skiprows=1, delimiter=",")
    y = X[:, -1]
    X = X[:, :-1]
    feature_names = ["NV", "PI", "EH"]
    return X, y, feature_names



# %%
from statsmodels.regression.linear_model import (
    RegressionModel,
    RegressionResults,
)

class Statsmodels_FirthLogisticRegression(RegressionModel):
    """
    Logistic regression with Firth's bias reduction method.

    """

    def __init__(self, endog, exog, column_names, **kwargs):
        super().__init__(endog, exog, **kwargs)

        self.model = FirthLogisticRegression(**kwargs, skip_pvals=False, skip_ci=True, wald=True, fit_intercept=False)
        self.column_names = column_names

    def fit(self, *args, **kwargs):
        self.model.fit(self.exog, self.endog)
        return self

    def whiten(self, data):
        return data

    def loglike(self, params):
        XW = _get_XW(self.exog, params)
        return _loglikelihood(self.exog, self.endog, preds=XW)

    # def score(self, params):
    #     pass
    #
    # def information(self, params):
    #     pass
    #
    # def hessian(self, params):
    #     pass

    @property
    def llf(self):
        return self.model.loglik_

    @property
    def mle_retvals(self):
        return {
            "converged": self.model.has_converged_
        }

    @property
    def null(self):
        """
        Fitted values of the null model
        """
        endog = self.endog
        model = self.model
        exog = np.ones((len(endog), 1))

        fitted = FirthLogisticRegression(skip_pvals=True, skip_ci=True, fit_intercept=False).fit(exog, endog)

        return fitted

    @property
    def llnull(self):
        """
        Log-likelihood of the model fit with a constant as the only regressor
        """
        return self.null.loglik_

    @property
    def prsquared(self):
        prsq = 1 - np.exp((self.llnull - self.llf) * (2 / self.nobs))
        return prsq
    
    @property
    def pvalues(self):
        retval = pd.Series(self.model.pvals_, index=self.column_names)
        return retval

    @property
    def params(self):
        retval = pd.Series(self.model.coef_, index=self.column_names)
        return retval



# %%
def firth_regression(formula, data):
    from patsy import dmatrices
    y, X = dmatrices(formula, data)
    
    fl = Statsmodels_FirthLogisticRegression(y, X, column_names=X.design_info.column_names, max_iter=200)
    
    return fl


# %% [raw]
# def _fit(
#     formatted_restricted_formula,
#     formatted_full_formula,
#     data_df,
#     boolean_regression: bool=False,
#     enable_firth_regression=True,
#     fitted_restricted_model=None,
#     fitted_full_model=None,
# ):
#     if boolean_regression:
#         if fitted_restricted_model is None:
#             # fit restricted model
#             restricted_model_converged = False
#             restricted_model = firth_regression(
#                 formatted_restricted_formula,
#                 data = data_df
#             )
#             try:
#                 restricted_model = restricted_model.fit()
#                 restricted_model_converged = restricted_model.mle_retvals["converged"]
#             except np.linalg.LinAlgError as e:
#                 pass
#         else:
#             restricted_model = fitted_restricted_model
#             restricted_model_converged = True
#
#         # fit full model
#         full_model_converged = False
#         full_model = firth_regression(
#             formatted_full_formula,
#             data = data_df
#         )
#         try:
#             if fitted_full_model is None:
#                 full_model = full_model.fit()
#             else:
#                 full_model = fitted_full_model
#             full_model = full_model.fit()
#             full_model_converged = full_model.mle_retvals["converged"]
#         except np.linalg.LinAlgError as e:
#             pass
#
# #             converged = False
# #             try:
# #                 restricted_model = smf.logit(
# #                     formatted_restricted_formula,
# #                     data = data_df
# #                 ).fit(maxiter=200, disp=0, warn_convergence=False)
#
# #                 full_model = smf.logit(
# #                     formatted_full_formula,
# #                     data = data_df
# #                 ).fit(maxiter=200, disp=0, warn_convergence=False)
#
# #                 converged = restricted_model.mle_retvals["converged"] & full_model.mle_retvals["converged"]
# #             except np.linalg.LinAlgError as linalg_error:
# #                 # newton failed
# #                 pass
#
# #             if not converged:
# #                 # retry with l-bfgs-b
# #                 print("retry with l-bfgs-b")
# #                 restricted_model = smf.logit(
# #                     formatted_restricted_formula,
# #                     data = data_df
# #                 ).fit(method="lbfgs", maxiter=1000, disp=0, warn_convergence=False)
#
# #                 full_model = smf.logit(
# #                     formatted_full_formula,
# #                     data = data_df
# #                 ).fit(method="lbfgs", maxiter=1000, disp=0, warn_convergence=False)
#
#
#         # calculate statistics
#         if restricted_model_converged and full_model_converged:
#             lr_stat, lr_pval, lr_df_diff = lr_test(restricted_model, full_model)
#         else:
#             lr_stat = np.nan
#             lr_pval = np.nan
#             lr_df_diff = 0
#
#         if restricted_model_converged:
#             rsquared_restricted = prsquared_adj(restricted_model)
#             rsquared_restricted_raw = restricted_model.prsquared
#             restricted_model_llf = restricted_model.llf
#         else:
#             rsquared_restricted = np.nan
#             rsquared_restricted_raw = np.nan
#             restricted_model_llf = np.nan
#
#         if full_model_converged:
#             rsquared = prsquared_adj(full_model)
#             rsquared_raw = full_model.prsquared
#             full_model_llf = full_model.llf
#             term_pvalues = full_model.pvalues.to_dict()
#             full_model_params = full_model.params.to_dict()
#         else:
#             rsquared = np.nan
#             rsquared_raw = np.nan
#             full_model_llf = np.nan
#             term_pvalues = {}
#             full_model_params = {}
#     else:
#         # fit restricted model
#         if fitted_restricted_model is None:
#             restricted_model_converged = False
#             restricted_model = smf.ols(
#                 formatted_restricted_formula,
#                 data = data_df
#             )
#             try:
#                 restricted_model = restricted_model.fit()
#                 restricted_model_converged = True
#             except np.linalg.LinAlgError as e:
#                 pass
#         else:
#             restricted_model = fitted_restricted_model
#             restricted_model_converged = True
#
#         # fit full model
#         full_model_converged = False
#         full_model = smf.ols(
#             formatted_full_formula,
#             data = data_df
#         )
#         try:
#             full_model = full_model.fit()
#             full_model_converged = True
#         except np.linalg.LinAlgError as e:
#             pass
#
#         # calculate statistics
#         if restricted_model_converged and full_model_converged:
#             lr_stat, lr_pval, lr_df_diff = full_model.compare_lr_test(restricted_model)
#         else:
#             lr_stat = np.nan
#             lr_pval = np.nan
#             lr_df_diff = 0
#
#         if restricted_model_converged:
#             rsquared_restricted = restricted_model.rsquared_adj
#             rsquared_restricted_raw = restricted_model.rsquared
#             restricted_model_llf = restricted_model.llf
#         else:
#             rsquared_restricted = np.nan
#             rsquared_restricted_raw = np.nan
#             restricted_model_llf = np.nan
#
#         if full_model_converged:
#             rsquared = full_model.rsquared_adj
#             rsquared_raw = full_model.rsquared
#             full_model_llf = full_model.llf
#             term_pvalues = full_model.pvalues.to_dict()
#             full_model_params = full_model.params.to_dict()
#         else:
#             rsquared = np.nan
#             rsquared_raw = np.nan
#             full_model_llf = np.nan
#             term_pvalues = {}
#             full_model_params = {}

# %%
def fit(
    formula,
    data_df,
    boolean_regression: bool=False,
):
    if boolean_regression:
        model_converged = False
        model = firth_regression(
            formula,
            data = data_df
        )
        try:
            model = model.fit()
            model_converged = model.mle_retvals["converged"]
        except np.linalg.LinAlgError as e:
            pass


#         converged = False
#         try:
#             model = smf.logit(
#                 formula,
#                 data = data_df
#             ).fit(maxiter=200, disp=0, warn_convergence=False)

#             converged = model.mle_retvals["converged"]
#         except np.linalg.LinAlgError as linalg_error:
#             # newton failed
#             pass

#         if not converged:
#             # retry with l-bfgs-b
#             print("retry with l-bfgs-b")
#             model = smf.logit(
#                 formula,
#                 data = data_df
#             ).fit(method="lbfgs", maxiter=1000, disp=0, warn_convergence=False)

        if model_converged:
            rsquared = prsquared_adj(model)
            rsquared_raw = model.prsquared
            llf = model.llf
            term_pvalues = model.pvalues.to_dict()
            params = model.params.to_dict()
        else:
            rsquared = np.nan
            rsquared_raw = np.nan
            llf = np.nan
            term_pvalues = {}
            params = {}
    else:
        # fit restricted model
        model_converged = False
        model = smf.ols(
            formula,
            data = data_df
        )
        try:
            model = model.fit()
            model_converged = True
        except np.linalg.LinAlgError as e:
            pass

        if model_converged:
            rsquared = model.rsquared_adj
            rsquared_raw = model.rsquared
            llf = model.llf
            term_pvalues = model.pvalues.to_dict()
            params = model.params.to_dict()
        else:
            rsquared = np.nan
            rsquared_raw = np.nan
            llf = np.nan
            term_pvalues = {}
            full_model_params = {}
    
    return model, {
        "converged": model_converged,
        "rsquared": rsquared,
        "rsquared_raw": rsquared_raw,
        "llf": llf,
        "term_pvalues": term_pvalues,
        "params": params,
    }

# %%
import statsmodels.formula.api as smf
import sklearn
import patsy

from threadpoolctl import threadpool_limits
from typing import List, Union

def test_association(
    pd_df,
    groupby_columns: List[str],
    full_formula: str,
    restricted_formula: str,
    covariates_df: pd.DataFrame,
    clumping_variants_df: pd.DataFrame = None,
    add_clumping: bool = config["covariates"]["add_clumping"],
    boolean_regression: bool = is_boolean_dtype,
    normalize_features=True,
    aggregate_variables=aggregate_variables,
    return_models=False
):
    keys=pd_df.loc[:, groupby_columns].iloc[0].to_dict()

    formatted_restricted_formula = format_formula(
        formula=restricted_formula,
        keys=keys,
        add_clumping=add_clumping,
        clumping_variants_df=clumping_variants_df,
    )
    formatted_full_formula = format_formula(
        formula=full_formula,
        keys=keys,
        add_clumping=add_clumping,
        clumping_variants_df=clumping_variants_df,
    )
    
    missing_terms_in_full_formula = [c for c in formatted_restricted_formula.rhs_termlist if c not in formatted_full_formula.rhs_termlist]
    assert len(missing_terms_in_full_formula) == 0, (
        f"Restricted formula is not a subset of the full formula. Missing terms: '{missing_terms_in_full_formula}'"
    )

    # get necessary columns
    restricted_variables = get_variables_from_formula(formatted_restricted_formula)
    full_variables = get_variables_from_formula(formatted_full_formula)
    # the trait(s) which we want to associate:
    target_variables = get_variables_from_formula(formatted_restricted_formula, rhs=False)
    necessary_columns = {
        *restricted_variables,
        *full_variables,
        "individual",
    }
    
    tested_terms = [
        c.name() for c in formatted_full_formula.rhs_termlist
        if c not in set(formatted_restricted_formula.rhs_termlist)
    ]
    
    assert len(tested_terms) > 0, (
        "there are no tested variables! `full_formula` needs to be a superset of `restricted_formula`!"
    )

    # subset to necessary columns
    covariates_df = covariates_df[
        [c for c in covariates_df_columns if c in necessary_columns]
    ]

    # merge with phenotype df to make sure that we have scores for all individuals with a measurement
    data_df = covariates_df.merge(pd_df, on=["individual"], how="left")
    # fill missing values
    data_df = data_df.fillna({
        c: 0 for c in pd_df.columns
    })
    if normalize_features:
        data_df = data_df.assign(**{
            col: sklearn.preprocessing.StandardScaler().fit_transform(data_df[[col]])[:, 0] for col in full_variables if col in pd_df.columns
        })
    # cast phenotype to float
    for tv in target_variables:
        data_df = data_df.assign(**{
            tv: data_df[tv].astype("float32")
        })
    
    assert necessary_columns.issubset(data_df.columns), f"Missing columns in merged dataframe: '{necessary_columns.difference(data_df.columns)}'"
    
    try:
        # fit models
        restricted_model, restricted_model_stats = fit(
            formula=formatted_restricted_formula,
            data_df=data_df,
            boolean_regression=boolean_regression,
        )
        
        full_model, full_model_stats = fit(
            formula=formatted_full_formula,
            data_df=data_df,
            boolean_regression=boolean_regression,
        )
        
        # get most associating variable
        most_associating_term = pd.DataFrame({
            "variable": tested_terms,
            "pval": [full_model_stats["term_pvalues"][v] for v in tested_terms],
            "param": [full_model_stats["params"][v] for v in tested_terms],
        }).sort_values(['pval', 'param'], ascending=[True, False]).iloc[0].to_dict()
        
        if aggregate_variables and (len(tested_terms) > 0):
            # re-fit full model with only the most associating variable
            most_associating_term_formula = patsy.ModelDesc(
                lhs_termlist=formatted_full_formula.lhs_termlist,
                rhs_termlist=[
                    t for t in formatted_full_formula.rhs_termlist if (
                        t in set(formatted_restricted_formula.rhs_termlist)
                    ) or (
                        t.name() == most_associating_term["variable"]
                    )
                ]
            )
            # most_associating_term_formula.rhs_termlist += [
            #     patsy.Term([patsy.EvalFactor(f"""Q('{most_associating_term["variable"]}')""")])
            # ]

            full_model, full_model_stats = fit(
                formula=most_associating_term_formula,
                data_df=data_df,
                boolean_regression=boolean_regression,
            )

        if not restricted_model_stats["converged"] or not full_model_stats["converged"]:
            lr_stat = np.nan
            lr_pval = np.nan
            lr_df_diff = 0
        else:
            # calculate statistics
            if boolean_regression:
                lr_stat, lr_pval, lr_df_diff = lr_test(restricted_model, full_model)
            else:
                # use the statsmodels built-in likelihood ratio test
                lr_stat, lr_pval, lr_df_diff = full_model.compare_lr_test(restricted_model)
    

    except Exception as e:
        print("---------------- error log -----------------")
        print("keys:")
        print(keys)
        print("-------------- error log end ---------------")
        raise e

    stats_df = (
        pd_df
        .loc[:, groupby_columns]
        .iloc[:1]
        .copy()
        .assign(**{
            "n_observations": [int(full_model.nobs)],
            "term_pvals": [full_model_stats["term_pvalues"]], 
            "params": [full_model_stats["params"]],
            "restricted_model_converged": [restricted_model_stats["converged"]],
            "full_model_converged": [full_model_stats["converged"]],
            "restricted_model_llf": [restricted_model_stats["llf"]],
            "full_model_llf": [full_model_stats["llf"]],
            "rsquared_restricted": [restricted_model_stats["rsquared"]],
            "rsquared_restricted_raw": [restricted_model_stats["rsquared_raw"]],
            "rsquared": [full_model_stats["rsquared"]],
            "rsquared_raw": [full_model_stats["rsquared_raw"]],
            "lr_stat": [lr_stat],
            "lr_pval": [lr_pval],
            "lr_df_diff": [lr_df_diff],
            "most_associating_term": [most_associating_term["variable"]],
        })
    )
    if return_models:
        return stats_df, restricted_model, full_model
    else:
        return stats_df

# %% [raw]
# _test_term = f"I(Q('{phenotype_col}') ** 2)"
# _test = test_association(
#     pd_df=covariates_df.select("individual").collect().to_pandas().assign(gene="ENSG00000084674"),
#     groupby_columns=["gene"],
#     full_formula=f"{restricted_formula} + {_test_term}",
#     restricted_formula=restricted_formula,
#     # covariates_df=broadcast_covariates_df,
#     covariates_df=normalized_covariates_df,
#     clumping_variants_df=clumping_variants_df,
# )
# display(_test)
# assert _test.iloc[0]["most_associating_term"] == _test_term
# del _test, _test_term

# %% [raw]
# test_association(
#     pd_df=covariates_df.select("individual").collect().to_pandas().assign(gene="ENSG00000084674"),
#     groupby_columns=["gene"],
#     restricted_formula='Asthma ~ 1\n + genetic_principal_components_f22009_0_1\n + genetic_principal_components_f22009_0_2\n + genetic_principal_components_f22009_0_3\n + genetic_principal_components_f22009_0_4\n + genetic_principal_components_f22009_0_5\n + genetic_principal_components_f22009_0_6\n + genetic_principal_components_f22009_0_7\n + genetic_principal_components_f22009_0_8\n + genetic_principal_components_f22009_0_9\n + genetic_principal_components_f22009_0_10\n + genetic_principal_components_f22009_0_11\n + genetic_principal_components_f22009_0_12\n + genetic_principal_components_f22009_0_13\n + genetic_principal_components_f22009_0_14\n + genetic_principal_components_f22009_0_15\n + genetic_principal_components_f22009_0_16\n + genetic_principal_components_f22009_0_17\n + genetic_principal_components_f22009_0_18\n + genetic_principal_components_f22009_0_19\n + genetic_principal_components_f22009_0_20\n + genetic_principal_components_f22009_0_21\n + genetic_principal_components_f22009_0_22\n + genetic_principal_components_f22009_0_23\n + genetic_principal_components_f22009_0_24\n + genetic_principal_components_f22009_0_25\n + genetic_principal_components_f22009_0_26\n + genetic_principal_components_f22009_0_27\n + genetic_principal_components_f22009_0_28\n + genetic_principal_components_f22009_0_29\n + genetic_principal_components_f22009_0_30\n + PGS001849',
#     full_formula=restricted_formula,
#     # covariates_df=broadcast_covariates_df,
#     covariates_df=covariates_df.collect().to_pandas(),
#     clumping_variants_df=clumping_variants_df,
# )

# %%
import statsmodels.formula.api as smf

from threadpoolctl import threadpool_limits
from typing import List, Union

def regression(
    dataframe: pyspark.sql.DataFrame,
    groupby_columns: Union[str, List[str]], 
    full_formula: str,
    restricted_formula: str,
    covariates_df: pd.DataFrame,
    clumping_variants_df: Union[pyspark.Broadcast, pd.DataFrame] = None,
    add_clumping: bool = config["covariates"]["add_clumping"],
    boolean_regression: bool = is_boolean_dtype,
):
    if isinstance(groupby_columns, str):
        groupby_columns = [groupby_columns]
    
    print("broadcasting...")
    
    if add_clumping:
        broadcast_clumping_variants_df = broadcast(clumping_variants_df, idempotent=True)
    else:
        broadcast_clumping_variants_df = None
    
    broadcast_covariates_df = spark_broadcast_feather(covariates_df)
    
    print("broadcasting done")
    
    def spark_fit_parallel(pd_df):
        with threadpool_limits(limits=1):
            # print("dereferencing...")
            
            # unpack the clumping variant annotation if needed
            if add_clumping:
                clumping_variants_df = deref(broadcast_clumping_variants_df)
            else:
                clumping_variants_df = None
            
            # print("dereferencing done")
            
            # Now compute all necessary columns that need to be loaded into memory.
            # Therefore we need to first get all required formula terms
            keys=pd_df.loc[:, groupby_columns].iloc[0].to_dict()
            
            formatted_restricted_formula = format_formula(
                formula=restricted_formula,
                keys=keys,
                add_clumping=add_clumping,
                clumping_variants_df=clumping_variants_df,
            )
            formatted_full_formula = format_formula(
                formula=full_formula,
                keys=keys,
                add_clumping=add_clumping,
                clumping_variants_df=clumping_variants_df,
            )
            
            # get necessary columns
            restricted_variables = get_variables_from_formula(formatted_restricted_formula)
            full_variables = get_variables_from_formula(formatted_full_formula)
            target_variables = get_variables_from_formula(formatted_restricted_formula, rhs=False)
            # Necessary are the union of all term variables and the individual
            necessary_columns = {
                *restricted_variables,
                *full_variables,
                "individual",
            }
            
            # Now load only the necessary columns into memory.
            # Make sure to use memory-mapping and disable threading.
            covariates_df = feather.read_feather(
                get_sparkfile_path(broadcast_covariates_df),
                columns=[c for c in covariates_df_columns if c in necessary_columns],
                use_threads=False,
                memory_map=True
            )
            
            # Finally, pass data to the fit function
            retval = test_association(
                pd_df,
                groupby_columns=groupby_columns,
                full_formula=full_formula,
                restricted_formula=restricted_formula,
                covariates_df=covariates_df,
                clumping_variants_df=clumping_variants_df,
                add_clumping=add_clumping,
                boolean_regression=boolean_regression,
            )

            return retval
    
    return dataframe.groupby(groupby_columns).applyInPandas(
        func=spark_fit_parallel,
        schema=t.StructType([
            *[dataframe.schema[k] for k in groupby_columns],
            t.StructField("n_observations", t.LongType()),
            t.StructField("term_pvals", t.MapType(t.StringType(), t.DoubleType())),
            t.StructField("params", t.MapType(t.StringType(), t.DoubleType())),
            t.StructField("restricted_model_converged", t.BooleanType()),
            t.StructField("full_model_converged", t.BooleanType()),
            t.StructField("restricted_model_llf", t.DoubleType()),
            t.StructField("full_model_llf", t.DoubleType()),
            t.StructField("rsquared_restricted", t.DoubleType()),
            t.StructField("rsquared_restricted_raw", t.DoubleType()),
            t.StructField("rsquared", t.DoubleType()),
            t.StructField("rsquared_raw", t.DoubleType()),
            t.StructField("lr_stat", t.DoubleType()),
            t.StructField("lr_pval", t.DoubleType()),
            t.StructField("lr_df_diff", t.DoubleType()),
            t.StructField("most_associating_term", t.StringType()),
        ])
    )

_test_term = f"I(Q('{phenotype_col}') ** 2)"
_test = regression(
    dataframe=spark.createDataFrame(covariates_df.select("individual").collect().to_pandas().assign(gene="ENSG00000084674")), 
    groupby_columns=["gene"],
    full_formula=f"{restricted_formula} + {_test_term}",
    restricted_formula=restricted_formula,
    covariates_df=normalized_covariates_df,
    clumping_variants_df=broadcast_clumping_variants_df,
).toPandas()

# %%
assert _test.iloc[0]["most_associating_term"] == _test_term
display(_test)
del _test, _test_term

# %% [markdown]
# # read features

# %%
groupby_columns = config["groupby_columns"]
groupby_columns

# %%
config["feature_sets"]

# %%
config["variables"]

# %% [markdown]
# ## protein-coding genes

# %%
protein_coding_genes_df = (
    spark.read.parquet(snakemake.input["protein_coding_genes_pq"])
    .withColumnRenamed("gene_id", "gene")
)
protein_coding_genes_df.printSchema()

# %% [markdown]
# ## feature dfs

# %%
feature_dfs = {}
for feature_name, path in config["feature_sets"].items():
    feature_dfs[feature_name] = (
        spark.read.parquet(path + "/data.parquet")
        .filter(f.col("individual").isNotNull())
        .filter(f.col("gene").isNotNull())
    )

# %%
len(feature_dfs)

# %%
from rep.data import join_featuresets

# %%
fill_values = config.get('fill_values')
if fill_values is not None and len(fill_values) > 0:
    # quote the column names with backticks
    fill_values = {"`" + k + "`": v for k, v in fill_values.items()}

features_df = join_featuresets(
    dataframes=feature_dfs,
    variables=config["variables"],
    index_cols=["individual", *groupby_columns],
    fill_values=fill_values,
    join="outer",
)
features_df = protein_coding_genes_df.join(
    features_df,
    on=["gene"],
    how="inner"
)
features_df.printSchema()

# %%
# features_df.select("gene").distinct().count()

# %%
renamed_features_df = features_df
features = []

for c in features_df.columns:
    if c.startswith("feature."):
        # rename columns because spark is stupid and fails with selecting in a pandas udf
        new_name = (
            c[8:]
            .replace(".", "_")
            # .replace("@", "__")
        )
        renamed_features_df = renamed_features_df.withColumnRenamed(c, new_name)
        
        features.append(new_name)

renamed_features_df.printSchema()

# %%
features

# %% [markdown]
# # perform regression

# %%
# scores_sdf = (
#     features_df
#     # .filter(f.col("subtissue") == "Cells - Cultured fibroblasts")
#     .filter(f.col("individual").isNotNull())
#     .select(
#         "individual",
#         *groupby_columns,
#         *[f"`{c}`" for c in features],
#     )
# )

# scores_sdf.printSchema()

# %%
# renamed_features_df = renamed_features_df.persist()
# renamed_features_df.count()

# %%
full_formula = "\n + ".join([
    restricted_formula,
    *[f"Q('{c}')" for c in features],
])
print(full_formula)

# %%
# TODO: Always use `covariates_df=normalized_covariates_df`?
regression_results_sdf = regression(
    renamed_features_df, 
    groupby_columns=groupby_columns, 
    full_formula=full_formula,
    restricted_formula=restricted_formula,
    covariates_df=normalized_covariates_df,
    clumping_variants_df=broadcast_clumping_variants_df,
)
regression_results_sdf.printSchema()

# %%
for k, v in snakemake.wildcards.items():
    regression_results_sdf = regression_results_sdf.withColumn(k, f.lit(v))

# %%
regression_results_sdf.printSchema()

# %%

# %%
# for debugging

# %% [raw]
# gene_list = renamed_features_df.select("gene").distinct().limit(10).toPandas()["gene"].tolist()
# gene_list

# %% [raw]
# gene_list = ['ENSG00000004059']
# gene_list = ['ENSG00000004059', 'ENSG00000171505']
# gene_list = [
#     'ENSG00000183706',
#     'ENSG00000185386',
#     'ENSG00000205238',
#     'ENSG00000285472',
#     'ENSG00000120280',
#     'ENSG00000176679',
# ]
# gene_list = [
#     'ENSG00000242515',
# ]

# %% [raw]
# %%time
# test = regression(
#     renamed_features_df.filter(f.col("gene").isin(
#         gene_list
#     )).repartition(2, "gene"), 
#     groupby_columns=groupby_columns, 
#     full_formula=full_formula,
#     restricted_formula=restricted_formula,
#     covariates_df=normalized_covariates_df,
#     clumping_variants_df=broadcast_clumping_variants_df,
# ).toPandas()
# test

# %% [raw]
# test_df = renamed_features_df.filter(f.col("gene").isin(
#     # ['ENSG00000171505']
#     ['ENSG00000205238']
# )).toPandas()
# test_df

# %% [raw]
# %%time
# test = test_association(
#     test_df, 
#     groupby_columns=groupby_columns, 
#     full_formula=full_formula,
#     restricted_formula=restricted_formula,
#     covariates_df=normalized_covariates_df,
#     clumping_variants_df=clumping_variants_df,
# )
# test

# %% [raw]
# normalized_test_df = test_df.assign(**{
#     col: sklearn.preprocessing.StandardScaler().fit_transform(test_df[[col]])[:, 0] for col in features
# })
# normalized_test_df

# %% [raw]
# %%time
# test = fit(
#     normalized_test_df, 
#     groupby_columns=groupby_columns, 
#     full_formula=full_formula,
#     restricted_formula=restricted_formula,
#     covariates_df=normalized_covariates_df,
#     clumping_variants_df=clumping_variants_df,
# )
# test

# %%

# %%
snakemake.output

# %%
regression_results_sdf.write.parquet(snakemake.output["associations_pq"], mode="overwrite")

# %%
# regression_results_sdf = spark.read.parquet(snakemake.output["associations_pq"])

# %% [raw]
# # read association results

# %% [raw]
# snakemake.output["associations_pq"]

# %% [raw]
# regression_results_df = (
#     spark.read.parquet(snakemake.output["associations_pq"])
#     .sort("rsquared", ascending=False)
#     # .withColumn("score_pval", f.col("term_pvals")["AbExp_DNA"])
#     .drop(
#         "term_pvals",
#         "params",
#     )
#     .toPandas()
# )
# regression_results_df

# %%
tmpdir.cleanup()


# %%
