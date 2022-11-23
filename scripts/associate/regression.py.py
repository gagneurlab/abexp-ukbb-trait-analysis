# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
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

# %% [raw] {"tags": []}
#

# %%
# we need ray to efficiently share the covariates df in memory
# use_ray = True
use_ray = False

# %% {"tags": []}
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
#del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'associate__regression',
        default_wildcards={
            "phenotype_col": "severe_LDL",
            # "phenotype_col": "Asthma",
            # "phenotype_col": "triglycerides_f30870_0_0",
            # "phenotype_col": "standing_height_f50_0_0",
            # "phenotype_col": "body_mass_index_bmi_f21001_0_0",
            # "phenotype_col": "systolic_blood_pressure_automated_reading_f4080_0_0",
            # "phenotype_col": "hdl_cholesterol_f30760_0_0",
            "feature_set": "LOFTEE_pLoF",
            # "feature_set": "AbExp_all_tissues",
            # "feature_set": "LOFTEE_pLoF",
            # "covariates": "sex_age_genPC_CLMP_PRS",
            # "covariates": "sex_age_genPC_CLMP",
            "covariates": "sex_age_genPC",
        }
    )

# %%
print(json.dumps(snakemake.__dict__, indent=2, default=str))

# %% [markdown] {"tags": []}
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

# %% [markdown]
# # Read covariates

# %%
snakemake.input["covariates_ipc"]

# %%
covariates_df = pl.scan_ipc(snakemake.input["covariates_ipc"])
covariates_df_columns = covariates_df.columns

# %% [raw]
# display(
#     "Size of 'covariates_df': %.3fGb" % (covariates_df.to_arrow().nbytes / 1024**3)
# )

# %%
covariates_sparkfilename = os.path.basename(snakemake.input["covariates_ipc"])
covariates_sparkfilename

# %%
# add file to all spark nodes
spark.sparkContext.addFile(snakemake.input["covariates_ipc"])

# %% [markdown] {"tags": []}
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
def format_formula(formula, keys, add_clumping=True, clumping_variants_df=clumping_variants_df):
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


# %%
# test_formula = format_formula(
#     formula=restricted_formula,
#     clumping_variants_df=clumping_variants_df,
#     add_clumping=True,
#     keys={
#         "gene": "ENSG00000160584",
#     }
# )

# %%
# print(test_formula)

# %%
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
def get_variables_from_formula(formula, lhs=True, rhs=True):
    if not isinstance(formula, patsy.ModelDesc):
        model_desc = patsy.ModelDesc.from_formula(formula)
    else:
        model_desc = formula

    covariate_cols = [
        *([factor.name() for term in model_desc.lhs_termlist for factor in term.factors] if lhs else []),
        *([factor.name() for term in model_desc.rhs_termlist for factor in term.factors] if rhs else []),
    ]
    # deduplicate
    covariate_cols = list(dict.fromkeys(covariate_cols))
    
    # replace Q($i) -> $i
    regex = re.compile(r"""Q\(['"](.*)['"]\)""")
    covariate_cols = [regex.sub(r"\1", c) for c in covariate_cols]
    
    return covariate_cols


# %%
# get_variables_from_formula(test_formula)

# %%
def broadcast(obj, idempotent=True):
    if idempotent:
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
def deref(obj):
    if isinstance(obj, pyspark.Broadcast):
        return obj.value
    elif ray is not None and isinstance(obj, ray.ObjectRef):
        with ray_context:
            return ray.get(obj)
    else:
        return obj


# %% {"tags": []}
# broadcast_covariates_df = broadcast(covariates_df.to_arrow())
broadcast_clumping_variants_df = broadcast(clumping_variants_df)

# %%
# covariates_pq_path = snakemake.input["covariates_pq"]
# covariates_df_columns = covariates_df.columns

# %%
from pyarrow import feather
import pyarrow as pa

# %%
data_df = feather.read_table(
    snakemake.input["covariates_ipc"],
    use_threads=False,
    memory_map=True
)

data_type = data_df.schema[data_df.schema.get_field_index(snakemake.wildcards["phenotype_col"])].type
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

# %% {"tags": []}
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
    
    df0, df1 = restricted_model.df_model, full_model.df_model
    df_diff = full_model.df_model - restricted_model.df_model
    
    chi2_stat = likelihood_ratio(
        restricted_model.llf,
        df_diff,
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


# %% {"tags": []}
import statsmodels.formula.api as smf

from threadpoolctl import threadpool_limits
from typing import List, Union

def regression(
    dataframe: pyspark.sql.DataFrame,
    groupby_columns: Union[str, List[str]], 
    full_formula: str,
    restricted_formula: str,
    # covariates_df: Union[pyspark.Broadcast, pd.DataFrame],
    clumping_variants_df: Union[pyspark.Broadcast, pd.DataFrame] = None,
    add_clumping: bool = config["covariates"]["add_clumping"],
    boolean_regression=is_boolean_dtype,
):
    if isinstance(groupby_columns, str):
        groupby_columns = [groupby_columns]
    
    print("broadcasting...")
    # make sure to use broadcasted variable
    # broadcast_covariates_df = broadcast(covariates_df, idempotent=True)
    
    if add_clumping:
        broadcast_clumping_variants_df = broadcast(clumping_variants_df, idempotent=True)
    else:
        broadcast_clumping_variants_df = None
    
    print("broadcasting done")
    
    def fit(pd_df):
        with threadpool_limits(limits=1):
            keys=pd_df.loc[:, groupby_columns].iloc[0].to_dict()
    
            # print("dereferencing...")
            
            # ## unpack all broadcast variables
            # covariates_df = deref(broadcast_covariates_df)
            
            # unpack only if needed
            if add_clumping:
                clumping_variants_df = deref(broadcast_clumping_variants_df)
            else:
                clumping_variants_df = None
            # ## unpacking done
            
            # print("dereferencing done")
            
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
            necessary_columns = {
                *restricted_variables,
                *full_variables,
                "individual",
            }
            
            # make sure to use memory-mapping and disable threading
            covariates_df = feather.read_feather(
                pyspark.SparkFiles.get(covariates_sparkfilename),
                columns=[c for c in covariates_df_columns if c in necessary_columns],
                use_threads=False,
                memory_map=True
            )
            
            # merge with phenotype df to make sure that we have all scores predicted
            data_df = covariates_df.merge(pd_df, on=["individual"], how="left")
            # fill missing values
            data_df = data_df.fillna({
                c: 0 for c in pd_df.columns
            })
            # cast phenotype to float
            for tv in target_variables:
                data_df = data_df.assign(**{
                    tv: data_df[tv].astype("float32")
                })
            
            
            if boolean_regression:
                try:
                    converged = False
                    try:
                        restricted_model = smf.logit(
                            formatted_restricted_formula,
                            data = data_df
                        ).fit(maxiter=200, disp=0, warn_convergence=False)

                        full_model = smf.logit(
                            formatted_full_formula,
                            data = data_df
                        ).fit(maxiter=200, disp=0, warn_convergence=False)
                        
                        converged = restricted_model.mle_retvals["converged"] & full_model.mle_retvals["converged"]
                    except np.linalg.LinAlgError as linalg_error:
                        # newton failed
                        pass
                    
                    if not converged:
                        # retry with l-bfgs-b
                        # print("retry with l-bfgs-b")
                        restricted_model = smf.logit(
                            formatted_restricted_formula,
                            data = data_df
                        ).fit(method="lbfgs", maxiter=1000, disp=0, warn_convergence=False)

                        full_model = smf.logit(
                            formatted_full_formula,
                            data = data_df
                        ).fit(method="lbfgs", maxiter=1000, disp=0, warn_convergence=False)
                except Exception as e:
                    print("---------------- error log -----------------")
                    print("keys:")
                    print(keys)
                    print("-------------- error log end ---------------")
                    raise e

                # calculate statistics
                lr_stat, lr_pval, lr_df_diff = lr_test(restricted_model, full_model)
                restricted_model_converged = restricted_model.mle_retvals["converged"]
                rsquared_restricted = prsquared_adj(restricted_model)
                rsquared_restricted_raw = restricted_model.prsquared
                full_model_converged = full_model.mle_retvals["converged"]
                rsquared = prsquared_adj(full_model)
                rsquared_raw = full_model.prsquared
            else:
                try:
                    converged = False
                    try:
                        restricted_model = smf.ols(
                            formatted_restricted_formula,
                            data = data_df
                        ).fit(maxiter=200, disp=0, warn_convergence=False)

                        full_model = smf.ols(
                            formatted_full_formula,
                            data = data_df
                        ).fit(maxiter=200, disp=0, warn_convergence=False)
                        
                        converged = restricted_model.mle_retvals["converged"] & full_model.mle_retvals["converged"]
                    except np.linalg.LinAlgError as linalg_error:
                        # newton failed
                        pass
                    
                    if not converged:
                        # retry with l-bfgs-b
                        # print("retry with l-bfgs-b")
                        restricted_model = smf.ols(
                            formatted_restricted_formula,
                            data = data_df
                        ).fit(method="lbfgs", maxiter=1000, disp=0, warn_convergence=False)

                        full_model = smf.ols(
                            formatted_full_formula,
                            data = data_df
                        ).fit(method="lbfgs", maxiter=1000, disp=0, warn_convergence=False)
                except Exception as e:
                    print("---------------- error log -----------------")
                    print("keys:")
                    print(keys)
                    print("-------------- error log end ---------------")
                    raise e

                # calculate statistics
                lr_stat, lr_pval, lr_df_diff = full_model.compare_lr_test(restricted_model)
                restricted_model_converged = restricted_model.mle_retvals["converged"]
                rsquared_restricted = restricted_model.rsquared_adj
                rsquared_restricted_raw = restricted_model.rsquared
                full_model_converged = full_model.mle_retvals["converged"]
                rsquared = full_model.rsquared_adj
                rsquared_raw = full_model.rsquared

            return (
                pd_df
                .loc[:, groupby_columns]
                .iloc[:1]
                .copy()
                .assign(**{
                    "n_observations": [int(full_model.nobs)],
                    "term_pvals": [full_model.pvalues.to_dict()], 
                    "params": [full_model.params.to_dict()],
                    "restricted_model_converged": [restricted_model_converged],
                    "full_model_converged": [full_model_converged],
                    "restricted_model_llf": [restricted_model.llf],
                    "full_model_llf": [full_model.llf],
                    "rsquared_restricted": [rsquared_restricted],
                    "rsquared_restricted_raw": [rsquared_restricted_raw],
                    "rsquared": [rsquared],
                    "rsquared_raw": [rsquared_raw],
                    "lr_stat": [lr_stat],
                    "lr_pval": [lr_pval],
                    "lr_df_diff": [lr_df_diff],
                })
            )
    
    return dataframe.groupby(groupby_columns).applyInPandas(
        func=fit,
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
        ])
    )

_test = regression(
    dataframe=spark.createDataFrame(covariates_df.select("individual").collect().to_pandas().assign(gene="ENSG00000084674")), 
    groupby_columns=["gene"],
    full_formula=restricted_formula,
    restricted_formula=restricted_formula,
    # covariates_df=broadcast_covariates_df,
    clumping_variants_df=broadcast_clumping_variants_df,
).toPandas()

# %%
assert _test.iloc[0]["lr_df_diff"] == 0
display(_test)
del _test

# %% [markdown] {"tags": []}
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
regression_results_sdf = regression(
    renamed_features_df, 
    groupby_columns=groupby_columns, 
    full_formula=full_formula,
    restricted_formula=restricted_formula,
    # covariates_df=broadcast_covariates_df,
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

# %% [raw]
# %%time
# test = regression(
#     renamed_features_df.filter(f.col("gene").isin(
#         gene_list
#     )).repartition(2, "gene"), 
#     groupby_columns=groupby_columns, 
#     full_formula=full_formula,
#     restricted_formula=restricted_formula,
#     # covariates_df=broadcast_covariates_df,
#     clumping_variants_df=broadcast_clumping_variants_df,
# ).toPandas()
# test

# %%

# %%
snakemake.output

# %% {"tags": []}
regression_results_sdf.write.parquet(snakemake.output["associations_pq"], mode="overwrite")

# %%
# regression_results_sdf = spark.read.parquet(snakemake.output["associations_pq"])

# %% [raw] {"tags": []}
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

