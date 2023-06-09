import utils
import warnings
import pandas as pd
from scipy import stats
import normalization_analysis
import imputation_analysis



warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.filterwarnings("ignore", category=RuntimeWarning) 


def transform_proteomics_edgelist(df, index_cols=['group', 'sample', 'subject'], drop_cols=['sample'], group='group', identifier='identifier', extra_identifier='name', value_col='LFQ_intensity'):
    """
    Transforms a long format proteomics matrix into a wide format

    :param df: long-format pandas dataframe with columns 'group', 'sample', 'subject', 'identifier' (protein), 'name' (gene) and 'LFQ_intensity'.
    :param list index_cols: column labels to be be kept as index identifiers.
    :param list drop_cols: column labels to be dropped from the dataframe.
    :param str group: column label containing group identifiers.
    :param str identifier: column label containing feature identifiers.
    :param str extra_identifier: column label containing additional protein identifiers (e.g. gene names).
    :param str value_col: column label containing expression values.

    :return: Pandas dataframe with samples as rows and protein identifiers (UniprotID~GeneName) as columns (with additional columns 'group', 'sample' and 'subject').

    Example:
        df = transform_proteomics_edgelist(original, index_cols=['group', 'sample', 'subject'], drop_cols=['sample'], group='group', identifier='identifier', value_col='LFQ_intensity')
    """
    wdf = None
    if df.columns.isin(index_cols).sum() == len(index_cols):
        wdf = df[df[group].notna()]
        wdf[index_cols] = wdf[index_cols].astype(str)
        wdf = wdf.set_index(index_cols)
        if extra_identifier is not None and extra_identifier in wdf.columns:
            wdf[identifier] = wdf[extra_identifier].map(str) + "~" + wdf[identifier].map(str)
        wdf = wdf.pivot_table(values=value_col, index=wdf.index, columns=identifier, aggfunc='first')
        wdf = wdf.reset_index()
        wdf[index_cols] = wdf["index"].apply(pd.Series)
        wdf = wdf.drop(["index"], axis=1)

    return wdf


def get_proteomics_measurements_ready(df, index_cols=['group', 'sample', 'subject'], drop_cols=['sample'], group='group', identifier='identifier', extra_identifier='name', filter_samples=False, filter_samples_percent=0.5, imputation=True, imputation_method='distribution', missing_method='percentage', missing_per_group=True, missing_max=0.3, min_valid=1, value_col='LFQ_intensity', shift=1.8, nstd=0.3, knn_cutoff=0.6, normalize=False, normalization_method='median', normalize_group=False, normalize_by=None):
    """
    Processes proteomics data extracted from the database: 1) filter proteins with high number of missing values (> missing_max or min_valid), 2) impute missing values.

    :param df: long-format pandas dataframe with columns 'group', 'sample', 'subject', 'identifier' (protein), 'name' (gene) and 'LFQ_intensity'.
    :param list index_cols: column labels to be be kept as index identifiers.
    :param list drop_cols: column labels to be dropped from the dataframe.
    :param str group: column label containing group identifiers.
    :param str identifier: column label containing feature identifiers.
    :param str extra_identifier: column label containing additional protein identifiers (e.g. gene names).
    :param bool filter_samples: if True filter samples with valid values below percentage (filter_samples_percent).
    :param float filter_samples_percent: defines the maximum percentage of missing values allowed in a sample.
    :param bool imputation: if True performs imputation of missing values.
    :param str imputation_method:  method for missing values imputation ('KNN', 'distribuition', or 'mixed')
    :param str missing_method: defines which expression rows are counted to determine if a column has enough valid values to survive the filtering process.
    :param bool missing_per_group: if True filter proteins based on valid values per group; if False filter across all samples.
    :param float missing_max: maximum ratio of missing/valid values to be filtered.
    :param int min_valid: minimum number of valid values to be filtered.
    :param str value_col: column label containing expression values.
    :param float shift: when using distribution imputation, the down-shift
    :param float nstd: when using distribution imputation, the width of the distribution
    :param float knn_cutoff: when using KNN imputation, the minimum percentage of valid values for which to use KNN imputation (i.e. 0.6 -> if 60% valid values use KNN, otherwise MinProb)
    :param bool normalize: whether or not to normalize the data
    :param str normalization_method: method to be used to normalize the data ('median', 'quantile', 'linear', 'zscore', 'median_polish') (only with normalize=True)
    :param bool normalize_group: normalize per group or not (only with normalize=True)
    :param str normalize_by: whether the normalization should be done by 'features' (columns) or 'samples' (rows) (only with normalize=True)
    :return: Pandas dataframe with samples as rows and protein identifiers as columns (with additional columns 'group', 'sample' and 'subject').

    Example 1::

        result = get_proteomics_measurements_ready(df, index_cols=['group', 'sample', 'subject'], drop_cols=['sample'], group='group', identifier='identifier', extra_identifier='name', imputation=True, method = 'distribution', missing_method = 'percentage', missing_per_group=True, missing_max = 0.3, value_col='LFQ_intensity')

    Example 2::

        result = get_proteomics_measurements_ready(df, index_cols=['group', 'sample', 'subject'], drop_cols=['sample'], group='group', identifier='identifier', extra_identifier='name', imputation = True, method = 'mixed', missing_method = 'at_least_x', missing_per_group=False, min_valid=5, value_col='LFQ_intensity')
    """
    df = transform_proteomics_edgelist(df, index_cols=index_cols, drop_cols=drop_cols, group=group, identifier=identifier, extra_identifier=extra_identifier, value_col=value_col)
    if df is not None:
        aux = []
        aux.extend(index_cols)
        g = group
        if not missing_per_group:
            g = None
        if missing_method == 'at_least_x':
            aux.extend(imputation_analysis.extract_number_missing(df, min_valid, drop_cols, group=g))
        elif missing_method == 'percentage':
            aux.extend(imputation_analysis.extract_percentage_missing(df,  missing_max, drop_cols, group=g))

        df = df[list(set(aux))]
        if filter_samples:
            df = df.loc[~(df.T.isna().mean() > filter_samples_percent)]

        if normalize:
            if not normalize_group:
                df = normalization_analysis.normalize_data(df, method=normalization_method, normalize=normalize_by)
            else:
                df = normalization_analysis.normalize_data_per_group(df, group=group, method=normalization_method, normalize=normalize_by)

        if imputation:
            if imputation_method.lower() == "knn":
                df = imputation_analysis.imputation_KNN(df, drop_cols=index_cols, group=group, cutoff=knn_cutoff, alone=True)
            elif imputation_method == "distribution":
                df = imputation_analysis.imputation_normal_distribution(df, index_cols=index_cols, shift=shift, nstd=nstd)
                df = df.reset_index()
            elif imputation_method == 'mixed':
                df = imputation_analysis.imputation_mixed_norm_KNN(df, index_cols=index_cols, shift=shift, nstd=nstd, group=group, cutoff=knn_cutoff)
                df = df.reset_index()
        else:
            df = df.set_index(index_cols).dropna(axis=1).reset_index()
    return df


def get_clinical_measurements_ready(df, subject_id='subject', sample_id='biological_sample', group_id='group', columns=['clinical_variable'], values='values', extra=['group'], imputation=True, imputation_method='KNN', missing_method='percentage', missing_max=0.3, min_valid=1):
    """
    Processes clinical data extracted from the database by converting dataframe to wide-format and imputing missing values.

    :param df: long-format pandas dataframe with columns 'group', 'biological_sample', 'subject', 'clinical_variable', 'value'.
    :param str subject_id: column label containing subject identifiers.
    :param str sample_id: column label containing biological sample identifiers.
    :param str group_id: column label containing group identifiers.
    :param list columns: column name whose unique values will become the new column names
    :param str values: column label containing clinical variable values.
    :param list extra: additional column labels to be kept as columns
    :param bool imputation: if True performs imputation of missing values.
    :param str imputation_method: method for missing values imputation ('KNN', 'distribuition', or 'mixed').
    :param str missing_method: defines which expression rows are counted to determine if a column has enough valid values to survive the filtering process.
    :param float missing_max: maximum ratio of missing/valid values to be filtered.
    :param int min_valid: minimum number of valid values to be filtered.
    :return: Pandas dataframe with samples as rows and clinical variables as columns (with additional columns 'group', 'subject' and 'biological_sample').

    Example::

        result = get_clinical_measurements_ready(df, subject_id='subject', sample_id='biological_sample', group_id='group', columns=['clinical_variable'], values='values', extra=['group'], imputation=True, imputation_method='KNN')
    """
    index = [subject_id, sample_id]
    aux = [group_id]
    aux.extend(index)
    drop_cols = index

    processed_df = df[df['rel_type'] == 'HAS_QUANTIFIED_CLINICAL']
    processed_df[values] = processed_df[values].astype('float')
    processed_df = utils.transform_into_wide_format(processed_df, index=index, columns=columns, values=values, extra=extra)
    if missing_method == 'at_least_x':
        aux.extend(imputation_analysis.extract_number_missing(processed_df, min_valid, drop_cols, group=group_id))
    elif missing_method == 'percentage':
        aux.extend(imputation_analysis.extract_percentage_missing(processed_df,  missing_max, drop_cols, group=group_id))

    processed_df = processed_df[list(set(aux))]
    drop_cols.append(group_id)
    if imputation:
        if imputation_method.lower() == "knn":
            df = imputation_analysis.imputation_KNN(processed_df, drop_cols=drop_cols, group=group_id, cutoff=0.0)
        elif imputation_method.lower() == "distribution":
            df = imputation_analysis.imputation_normal_distribution(processed_df, index_cols=index)
        elif imputation_method.lower() == 'mixed':
            df = imputation_analysis.imputation_mixed_norm_KNN(processed_df,index_cols=index, group=group_id)

    return df


def get_ranking_with_markers(data, drop_columns, group, columns, list_markers, annotation={}):
    """
    This function creates a long-format dataframe with features and values to be plotted together with disease biomarker annotations.

    :param data: wide-format Pandas DataFrame with samples as rows and features as columns
    :param list drop_columns: columns to be deleted
    :param str group: column to use as identifier variables
    :param list columns: names to use for the 1)variable column, and for the 2)value column
    :param list list_markers: list of features from data, known to be markers associated to disease.
    :param dict annotation: markers, from list_markers, and associated diseases.
    :return: Long-format pandas DataFrame with group identifiers as rows and columns: 'name' (identifier), 'y' (LFQ intensity), 'symbol' and 'size'.

    Example::

        result = get_ranking_with_markers(data, drop_columns=['sample', 'subject'], group='group', columns=['name', 'y'], list_markers, annotation={})
    """
    long_data = pd.DataFrame()
    if data is not None:
        long_data = utils.transform_into_long_format(data, drop_columns, group, columns)
        if len(set(long_data['name'].values.tolist()).intersection(list_markers)) > 0:
            long_data = long_data.drop_duplicates()
            long_data['symbol'] = [17 if p in list_markers else 0 for p in long_data['name'].tolist()]
            long_data['size'] = [25 if p in list_markers else 7 for p in long_data['name'].tolist()]
            long_data['name'] = [p+' marker in '+annotation[p] if p in annotation else p for p in long_data['name'].tolist()]

    return long_data

def get_summary_data_matrix(data):
    """
    Returns some statistics on the data matrix provided.

    :param data: pandas dataframe.
    :return: dictionary with the type of statistics as key and the statistic as value in the shape of a pandas data frame

    Example::

        result = get_summary_data_matrix(data)
    """
    summary = {}
    summary["Data Matrix Shape"] = pd.DataFrame(data=[data.shape], columns=["Rows", "Columns"])
    summary["Stats"] = data.describe().transpose().reset_index()

    return summary


def check_equal_variances(data, drop_cols=['group','sample', 'subject'], group_col='group', alpha=0.05):
    levene_results = []
    for c in data.columns.drop(drop_cols):
        values = []
        for group, g in data.groupby(group_col):
            values.append(g[c].values)
        levene_results.append(stats.levene(*values)+(c,))
    levene_results = pd.DataFrame(levene_results, columns=['test', 'pvalue','identifier'])
    levene_results['pass'] = levene_results['pvalue'] > alpha

    return levene_results


def check_normality(data, drop_cols=['group','sample', 'subject'], group_col='group', alpha=0.05):
    shapiro_wilk_results = []
    for group in data[group_col].unique().tolist():
        for c in data.columns.drop(drop_cols):
            shapiro_wilk_results.append(stats.shapiro(data.loc[data[group_col] == group, c].values)+(group, c))
    shapiro_wilk_results = pd.DataFrame(shapiro_wilk_results, columns=['test', 'pvalue', 'group', 'identifier'])
    shapiro_wilk_results['pass'] = shapiro_wilk_results['pvalue'] > alpha

    return shapiro_wilk_results


def merge_for_polar(regulation_data, regulators, identifier_col='identifier', group_col='group', theta_col='modifier', aggr_func='mean', normalize=True):
    aggr_df = pd.DataFrame()
    if normalize:
        regulation_data = normalization_analysis.normalize_data(regulation_data, method='zscore')

    df = regulation_data.groupby(group_col)
    list_cols = []
    for i, group in df:
        if aggr_func == 'mean':
            list_cols.append(group.mean())
        elif aggr_func == 'median':
            list_cols.append(group.median())
        elif aggr_func == 'sum':
            list_cols.append(group.sum())
        else:
            break

    if len(list_cols) > 0:
        aggr_df = pd.DataFrame(list_cols, index=regulation_data[group_col].unique(), columns=regulation_data.set_index(group_col).columns).stack().reset_index()
        aggr_df.columns = [group_col, identifier_col, 'value']
        aggr_df = pd.merge(aggr_df, regulators, on=identifier_col)

    if not aggr_df.empty:
        list_cols = []
        for i, group in aggr_df.groupby([group_col, theta_col]):
            if aggr_func == 'mean':
                value = group['value'].mean()
            elif aggr_func == 'median':
                value = group['value'].median()
            elif aggr_func == 'sum':
                value = group['value'].sum()
            else:
                break

            list_cols.append((*i, value))

        if len(list_cols) > 0:
            aggr_df = pd.DataFrame(list_cols, columns=[group_col, theta_col, 'value'])

    return aggr_df

def run_qc_markers_analysis(data, qc_markers, sample_col='sample', group_col='group', drop_cols=['subject'], identifier_col='identifier', qcidentifier_col='identifier', qcclass_col='class'):
    """
    """
    bdf = None
    if data is not None and qc_markers is not None:
        if not data.empty and not qc_markers.empty:
            cols = list(set(data.columns.tolist()).intersection(qc_markers[qcidentifier_col].tolist()))

            if len(cols) > 0:
                nd = data.set_index([sample_col]).drop(drop_cols + [group_col], axis=1)
                z = pd.DataFrame(normalization_analysis.normalize_data(nd, method='zscore'), index=nd.index, columns=nd.columns)
                nd = z.unstack()
                nd = nd.reset_index().set_index(sample_col).join(data[[sample_col, group_col]].set_index(sample_col)).reset_index()
                bdf = pd.DataFrame()
                for i, group in qc_markers.groupby(qcclass_col):
                    c = group[group[qcidentifier_col].isin(cols)][qcidentifier_col].tolist()
                    if bdf.empty:
                        bdf = nd.set_index(identifier_col).loc[c].reset_index()
                        bdf.columns = [identifier_col, sample_col, 'z-score', group_col]
                        bdf[qcclass_col] = i
                    else:
                        aux = nd.set_index(identifier_col).loc[c].reset_index()
                        aux.columns = [identifier_col, sample_col, 'z-score', group_col]
                        aux[qcclass_col] = i
                        bdf = bdf.append(aux)

    return bdf


