# acore package

## Subpackages

* [acore.batch_correction package](acore.batch_correction.md)
  * [`combat_batch_correction()`](acore.batch_correction.md#acore.batch_correction.combat_batch_correction)
* [acore.correlation_analysis package](acore.correlation_analysis.md)
  * [`CorrelationCoefficient`](acore.correlation_analysis.md#acore.correlation_analysis.CorrelationCoefficient)
    * [`CorrelationCoefficient.coefficient`](acore.correlation_analysis.md#acore.correlation_analysis.CorrelationCoefficient.coefficient)
    * [`CorrelationCoefficient.pvalue`](acore.correlation_analysis.md#acore.correlation_analysis.CorrelationCoefficient.pvalue)
    * [`CorrelationCoefficient.count()`](acore.correlation_analysis.md#acore.correlation_analysis.CorrelationCoefficient.count)
    * [`CorrelationCoefficient.index()`](acore.correlation_analysis.md#acore.correlation_analysis.CorrelationCoefficient.index)
  * [`corr_lower_triangle()`](acore.correlation_analysis.md#acore.correlation_analysis.corr_lower_triangle)
  * [`calculate_correlations()`](acore.correlation_analysis.md#acore.correlation_analysis.calculate_correlations)
  * [`run_correlation()`](acore.correlation_analysis.md#acore.correlation_analysis.run_correlation)
  * [`run_multi_correlation()`](acore.correlation_analysis.md#acore.correlation_analysis.run_multi_correlation)
  * [`calculate_rm_correlation()`](acore.correlation_analysis.md#acore.correlation_analysis.calculate_rm_correlation)
  * [`run_rm_correlation()`](acore.correlation_analysis.md#acore.correlation_analysis.run_rm_correlation)
  * [`calculate_pvalue_correlation_old()`](acore.correlation_analysis.md#acore.correlation_analysis.calculate_pvalue_correlation_old)
  * [`calculate_pvalue_correlation()`](acore.correlation_analysis.md#acore.correlation_analysis.calculate_pvalue_correlation)
  * [`run_efficient_correlation()`](acore.correlation_analysis.md#acore.correlation_analysis.run_efficient_correlation)
* [acore.decomposition package](acore.decomposition.md)
  * [Submodules](acore.decomposition.md#submodules)
  * [acore.decomposition.pca module](acore.decomposition.md#module-acore.decomposition.pca)
    * [`run_pca()`](acore.decomposition.md#acore.decomposition.pca.run_pca)
  * [acore.decomposition.umap module](acore.decomposition.md#module-acore.decomposition.umap)
    * [`run_umap()`](acore.decomposition.md#acore.decomposition.umap.run_umap)
* [acore.differential_regulation package](acore.differential_regulation.md)
  * [`run_anova()`](acore.differential_regulation.md#acore.differential_regulation.run_anova)
  * [`run_ancova()`](acore.differential_regulation.md#acore.differential_regulation.run_ancova)
  * [`run_diff_analysis()`](acore.differential_regulation.md#acore.differential_regulation.run_diff_analysis)
  * [`run_mixed_anova()`](acore.differential_regulation.md#acore.differential_regulation.run_mixed_anova)
  * [`run_repeated_measurements_anova()`](acore.differential_regulation.md#acore.differential_regulation.run_repeated_measurements_anova)
  * [`run_ttest()`](acore.differential_regulation.md#acore.differential_regulation.run_ttest)
  * [`run_two_way_anova()`](acore.differential_regulation.md#acore.differential_regulation.run_two_way_anova)
  * [Submodules](acore.differential_regulation.md#submodules)
  * [acore.differential_regulation.format_test_table module](acore.differential_regulation.md#module-acore.differential_regulation.format_test_table)
  * [acore.differential_regulation.tests module](acore.differential_regulation.md#module-acore.differential_regulation.tests)
    * [`calc_means_between_groups()`](acore.differential_regulation.md#acore.differential_regulation.tests.calc_means_between_groups)
    * [`calc_ttest()`](acore.differential_regulation.md#acore.differential_regulation.tests.calc_ttest)
    * [`calculate_ttest()`](acore.differential_regulation.md#acore.differential_regulation.tests.calculate_ttest)
    * [`calculate_thsd()`](acore.differential_regulation.md#acore.differential_regulation.tests.calculate_thsd)
    * [`calculate_pairwise_ttest()`](acore.differential_regulation.md#acore.differential_regulation.tests.calculate_pairwise_ttest)
    * [`complement_posthoc()`](acore.differential_regulation.md#acore.differential_regulation.tests.complement_posthoc)
    * [`calculate_anova()`](acore.differential_regulation.md#acore.differential_regulation.tests.calculate_anova)
    * [`calculate_ancova()`](acore.differential_regulation.md#acore.differential_regulation.tests.calculate_ancova)
    * [`calculate_repeated_measures_anova()`](acore.differential_regulation.md#acore.differential_regulation.tests.calculate_repeated_measures_anova)
    * [`calculate_mixed_anova()`](acore.differential_regulation.md#acore.differential_regulation.tests.calculate_mixed_anova)
    * [`pairwise_ttest_with_covariates()`](acore.differential_regulation.md#acore.differential_regulation.tests.pairwise_ttest_with_covariates)
    * [`format_anova_table()`](acore.differential_regulation.md#acore.differential_regulation.tests.format_anova_table)
    * [`calculate_pvalue_from_tstats()`](acore.differential_regulation.md#acore.differential_regulation.tests.calculate_pvalue_from_tstats)
    * [`eta_squared()`](acore.differential_regulation.md#acore.differential_regulation.tests.eta_squared)
    * [`omega_squared()`](acore.differential_regulation.md#acore.differential_regulation.tests.omega_squared)
* [acore.drift_correction package](acore.drift_correction.md)
  * [`check_missingness()`](acore.drift_correction.md#acore.drift_correction.check_missingness)
  * [`cpca_centroid()`](acore.drift_correction.md#acore.drift_correction.cpca_centroid)
  * [`run_cpca_drift_correction()`](acore.drift_correction.md#acore.drift_correction.run_cpca_drift_correction)
  * [`run_loess_drift_correction()`](acore.drift_correction.md#acore.drift_correction.run_loess_drift_correction)
  * [`qc_rlsc_loess()`](acore.drift_correction.md#acore.drift_correction.qc_rlsc_loess)
  * [Submodules](acore.drift_correction.md#submodules)
  * [acore.drift_correction.cpca_drift_correction module](acore.drift_correction.md#module-acore.drift_correction.cpca_drift_correction)
    * [`check_missingness()`](acore.drift_correction.md#acore.drift_correction.cpca_drift_correction.check_missingness)
    * [`run_cpca_drift_correction()`](acore.drift_correction.md#acore.drift_correction.cpca_drift_correction.run_cpca_drift_correction)
    * [`cpca_centroid()`](acore.drift_correction.md#acore.drift_correction.cpca_drift_correction.cpca_centroid)
  * [acore.drift_correction.loess_drift_correction module](acore.drift_correction.md#module-acore.drift_correction.loess_drift_correction)
    * [`filter_features_by_qc()`](acore.drift_correction.md#acore.drift_correction.loess_drift_correction.filter_features_by_qc)
    * [`qc_rlsc_loess()`](acore.drift_correction.md#acore.drift_correction.loess_drift_correction.qc_rlsc_loess)
    * [`run_loess_drift_correction()`](acore.drift_correction.md#acore.drift_correction.loess_drift_correction.run_loess_drift_correction)
* [acore.enrichment_analysis package](acore.enrichment_analysis.md)
  * [`run_site_regulation_enrichment()`](acore.enrichment_analysis.md#acore.enrichment_analysis.run_site_regulation_enrichment)
  * [`run_up_down_regulation_enrichment()`](acore.enrichment_analysis.md#acore.enrichment_analysis.run_up_down_regulation_enrichment)
  * [`run_fisher()`](acore.enrichment_analysis.md#acore.enrichment_analysis.run_fisher)
  * [`run_kolmogorov_smirnov()`](acore.enrichment_analysis.md#acore.enrichment_analysis.run_kolmogorov_smirnov)
  * [Subpackages](acore.enrichment_analysis.md#subpackages)
    * [acore.enrichment_analysis.statistical_tests namespace](acore.enrichment_analysis.statistical_tests.md)
      * [Submodules](acore.enrichment_analysis.statistical_tests.md#submodules)
      * [acore.enrichment_analysis.statistical_tests.fisher module](acore.enrichment_analysis.statistical_tests.md#module-acore.enrichment_analysis.statistical_tests.fisher)
      * [acore.enrichment_analysis.statistical_tests.kolmogorov_smirnov module](acore.enrichment_analysis.statistical_tests.md#module-acore.enrichment_analysis.statistical_tests.kolmogorov_smirnov)
  * [Submodules](acore.enrichment_analysis.md#submodules)
  * [acore.enrichment_analysis.annotate module](acore.enrichment_analysis.md#module-acore.enrichment_analysis.annotate)
    * [`annotate_features()`](acore.enrichment_analysis.md#acore.enrichment_analysis.annotate.annotate_features)
* [acore.exploratory_analysis package](acore.exploratory_analysis.md)
  * [`get_histogram_series()`](acore.exploratory_analysis.md#acore.exploratory_analysis.get_histogram_series)
  * [`calculate_coefficient_variation()`](acore.exploratory_analysis.md#acore.exploratory_analysis.calculate_coefficient_variation)
  * [`calculate_coef_of_var_and_mean()`](acore.exploratory_analysis.md#acore.exploratory_analysis.calculate_coef_of_var_and_mean)
  * [`get_coefficient_variation()`](acore.exploratory_analysis.md#acore.exploratory_analysis.get_coefficient_variation)
  * [`extract_number_missing()`](acore.exploratory_analysis.md#acore.exploratory_analysis.extract_number_missing)
  * [`extract_percentage_missing()`](acore.exploratory_analysis.md#acore.exploratory_analysis.extract_percentage_missing)
  * [`run_pca()`](acore.exploratory_analysis.md#acore.exploratory_analysis.run_pca)
  * [`run_tsne()`](acore.exploratory_analysis.md#acore.exploratory_analysis.run_tsne)
  * [`run_umap()`](acore.exploratory_analysis.md#acore.exploratory_analysis.run_umap)
* [acore.filter_metabolomics package](acore.filter_metabolomics.md)
  * [`filter_mz_rt()`](acore.filter_metabolomics.md#acore.filter_metabolomics.filter_mz_rt)
  * [Submodules](acore.filter_metabolomics.md#submodules)
  * [acore.filter_metabolomics.filter_data module](acore.filter_metabolomics.md#module-acore.filter_metabolomics.filter_data)
    * [`filter_biological_relevance()`](acore.filter_metabolomics.md#acore.filter_metabolomics.filter_data.filter_biological_relevance)
  * [acore.filter_metabolomics.make_numeric module](acore.filter_metabolomics.md#module-acore.filter_metabolomics.make_numeric)
    * [`parse_average()`](acore.filter_metabolomics.md#acore.filter_metabolomics.make_numeric.parse_average)
    * [`convert_to_numeric()`](acore.filter_metabolomics.md#acore.filter_metabolomics.make_numeric.convert_to_numeric)
* [acore.imputation_analysis package](acore.imputation_analysis.md)
  * [`imputation_KNN()`](acore.imputation_analysis.md#acore.imputation_analysis.imputation_KNN)
  * [`imputation_mixed_norm_KNN()`](acore.imputation_analysis.md#acore.imputation_analysis.imputation_mixed_norm_KNN)
  * [`imputation_normal_distribution()`](acore.imputation_analysis.md#acore.imputation_analysis.imputation_normal_distribution)
* [acore.io package](acore.io.md)
  * [`download_PRIDE_data()`](acore.io.md#acore.io.download_PRIDE_data)
  * [`unrar()`](acore.io.md#acore.io.unrar)
  * [`download_file()`](acore.io.md#acore.io.download_file)
  * [Subpackages](acore.io.md#subpackages)
    * [acore.io.uniprot package](acore.io.uniprot.md)
      * [`fetch_annotations()`](acore.io.uniprot.md#acore.io.uniprot.fetch_annotations)
      * [`process_annotations()`](acore.io.uniprot.md#acore.io.uniprot.process_annotations)
      * [Submodules](acore.io.uniprot.md#submodules)
      * [acore.io.uniprot.uniprot module](acore.io.uniprot.md#module-acore.io.uniprot.uniprot)
  * [Submodules](acore.io.md#submodules)
  * [acore.io.ftp module](acore.io.md#module-acore.io.ftp)
    * [`download_from_ftp()`](acore.io.md#acore.io.ftp.download_from_ftp)
  * [acore.io.http module](acore.io.md#module-acore.io.http)
    * [`download_file()`](acore.io.md#acore.io.http.download_file)
  * [acore.io.pride module](acore.io.md#module-acore.io.pride)
    * [`download_PRIDE_data()`](acore.io.md#acore.io.pride.download_PRIDE_data)
  * [acore.io.uncompress module](acore.io.md#module-acore.io.uncompress)
    * [`unrar()`](acore.io.md#acore.io.uncompress.unrar)
* [acore.kaplan_meier_analysis package](acore.kaplan_meier_analysis.md)
  * [`get_data_ready_for_km()`](acore.kaplan_meier_analysis.md#acore.kaplan_meier_analysis.get_data_ready_for_km)
  * [`group_data_based_on_marker()`](acore.kaplan_meier_analysis.md#acore.kaplan_meier_analysis.group_data_based_on_marker)
  * [`run_km()`](acore.kaplan_meier_analysis.md#acore.kaplan_meier_analysis.run_km)
  * [`get_km_results()`](acore.kaplan_meier_analysis.md#acore.kaplan_meier_analysis.get_km_results)
  * [`get_hazard_ratio_results()`](acore.kaplan_meier_analysis.md#acore.kaplan_meier_analysis.get_hazard_ratio_results)
* [acore.multiple_testing package](acore.multiple_testing.md)
  * [`apply_pvalue_correction()`](acore.multiple_testing.md#acore.multiple_testing.apply_pvalue_correction)
  * [`apply_pvalue_fdrcorrection()`](acore.multiple_testing.md#acore.multiple_testing.apply_pvalue_fdrcorrection)
  * [`apply_pvalue_twostage_fdrcorrection()`](acore.multiple_testing.md#acore.multiple_testing.apply_pvalue_twostage_fdrcorrection)
  * [`apply_pvalue_permutation_fdrcorrection()`](acore.multiple_testing.md#acore.multiple_testing.apply_pvalue_permutation_fdrcorrection)
  * [`calculate_anova()`](acore.multiple_testing.md#acore.multiple_testing.calculate_anova)
  * [`get_counts_permutation_fdr()`](acore.multiple_testing.md#acore.multiple_testing.get_counts_permutation_fdr)
  * [`get_max_permutations()`](acore.multiple_testing.md#acore.multiple_testing.get_max_permutations)
  * [`correct_pairwise_ttest()`](acore.multiple_testing.md#acore.multiple_testing.correct_pairwise_ttest)
* [acore.network_analysis package](acore.network_analysis.md)
  * [`get_network_communities()`](acore.network_analysis.md#acore.network_analysis.get_network_communities)
  * [`get_snf_clusters()`](acore.network_analysis.md#acore.network_analysis.get_snf_clusters)
  * [`most_central_edge()`](acore.network_analysis.md#acore.network_analysis.most_central_edge)
  * [`get_louvain_partitions()`](acore.network_analysis.md#acore.network_analysis.get_louvain_partitions)
  * [`run_snf()`](acore.network_analysis.md#acore.network_analysis.run_snf)
* [acore.normalization package](acore.normalization.md)
  * [`normalize_data()`](acore.normalization.md#acore.normalization.normalize_data)
  * [`normalize_data_per_group()`](acore.normalization.md#acore.normalization.normalize_data_per_group)
  * [Submodules](acore.normalization.md#submodules)
  * [acore.normalization.strategies module](acore.normalization.md#module-acore.normalization.strategies)
    * [`median_zero_normalization()`](acore.normalization.md#acore.normalization.strategies.median_zero_normalization)
    * [`median_normalization()`](acore.normalization.md#acore.normalization.strategies.median_normalization)
    * [`zscore_normalization()`](acore.normalization.md#acore.normalization.strategies.zscore_normalization)
    * [`median_polish_normalization()`](acore.normalization.md#acore.normalization.strategies.median_polish_normalization)
    * [`quantile_normalization()`](acore.normalization.md#acore.normalization.strategies.quantile_normalization)
    * [`linear_normalization()`](acore.normalization.md#acore.normalization.strategies.linear_normalization)
* [acore.power_analysis package](acore.power_analysis.md)
  * [`power_analysis()`](acore.power_analysis.md#acore.power_analysis.power_analysis)
* [acore.publications_analysis package](acore.publications_analysis.md)
  * [`getMedlineAbstracts()`](acore.publications_analysis.md#acore.publications_analysis.getMedlineAbstracts)
  * [`get_publications_abstracts()`](acore.publications_analysis.md#acore.publications_analysis.get_publications_abstracts)
* [acore.sklearn package](acore.sklearn.md)
  * [`transform_DataFrame()`](acore.sklearn.md#acore.sklearn.transform_DataFrame)
* [acore.tda_analysis package](acore.tda_analysis.md)
  * [`run_mapper()`](acore.tda_analysis.md#acore.tda_analysis.run_mapper)
* [acore.types package](acore.types.md)
  * [`check_numeric_dataframe()`](acore.types.md#acore.types.check_numeric_dataframe)
  * [`select_numeric_columns()`](acore.types.md#acore.types.select_numeric_columns)
  * [`build_schema_all_floats()`](acore.types.md#acore.types.build_schema_all_floats)
  * [Submodules](acore.types.md#submodules)
  * [acore.types.differential_analysis module](acore.types.md#module-acore.types.differential_analysis)
    * [`AnovaSchema`](acore.types.md#acore.types.differential_analysis.AnovaSchema)
      * [`AnovaSchema.group1`](acore.types.md#acore.types.differential_analysis.AnovaSchema.group1)
      * [`AnovaSchema.group2`](acore.types.md#acore.types.differential_analysis.AnovaSchema.group2)
      * [`AnovaSchema.mean_group1`](acore.types.md#acore.types.differential_analysis.AnovaSchema.mean_group1)
      * [`AnovaSchema.std_group1`](acore.types.md#acore.types.differential_analysis.AnovaSchema.std_group1)
      * [`AnovaSchema.mean_group2`](acore.types.md#acore.types.differential_analysis.AnovaSchema.mean_group2)
      * [`AnovaSchema.std_group2`](acore.types.md#acore.types.differential_analysis.AnovaSchema.std_group2)
      * [`AnovaSchema.t_statistics`](acore.types.md#acore.types.differential_analysis.AnovaSchema.t_statistics)
      * [`AnovaSchema.pvalue`](acore.types.md#acore.types.differential_analysis.AnovaSchema.pvalue)
      * [`AnovaSchema.log2FC`](acore.types.md#acore.types.differential_analysis.AnovaSchema.log2FC)
      * [`AnovaSchema.FC`](acore.types.md#acore.types.differential_analysis.AnovaSchema.FC)
      * [`AnovaSchema.padj`](acore.types.md#acore.types.differential_analysis.AnovaSchema.padj)
      * [`AnovaSchema.correction`](acore.types.md#acore.types.differential_analysis.AnovaSchema.correction)
      * [`AnovaSchema.rejected`](acore.types.md#acore.types.differential_analysis.AnovaSchema.rejected)
      * [`AnovaSchema.neg_log10_p_value`](acore.types.md#acore.types.differential_analysis.AnovaSchema.neg_log10_p_value)
      * [`AnovaSchema.Method`](acore.types.md#acore.types.differential_analysis.AnovaSchema.Method)
      * [`AnovaSchema.Config`](acore.types.md#acore.types.differential_analysis.AnovaSchema.Config)
      * [`AnovaSchema.build_schema_()`](acore.types.md#acore.types.differential_analysis.AnovaSchema.build_schema_)
      * [`AnovaSchema.empty()`](acore.types.md#acore.types.differential_analysis.AnovaSchema.empty)
      * [`AnovaSchema.example()`](acore.types.md#acore.types.differential_analysis.AnovaSchema.example)
      * [`AnovaSchema.from_json()`](acore.types.md#acore.types.differential_analysis.AnovaSchema.from_json)
      * [`AnovaSchema.from_yaml()`](acore.types.md#acore.types.differential_analysis.AnovaSchema.from_yaml)
      * [`AnovaSchema.get_metadata()`](acore.types.md#acore.types.differential_analysis.AnovaSchema.get_metadata)
      * [`AnovaSchema.pydantic_validate()`](acore.types.md#acore.types.differential_analysis.AnovaSchema.pydantic_validate)
      * [`AnovaSchema.strategy()`](acore.types.md#acore.types.differential_analysis.AnovaSchema.strategy)
      * [`AnovaSchema.to_json()`](acore.types.md#acore.types.differential_analysis.AnovaSchema.to_json)
      * [`AnovaSchema.to_json_schema()`](acore.types.md#acore.types.differential_analysis.AnovaSchema.to_json_schema)
      * [`AnovaSchema.to_schema()`](acore.types.md#acore.types.differential_analysis.AnovaSchema.to_schema)
      * [`AnovaSchema.to_yaml()`](acore.types.md#acore.types.differential_analysis.AnovaSchema.to_yaml)
      * [`AnovaSchema.validate()`](acore.types.md#acore.types.differential_analysis.AnovaSchema.validate)
    * [`AnovaSchemaMultiGroup`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup)
      * [`AnovaSchemaMultiGroup.t_statistics`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.t_statistics)
      * [`AnovaSchemaMultiGroup.posthoc_pvalue`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.posthoc_pvalue)
      * [`AnovaSchemaMultiGroup.f_statistics`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.f_statistics)
      * [`AnovaSchemaMultiGroup.posthoc_padj`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.posthoc_padj)
      * [`AnovaSchemaMultiGroup.posthoc_paired`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.posthoc_paired)
      * [`AnovaSchemaMultiGroup.posthoc_parametric`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.posthoc_parametric)
      * [`AnovaSchemaMultiGroup.posthoc_dof`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.posthoc_dof)
      * [`AnovaSchemaMultiGroup.posthoc_tail`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.posthoc_tail)
      * [`AnovaSchemaMultiGroup.posthoc_BF10`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.posthoc_BF10)
      * [`AnovaSchemaMultiGroup.posthoc_effsize`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.posthoc_effsize)
      * [`AnovaSchemaMultiGroup.efftype`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.efftype)
      * [`AnovaSchemaMultiGroup.Config`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.Config)
      * [`AnovaSchemaMultiGroup.FC`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.FC)
      * [`AnovaSchemaMultiGroup.Method`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.Method)
      * [`AnovaSchemaMultiGroup.build_schema_()`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.build_schema_)
      * [`AnovaSchemaMultiGroup.correction`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.correction)
      * [`AnovaSchemaMultiGroup.empty()`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.empty)
      * [`AnovaSchemaMultiGroup.example()`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.example)
      * [`AnovaSchemaMultiGroup.from_json()`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.from_json)
      * [`AnovaSchemaMultiGroup.from_yaml()`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.from_yaml)
      * [`AnovaSchemaMultiGroup.get_metadata()`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.get_metadata)
      * [`AnovaSchemaMultiGroup.group1`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.group1)
      * [`AnovaSchemaMultiGroup.group2`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.group2)
      * [`AnovaSchemaMultiGroup.log2FC`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.log2FC)
      * [`AnovaSchemaMultiGroup.mean_group1`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.mean_group1)
      * [`AnovaSchemaMultiGroup.mean_group2`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.mean_group2)
      * [`AnovaSchemaMultiGroup.neg_log10_p_value`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.neg_log10_p_value)
      * [`AnovaSchemaMultiGroup.padj`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.padj)
      * [`AnovaSchemaMultiGroup.pvalue`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.pvalue)
      * [`AnovaSchemaMultiGroup.pydantic_validate()`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.pydantic_validate)
      * [`AnovaSchemaMultiGroup.rejected`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.rejected)
      * [`AnovaSchemaMultiGroup.std_group1`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.std_group1)
      * [`AnovaSchemaMultiGroup.std_group2`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.std_group2)
      * [`AnovaSchemaMultiGroup.strategy()`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.strategy)
      * [`AnovaSchemaMultiGroup.to_json()`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.to_json)
      * [`AnovaSchemaMultiGroup.to_json_schema()`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.to_json_schema)
      * [`AnovaSchemaMultiGroup.to_schema()`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.to_schema)
      * [`AnovaSchemaMultiGroup.to_yaml()`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.to_yaml)
      * [`AnovaSchemaMultiGroup.validate()`](acore.types.md#acore.types.differential_analysis.AnovaSchemaMultiGroup.validate)
    * [`AncovaSchema`](acore.types.md#acore.types.differential_analysis.AncovaSchema)
      * [`AncovaSchema.t_statistics`](acore.types.md#acore.types.differential_analysis.AncovaSchema.t_statistics)
      * [`AncovaSchema.posthoc_pvalue`](acore.types.md#acore.types.differential_analysis.AncovaSchema.posthoc_pvalue)
      * [`AncovaSchema.coef`](acore.types.md#acore.types.differential_analysis.AncovaSchema.coef)
      * [`AncovaSchema.std_err`](acore.types.md#acore.types.differential_analysis.AncovaSchema.std_err)
      * [`AncovaSchema.conf_int_low`](acore.types.md#acore.types.differential_analysis.AncovaSchema.conf_int_low)
      * [`AncovaSchema.conf_int_upp`](acore.types.md#acore.types.differential_analysis.AncovaSchema.conf_int_upp)
      * [`AncovaSchema.f_statistics`](acore.types.md#acore.types.differential_analysis.AncovaSchema.f_statistics)
      * [`AncovaSchema.posthoc_padj`](acore.types.md#acore.types.differential_analysis.AncovaSchema.posthoc_padj)
      * [`AncovaSchema.Config`](acore.types.md#acore.types.differential_analysis.AncovaSchema.Config)
      * [`AncovaSchema.build_schema_()`](acore.types.md#acore.types.differential_analysis.AncovaSchema.build_schema_)
      * [`AncovaSchema.empty()`](acore.types.md#acore.types.differential_analysis.AncovaSchema.empty)
      * [`AncovaSchema.example()`](acore.types.md#acore.types.differential_analysis.AncovaSchema.example)
      * [`AncovaSchema.from_json()`](acore.types.md#acore.types.differential_analysis.AncovaSchema.from_json)
      * [`AncovaSchema.from_yaml()`](acore.types.md#acore.types.differential_analysis.AncovaSchema.from_yaml)
      * [`AncovaSchema.get_metadata()`](acore.types.md#acore.types.differential_analysis.AncovaSchema.get_metadata)
      * [`AncovaSchema.pydantic_validate()`](acore.types.md#acore.types.differential_analysis.AncovaSchema.pydantic_validate)
      * [`AncovaSchema.strategy()`](acore.types.md#acore.types.differential_analysis.AncovaSchema.strategy)
      * [`AncovaSchema.to_json()`](acore.types.md#acore.types.differential_analysis.AncovaSchema.to_json)
      * [`AncovaSchema.to_json_schema()`](acore.types.md#acore.types.differential_analysis.AncovaSchema.to_json_schema)
      * [`AncovaSchema.to_schema()`](acore.types.md#acore.types.differential_analysis.AncovaSchema.to_schema)
      * [`AncovaSchema.to_yaml()`](acore.types.md#acore.types.differential_analysis.AncovaSchema.to_yaml)
      * [`AncovaSchema.validate()`](acore.types.md#acore.types.differential_analysis.AncovaSchema.validate)
  * [acore.types.enrichment_analysis module](acore.types.md#module-acore.types.enrichment_analysis)
    * [`EnrichmentAnalysisSchema`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema)
      * [`EnrichmentAnalysisSchema.terms`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.terms)
      * [`EnrichmentAnalysisSchema.identifiers`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.identifiers)
      * [`EnrichmentAnalysisSchema.foreground`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.foreground)
      * [`EnrichmentAnalysisSchema.background`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.background)
      * [`EnrichmentAnalysisSchema.foreground_pop`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.foreground_pop)
      * [`EnrichmentAnalysisSchema.background_pop`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.background_pop)
      * [`EnrichmentAnalysisSchema.pvalue`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.pvalue)
      * [`EnrichmentAnalysisSchema.padj`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.padj)
      * [`EnrichmentAnalysisSchema.rejected`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.rejected)
      * [`EnrichmentAnalysisSchema.direction`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.direction)
      * [`EnrichmentAnalysisSchema.comparison`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.comparison)
      * [`EnrichmentAnalysisSchema.Config`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.Config)
      * [`EnrichmentAnalysisSchema.build_schema_()`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.build_schema_)
      * [`EnrichmentAnalysisSchema.empty()`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.empty)
      * [`EnrichmentAnalysisSchema.example()`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.example)
      * [`EnrichmentAnalysisSchema.from_json()`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.from_json)
      * [`EnrichmentAnalysisSchema.from_yaml()`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.from_yaml)
      * [`EnrichmentAnalysisSchema.get_metadata()`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.get_metadata)
      * [`EnrichmentAnalysisSchema.pydantic_validate()`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.pydantic_validate)
      * [`EnrichmentAnalysisSchema.strategy()`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.strategy)
      * [`EnrichmentAnalysisSchema.to_json()`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.to_json)
      * [`EnrichmentAnalysisSchema.to_json_schema()`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.to_json_schema)
      * [`EnrichmentAnalysisSchema.to_schema()`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.to_schema)
      * [`EnrichmentAnalysisSchema.to_yaml()`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.to_yaml)
      * [`EnrichmentAnalysisSchema.validate()`](acore.types.md#acore.types.enrichment_analysis.EnrichmentAnalysisSchema.validate)
  * [acore.types.exploratory_analysis module](acore.types.md#module-acore.types.exploratory_analysis)
    * [`TwoComponentSchema`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema)
      * [`TwoComponentSchema.group`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.group)
      * [`TwoComponentSchema.x`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.x)
      * [`TwoComponentSchema.y`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.y)
      * [`TwoComponentSchema.Config`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.Config)
      * [`TwoComponentSchema.build_schema_()`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.build_schema_)
      * [`TwoComponentSchema.empty()`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.empty)
      * [`TwoComponentSchema.example()`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.example)
      * [`TwoComponentSchema.from_json()`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.from_json)
      * [`TwoComponentSchema.from_yaml()`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.from_yaml)
      * [`TwoComponentSchema.get_metadata()`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.get_metadata)
      * [`TwoComponentSchema.pydantic_validate()`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.pydantic_validate)
      * [`TwoComponentSchema.strategy()`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.strategy)
      * [`TwoComponentSchema.to_json()`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.to_json)
      * [`TwoComponentSchema.to_json_schema()`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.to_json_schema)
      * [`TwoComponentSchema.to_schema()`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.to_schema)
      * [`TwoComponentSchema.to_yaml()`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.to_yaml)
      * [`TwoComponentSchema.validate()`](acore.types.md#acore.types.exploratory_analysis.TwoComponentSchema.validate)
    * [`TwoLoadingsSchema`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema)
      * [`TwoLoadingsSchema.value`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.value)
      * [`TwoLoadingsSchema.Config`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.Config)
      * [`TwoLoadingsSchema.build_schema_()`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.build_schema_)
      * [`TwoLoadingsSchema.empty()`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.empty)
      * [`TwoLoadingsSchema.example()`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.example)
      * [`TwoLoadingsSchema.from_json()`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.from_json)
      * [`TwoLoadingsSchema.from_yaml()`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.from_yaml)
      * [`TwoLoadingsSchema.get_metadata()`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.get_metadata)
      * [`TwoLoadingsSchema.group`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.group)
      * [`TwoLoadingsSchema.pydantic_validate()`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.pydantic_validate)
      * [`TwoLoadingsSchema.strategy()`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.strategy)
      * [`TwoLoadingsSchema.to_json()`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.to_json)
      * [`TwoLoadingsSchema.to_json_schema()`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.to_json_schema)
      * [`TwoLoadingsSchema.to_schema()`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.to_schema)
      * [`TwoLoadingsSchema.to_yaml()`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.to_yaml)
      * [`TwoLoadingsSchema.validate()`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.validate)
      * [`TwoLoadingsSchema.x`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.x)
      * [`TwoLoadingsSchema.y`](acore.types.md#acore.types.exploratory_analysis.TwoLoadingsSchema.y)
    * [`AnnotationResult`](acore.types.md#acore.types.exploratory_analysis.AnnotationResult)
      * [`AnnotationResult.x_title`](acore.types.md#acore.types.exploratory_analysis.AnnotationResult.x_title)
      * [`AnnotationResult.y_title`](acore.types.md#acore.types.exploratory_analysis.AnnotationResult.y_title)
      * [`AnnotationResult.group`](acore.types.md#acore.types.exploratory_analysis.AnnotationResult.group)
      * [`AnnotationResult.model_config`](acore.types.md#acore.types.exploratory_analysis.AnnotationResult.model_config)

## Submodules

## acore.utils module

### convertToEdgeList(data, cols)

This function converts a pandas dataframe to an edge list where index becomes the source nodes and columns the target nodes.

* **Parameters:**
  * **data** – pandas dataframe.
  * **cols** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – names for dataframe columns.
* **Returns:**
  Pandas dataframe with columns cols.

### check_is_paired(df, subject, group)

Check if samples are paired.

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and protein identifiers as columns (with additional columns ‘group’, ‘sample’ and ‘subject’).
  * **subject** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with subject identifiers
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with group identifiers
* **Returns:**
  True if paired samples.
* **Return type:**
  [bool](https://docs.python.org/3/library/functions.html#bool)

### transform_into_wide_format(data, index, columns, values, extra=[])

This function converts a Pandas DataFrame from long to wide format using
pandas pivot_table() function.

* **Parameters:**
  * **data** – long-format Pandas DataFrame
  * **index** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – columns that will be converted into the index
  * **columns** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column name whose unique values will become the new column names
  * **values** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column to aggregate
  * **extra** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – additional columns to be kept as columns
* **Returns:**
  Wide-format pandas DataFrame

Example:

```default
result = transform_into_wide_format(df, index='index', columns='x', values='y', extra='group')
```

### transform_into_long_format(data, drop_columns, group, columns=['name', 'y'])

Converts a Pandas DataDrame from wide to long format using pd.melt()
function.

* **Parameters:**
  * **data** – wide-format Pandas DataFrame
  * **drop_columns** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – columns to be deleted
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str) *or* [*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column(s) to use as identifier variables
  * **columns** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – names to use for the 1)variable column, and for the 2)value column
* **Returns:**
  Long-format Pandas DataFrame.

Example:

```default
result = transform_into_long_format(df, drop_columns=['sample', 'subject'], group='group', columns=['name','y'])
```

### remove_group(data)

Removes column with label ‘group’.

* **Parameters:**
  **data** – pandas dataframe with one column labelled ‘group’
* **Returns:**
  Pandas dataframe

Example:

```default
result = remove_group(data)
```

### calculate_fold_change(df, condition1, condition2)

Calculates fold-changes between two groups for all proteins in a dataframe.

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and protein identifiers as columns.
  * **condition1** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – identifier of first group.
  * **condition2** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – identifier of second group.
* **Returns:**
  Numpy array.

Example:

```default
result = calculate_fold_change(data, 'group1', 'group2')
```

### pooled_standard_deviation(sample1, sample2, ddof)

Calculates the pooled standard deviation.
For more information visit [https://www.hackdeploy.com/learn-what-is-statistical-power-with-python/](https://www.hackdeploy.com/learn-what-is-statistical-power-with-python/).

* **Parameters:**
  * **sample1** (*array*) – numpy array with values for first group
  * **sample2** (*array*) – numpy array with values for second group
  * **ddof** ([*int*](https://docs.python.org/3/library/functions.html#int)) – degrees of freedom

### cohens_d(sample1, sample2, ddof)

Calculates Cohen’s d effect size based on the distance between two means, measured in standard deviations.
For more information visit [https://www.hackdeploy.com/learn-what-is-statistical-power-with-python/](https://www.hackdeploy.com/learn-what-is-statistical-power-with-python/).

* **Parameters:**
  * **sample1** (*array*) – numpy array with values for first group
  * **sample2** (*array*) – numpy array with values for second group
  * **ddof** ([*int*](https://docs.python.org/3/library/functions.html#int)) – degrees of freedom

### hedges_g(df, condition1, condition2, ddof=0)

Calculates Hedges’ g effect size (more accurate for sample sizes below 20 than Cohen’s d).
For more information visit [https://docs.scipy.org/doc/numpy/reference/generated/numpy.nanstd.html](https://docs.scipy.org/doc/numpy/reference/generated/numpy.nanstd.html).

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and protein identifiers as columns.
  * **condition1** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – identifier of first group.
  * **condition2** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – identifier of second group.
  * **ddof** ([*int*](https://docs.python.org/3/library/functions.html#int)) – means Delta Degrees of Freedom.
* **Returns:**
  Numpy array.

Example:

```default
result = hedges_g(data, 'group1', 'group2', ddof=0)
```

### unit_vector(vector)

Returns the unit vector of the vector.
:param tuple vector: vector
:return tuple unit_vector: unit vector

### flatten(t, my_list=[])

Code from: [https://gist.github.com/shaxbee/0ada767debf9eefbdb6e](https://gist.github.com/shaxbee/0ada767debf9eefbdb6e)
Acknowledgements: Zbigniew Mandziejewicz (shaxbee)
Generator flattening the structure

```pycon
>>> list(flatten([2, [2, (4, 5, [7], [2, [6, 2, 6, [6], 4]], 6)]]))
[2, 2, 4, 5, 7, 2, 6, 2, 6, 6, 4, 6]
```

### angle_between(v1, v2)

Returns the angle in radians between vectors ‘v1’ and ‘v2’

* **Parameters:**
  * **v1** ([*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – vector 1
  * **v2** ([*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – vector 2
* **Return float angle:**
  angle between two vectors in radians

Example::
: angle = angle_between((1, 0, 0), (0, 1, 0))

### append_to_list(mylist, myappend)

### generator_to_dict(genvar)

## acore.wgcna_analysis module
