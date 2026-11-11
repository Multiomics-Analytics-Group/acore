"""A column name that is not in the data is reported with the columns that are.

Callers of acore, such as scripts and LLM tool-calling agents, pass column
names as arguments. When one is wrong, pandas raises a bare ``KeyError`` such as
``['subject'] not found in axis``, which names neither the argument at fault nor
the columns they could have used instead.

A column that is only listed to be dropped and is already missing is no error:
the data is as the caller wants it, so it is a warning.
"""

import warnings

import pandas as pd
import pytest

from acore import (
    batch_correction,
    correlation_analysis,
    differential_regulation,
    exploratory_analysis,
    imputation_analysis,
)
from acore.utils import check_columns_in_data, columns_to_drop


@pytest.fixture
def measurements() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample": [f"S{number}" for number in range(1, 13)],
            "subject": [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4],
            "group": ["A", "B", "C"] * 4,
            "visit": ["v1", "v2", "v3"] * 4,
            "batch": ["b1", "b1", "b2", "b2"] * 3,
            "protein_a": [
                1.4,
                6.2,
                9.0,
                2.3,
                7.4,
                14.0,
                4.5,
                9.6,
                10.4,
                7.5,
                11.6,
                16.4,
            ],
            "protein_b": [3.1, 2.2, 8.1, 4.0, 5.5, 9.9, 6.1, 7.2, 3.3, 8.8, 1.9, 2.7],
            "protein_c": [5.0, 4.1, 3.2, 6.3, 1.5, 2.6, 7.7, 8.1, 9.2, 0.4, 5.9, 6.8],
        }
    )


@pytest.fixture
def two_group_measurements(measurements: pd.DataFrame) -> pd.DataFrame:
    return measurements[measurements["group"].isin(["A", "B"])].copy()


# check_columns_in_data
class TestCheckColumnsInData:
    def test_names_the_parameter_the_missing_column_and_the_available_columns(
        self, measurements
    ):
        # Given a table without a 'treatment' column

        # When that column is requested through the 'group' parameter
        with pytest.raises(ValueError) as error:
            check_columns_in_data(measurements, group="treatment")

        # Then the message names the parameter, the column and the columns to use
        message = str(error.value)
        assert "group='treatment'" in message
        assert "Available columns: ['sample', 'subject', 'group'" in message

    def test_names_only_the_missing_columns_of_a_list(self, measurements):
        # Given a list in which one column exists and one does not
        columns_to_check = ["sample", "sample-id"]

        # When the list is checked
        with pytest.raises(ValueError) as error:
            check_columns_in_data(measurements, subject_ids=columns_to_check)

        # Then only the missing column is named
        message = str(error.value)
        assert "subject_ids=['sample-id']" in message
        assert "'sample', 'sample-id'" not in message

    def test_accepts_columns_that_exist_and_parameters_set_to_none(self, measurements):
        # Given parameters that exist or are not set

        # When they are checked
        check_columns_in_data(measurements, group="group", subject=None)

        # Then nothing is raised

    def test_accepts_a_tuple_that_is_one_column_label(self):
        # Given a table with MultiIndex columns
        multi_index_table = pd.DataFrame(
            [[1, "A"], [2, "B"]],
            columns=pd.MultiIndex.from_tuples([("meta", "id"), ("meta", "group")]),
        )

        # When the tuple naming one of its columns is checked
        check_columns_in_data(multi_index_table, group=("meta", "group"))

        # Then nothing is raised

    def test_shortens_a_long_list_of_available_columns(self):
        # Given a table with 50 columns
        wide_table = pd.DataFrame({f"protein_{number}": [1.0] for number in range(50)})

        # When a column that does not exist is requested
        with pytest.raises(ValueError) as error:
            check_columns_in_data(wide_table, group="group")

        # Then only the first columns are listed, and the rest are counted
        message = str(error.value)
        assert "'protein_0'" in message
        assert "'protein_49'" not in message
        assert "30 more" in message


# columns_to_drop
class TestColumnsToDrop:
    def test_returns_the_columns_that_exist_and_warns_about_the_others(
        self, measurements
    ):
        # Given a list in which one column exists and one does not
        columns = ["sample", "sample-id"]

        # When the columns to drop are resolved
        with pytest.warns(UserWarning) as recorded:
            existing_columns = columns_to_drop(measurements, columns, "drop_cols")

        # Then only the existing column is returned, and the other is named
        assert existing_columns == ["sample"]
        assert "drop_cols=['sample-id']" in str(recorded[0].message)

    def test_reads_an_iterator_only_once(self, measurements):
        # Given an iterator, which can be read only once
        columns = (column for column in ["sample", "subject"])

        # When the columns to drop are resolved
        existing_columns = columns_to_drop(measurements, columns, "drop_cols")

        # Then the result is a list that can be read again
        assert existing_columns == ["sample", "subject"]
        assert existing_columns == ["sample", "subject"]

    def test_does_not_warn_when_every_column_exists(self, measurements):
        # Given columns that all exist

        # When the columns to drop are resolved with warnings turned into errors
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            existing_columns = columns_to_drop(measurements, ["sample"], "drop_cols")

        # Then nothing is raised
        assert existing_columns == ["sample"]

    def test_returns_no_columns_for_none(self, measurements):
        # Given no columns to drop

        # When they are resolved
        existing_columns = columns_to_drop(measurements, None, "drop_cols")

        # Then the list is empty
        assert existing_columns == []


# differential_regulation
class TestDifferentialRegulationColumnErrors:
    def test_run_ttest_names_a_missing_subject(self, two_group_measurements):
        # Given a table without a 'subject_id' column

        # When a paired t-test is asked to pair on it
        with pytest.raises(ValueError) as error:
            differential_regulation.run_ttest(
                two_group_measurements,
                condition1="A",
                condition2="B",
                subject="subject_id",
                group="group",
                drop_cols=["sample"],
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "subject='subject_id'" in message
        assert "Available columns:" in message

    def test_run_ttest_names_a_missing_group(self, two_group_measurements):
        # Given a table without a 'condition' column

        # When it is used to label the groups
        with pytest.raises(ValueError) as error:
            differential_regulation.run_ttest(
                two_group_measurements,
                condition1="A",
                condition2="B",
                subject=None,
                group="condition",
                drop_cols=["sample"],
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "group='condition'" in message
        assert "Available columns:" in message

    def test_run_ttest_warns_about_a_missing_column_to_drop(
        self, two_group_measurements
    ):
        # Given a table without a 'sample_id' column

        # When it is listed in drop_cols
        with pytest.warns(UserWarning, match=r"drop_cols=\['sample_id'\]"):
            result = differential_regulation.run_ttest(
                two_group_measurements,
                condition1="A",
                condition2="B",
                subject=None,
                group="group",
                drop_cols=["sample", "subject", "visit", "batch", "sample_id"],
            )

        # Then the test still runs on the three proteins
        assert len(result) == 3

    def test_run_anova_names_a_missing_group(self, measurements):
        # Given a table without a 'condition' column

        # When it is used to label the groups
        with pytest.raises(ValueError) as error:
            differential_regulation.run_anova(
                measurements,
                subject=None,
                group="condition",
                drop_cols=["sample", "subject", "visit", "batch"],
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "group='condition'" in message
        assert "Available columns:" in message

    def test_run_ancova_names_a_missing_covariate(self, measurements):
        # Given a table without an 'age' column

        # When it is used as a covariate
        with pytest.raises(ValueError) as error:
            differential_regulation.run_ancova(
                measurements,
                covariates=["age"],
                drop_cols=["sample", "subject", "visit", "batch"],
                group="group",
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "covariates=['age']" in message
        assert "Available columns:" in message

    def test_run_repeated_measurements_anova_names_a_missing_subject(
        self, measurements
    ):
        # Given a table without a 'subject_id' column

        # When repeated measurements are linked through it
        with pytest.raises(ValueError) as error:
            differential_regulation.run_repeated_measurements_anova(
                measurements,
                subject="subject_id",
                within="visit",
                drop_cols=["sample", "group", "batch"],
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "subject='subject_id'" in message
        assert "Available columns:" in message

    def test_run_mixed_anova_names_a_missing_between_column(self, measurements):
        # Given a table without a 'treatment' column

        # When it is used as the between-subjects factor
        with pytest.raises(ValueError) as error:
            differential_regulation.run_mixed_anova(
                measurements,
                subject="subject",
                within="visit",
                between="treatment",
                drop_cols=["sample", "group", "batch"],
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "between='treatment'" in message
        assert "Available columns:" in message


# exploratory_analysis
class TestExploratoryAnalysisColumnErrors:
    def test_run_pca_names_a_missing_group(self, measurements):
        # Given a table without a 'collection_site' column

        # When samples are coloured by it
        with pytest.raises(ValueError) as error:
            exploratory_analysis.run_pca(
                measurements,
                group="collection_site",
                drop_cols=["sample", "subject", "visit", "batch"],
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "group='collection_site'" in message
        assert "Available columns:" in message

    def test_run_pca_warns_about_a_missing_column_to_drop(self, measurements):
        # Given a table without a 'sample_id' column

        # When it is listed in drop_cols
        with pytest.warns(UserWarning, match=r"drop_cols=\['sample_id'\]"):
            (components, _, _), _ = exploratory_analysis.run_pca(
                measurements,
                group="group",
                drop_cols=["sample", "subject", "visit", "batch", "sample_id"],
            )

        # Then the analysis still runs and returns a row per sample
        assert len(components) == 12

    def test_run_umap_names_a_missing_group(self, measurements):
        # Given a table without a 'diagnosis' column

        # When samples are coloured by it
        with pytest.raises(ValueError) as error:
            exploratory_analysis.run_umap(
                measurements,
                group="diagnosis",
                drop_cols=["sample", "subject", "visit", "batch"],
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "group='diagnosis'" in message
        assert "Available columns:" in message

    def test_run_tsne_names_a_missing_group(self, measurements):
        # Given a table without a 'diagnosis' column

        # When samples are coloured by it
        with pytest.raises(ValueError) as error:
            exploratory_analysis.run_tsne(
                measurements,
                group="diagnosis",
                drop_cols=["sample", "subject", "visit", "batch"],
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "group='diagnosis'" in message
        assert "Available columns:" in message

    def test_extract_number_missing_names_a_missing_group(self, measurements):
        # Given a table without a 'diagnosis' column

        # When valid values are counted per diagnosis
        with pytest.raises(ValueError) as error:
            exploratory_analysis.extract_number_missing(
                measurements, min_valid=3, group="diagnosis"
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "group='diagnosis'" in message
        assert "Available columns:" in message

    def test_extract_number_missing_drops_the_columns_that_exist(self, measurements):
        # Given a table without a 'sample_id' column

        # When it is listed in drop_cols beside columns that exist
        with pytest.warns(UserWarning, match=r"drop_cols=\['sample_id'\]"):
            kept_columns = exploratory_analysis.extract_number_missing(
                measurements,
                min_valid=1,
                group="group",
                drop_cols=["sample", "subject", "visit", "batch", "sample_id"],
            )

        # Then the columns that exist are dropped and the proteins are kept
        assert kept_columns == ["protein_a", "protein_b", "protein_c"]

    def test_extract_percentage_missing_names_a_missing_group(self, measurements):
        # Given a table without a 'diagnosis' column

        # When missing values are measured per diagnosis
        with pytest.raises(ValueError) as error:
            exploratory_analysis.extract_percentage_missing(
                measurements, missing_max=0.3, group="diagnosis"
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "group='diagnosis'" in message
        assert "Available columns:" in message

    def test_extract_percentage_missing_warns_when_its_default_column_is_missing(
        self, measurements
    ):
        # Given a table without a 'sample' column, which drop_cols defaults to
        table_without_sample = measurements.drop(
            columns=["sample", "subject", "visit", "batch"]
        )

        # When drop_cols is not given
        with pytest.warns(UserWarning, match=r"drop_cols=\['sample'\]"):
            kept_columns = exploratory_analysis.extract_percentage_missing(
                table_without_sample, missing_max=0.3, group="group"
            )

        # Then the proteins are kept
        assert kept_columns == ["protein_a", "protein_b", "protein_c"]

    def test_get_coefficient_variation_names_a_missing_group(self, measurements):
        # Given a table without a 'diagnosis' column

        # When the coefficient of variation is taken per diagnosis
        with pytest.raises(ValueError) as error:
            exploratory_analysis.get_coefficient_variation(
                measurements, drop_columns=["sample"], group="diagnosis"
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "group='diagnosis'" in message
        assert "Available columns:" in message

    def test_get_coefficient_variation_warns_about_a_missing_column_to_drop(
        self, measurements
    ):
        # Given a table without a 'sample_id' column
        numeric_measurements = measurements.drop(
            columns=["sample", "subject", "visit", "batch"]
        )

        # When it is listed in drop_columns
        with pytest.warns(UserWarning, match=r"drop_columns=\['sample_id'\]"):
            result = exploratory_analysis.get_coefficient_variation(
                numeric_measurements, drop_columns=["sample_id"], group="group"
            )

        # Then the coefficients of variation are still returned
        assert not result.empty


# batch_correction, correlation_analysis and imputation_analysis
class TestOtherModulesColumnErrors:
    def test_combat_batch_correction_names_a_missing_batch_column(self, measurements):
        # Given a table without a 'collection_site' column

        # When it is used to identify the batches
        with pytest.raises(ValueError) as error:
            batch_correction.combat_batch_correction(
                measurements, batch_col="collection_site"
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "batch_col='collection_site'" in message
        assert "Available columns:" in message

    def test_run_correlation_names_a_missing_subject(self, measurements):
        # Given a table without a 'subject_id' column

        # When repeated measurements are linked through it
        with pytest.raises(ValueError) as error:
            correlation_analysis.run_correlation(
                measurements, subject="subject_id", group="group"
            )

        # Then the message names the argument and lists the available columns
        message = str(error.value)
        assert "subject='subject_id'" in message
        assert "Available columns:" in message

    def test_imputation_zeros_warns_about_a_missing_column_to_drop(self, measurements):
        # Given a table without a 'sample-id' column

        # When it is listed in drop_cols
        with pytest.warns(UserWarning, match=r"drop_cols=\['sample-id'\]"):
            result = imputation_analysis.imputation_zeros(
                measurements, drop_cols=["sample", "sample-id"], group="group"
            )

        # Then the columns that exist are dropped
        assert "sample" not in result.columns

    def test_imputation_zeros_accepts_an_iterator_of_columns_to_drop(
        self, measurements
    ):
        # Given the columns to drop as an iterator, which can be read only once
        columns_to_drop_once = (column for column in ["sample", "subject"])

        # When the table is imputed
        result = imputation_analysis.imputation_zeros(
            measurements, drop_cols=columns_to_drop_once, group="group"
        )

        # Then both columns are dropped
        assert "sample" not in result.columns
        assert "subject" not in result.columns
