import unittest

import pandas as pd
import pytest
from scipy import stats

import acore.enrichment_analysis as ea


class TestRunFisher(unittest.TestCase):
    def test_run_fisher(self):
        group1 = [10, 5]
        group2 = [8, 12]
        alternative = "two-sided"

        expected_odds, expected_pvalue = stats.fisher_exact(
            [[10, 5], [8, 12]], alternative
        )

        result = ea.run_fisher(group1, group2, alternative=alternative)

        self.assertEqual(result[0], expected_odds)
        self.assertEqual(result[1], expected_pvalue)


class TestRunKolmogorovSmirnov(unittest.TestCase):
    def test_run_kolmogorov_smirnov(self):
        dist1 = [1, 2, 3, 4, 5]
        dist2 = [1, 2, 3, 4, 6]
        alternative = "two-sided"

        expected_result = stats.ks_2samp(
            dist1, dist2, alternative=alternative, mode="auto"
        )

        result = ea.run_kolmogorov_smirnov(dist1, dist2, alternative=alternative)

        self.assertEqual(result[0], expected_result.statistic)
        self.assertEqual(result[1], expected_result.pvalue)


def test_run_regulation_enrichment():
    """Integration test for run_regulation_enrichment. Indirectly tests
    run_enrichment from enrichment_analysis module."""
    annotation = {
        "annotation": ["path1", "path1", "path1", "path2", "path2", "path3", "path3"],
        "identifier": ["gene1", "gene2", "gene3", "gene1", "gene5", "gene6", "gene9"],
        "source": ["GO", "GO", "GO", "GO_P", "GO_P", "GO_P", "GO_P"],
    }
    annotation = pd.DataFrame(annotation)
    regulation_res = {
        "identifier": ["gene1", "gene2", "gene3", "gene4", "gene5", "gene6"],
        "rejected": [True, True, False, False, True, True],
    }
    regulation_res = pd.DataFrame(regulation_res)

    actual = ea.run_regulation_enrichment(
        regulation_data=regulation_res,
        annotation=annotation,
        min_detected_in_set=1,
    )

    expected = pd.DataFrame(
        {
            "terms": ["path1", "path2", "path3"],
            "identifiers": ["gene1,gene2", "gene1,gene5", "gene6"],
            "foreground": [2, 2, 1],
            "background": [1, 0, 0],
            "foreground_pop": [4, 4, 4],
            "background_pop": [6, 6, 6],
            "pvalue": [1.0, 0.4666666666666667, 1.0],
            "padj": [1.0, 1.0, 1.0],
            "rejected": [False, False, False],
        }
    )
    assert expected.equals(actual)


def test_run_regulation_enrichment_pep():
    """Integration test for run_regulation_enrichment on peptides level data.
    Indirectly tests run_enrichment from enrichment_analysis module."""
    annotation = {
        # annotations from UNIPROT on protein level have to be exploded to peptides
        # e.g. by matching the gene identifiers to the differential analysis
        "annotation": [
            "path1",
            "path1",
            "path1",
            "path1",
            "path2",
            "path2",
            "path3",
            "path3",
        ],
        "identifier": [
            "gene1_pep1",
            "gene2_pep1",
            "gene2_pep2",
            "gene3_pep1",
            "gene1_pep3",
            "gene5_pep1",
            "gene6_pep1",
            "gene9_pep1",
        ],
        "gene": [
            "gene1",
            "gene2",
            "gene2",  # duplicated for pathway 1
            "gene3",
            "gene1",
            "gene5",
            "gene6",
            "gene9",
        ],
        "source": ["GO", "GO", "GO", "GO", "GO_P", "GO_P", "GO_P", "GO_P"],
    }
    annotation = pd.DataFrame(annotation)
    regulation_res = {
        # one peptide in forground, one in background for gene2
        "identifier": [
            "gene1_pep1",
            "gene2_pep1",
            "gene2_pep2",
            "gene3_pep1",
            "gene1_pep3",
            "gene5_pep1",
            "gene6_pep1",
            "gene9_pep1",
        ],
        "rejected": [True, True, False, True, False, False, True, True],
    }
    regulation_res = pd.DataFrame(regulation_res)

    actual = ea.run_regulation_enrichment(
        regulation_data=regulation_res,
        annotation=annotation,
        min_detected_in_set=1,
    ).reset_index(drop=True)

    expected = pd.DataFrame(
        {
            "terms": ["path3", "path1"],  # path2 has no peptides in foreground
            "identifiers": [
                # forground peptides concatenated by comma
                "gene6_pep1,gene9_pep1",
                "gene1_pep1,gene2_pep1,gene3_pep1",
            ],
            "foreground": [2, 3],
            "background": [0, 1],
            "foreground_pop": [5, 5],
            "background_pop": [8, 8],
            "pvalue": [0.4642857142857143, 1.0],
            "padj": [0.9285714285714286, 1.0],
            "rejected": [False, False],
        }
    )
    assert expected.equals(actual)


def test_run_up_down_regulation_enrichment_large():
    reg_df = pd.DataFrame.from_dict(
        {
            "leading_protein": {
                0: "prot1",
                1: "prot1",
                2: "prot1",
                3: "prot1",
                4: "prot1",
                5: "prot1",
                6: "prot2",
                7: "prot2",
                8: "prot2",
                9: "prot2",
                10: "prot2",
                11: "prot2",
                12: "prot2",
                13: "prot3",
                14: "prot3",
                15: "prot3",
                16: "prot3",
                17: "prot3",
                18: "prot3",
            },
            "padj": {
                0: 0.07880597878376583,
                1: 0.05990102335987204,
                2: 0.04416869805802626,
                3: 0.02462350123927702,
                4: 0.04416869805802626,
                5: 0.3168274635257489,
                6: 0.14417390797088134,
                7: 0.18583494940817125,
                8: 0.15951739880605717,
                9: 0.010309966852553444,
                10: 0.7131015411321412,
                11: 0.3745206530801302,
                12: 0.4850070137907755,
                13: 0.1022458686972483,
                14: 0.08669167880093272,
                15: 0.406777908314452,
                16: 0.35251153228723875,
                17: 0.02668745279847266,
                18: 0.01685484894274798,
            },
            "log2FC": {
                0: 1.5332710592773324,
                1: 3.0236567380723844,
                2: 2.911955797241669,
                3: 4.150443014180416,
                4: 3.5324234756094217,
                5: 1.2368642314001477,
                6: -2.2585678954343216,
                7: -8.54965231217265,
                8: 2.3163401651964417,
                9: 2.548886849610412,
                10: 10.897935053458644,
                11: -1.7617072017167672,
                12: -0.9393657125286836,
                13: -3.2248710758791166,
                14: 1.2705583236213729,
                15: 0.5485793599314555,
                16: -1.826508586620303,
                17: -1.180310353128581,
                18: -2.625593343445475,
            },
            "rejected": {
                0: False,
                1: False,
                2: True,
                3: True,
                4: True,
                5: False,
                6: False,
                7: False,
                8: False,
                9: True,
                10: False,
                11: False,
                12: False,
                13: False,
                14: False,
                15: False,
                16: False,
                17: True,
                18: True,
            },
            "group1": {
                0: "timepoint1",
                1: "timepoint1",
                2: "timepoint1",
                3: "timepoint1",
                4: "timepoint1",
                5: "timepoint1",
                6: "timepoint1",
                7: "timepoint1",
                8: "timepoint1",
                9: "timepoint1",
                10: "timepoint1",
                11: "timepoint1",
                12: "timepoint1",
                13: "timepoint1",
                14: "timepoint1",
                15: "timepoint1",
                16: "timepoint1",
                17: "timepoint1",
                18: "timepoint1",
            },
            "group2": {
                0: "timepoint2",
                1: "timepoint2",
                2: "timepoint2",
                3: "timepoint2",
                4: "timepoint2",
                5: "timepoint2",
                6: "timepoint2",
                7: "timepoint2",
                8: "timepoint2",
                9: "timepoint2",
                10: "timepoint2",
                11: "timepoint2",
                12: "timepoint2",
                13: "timepoint2",
                14: "timepoint2",
                15: "timepoint2",
                16: "timepoint2",
                17: "timepoint2",
                18: "timepoint2",
            },
            "identifier": {
                0: "prot1_pep1",
                1: "prot1_pep2",
                2: "prot1_pep3",
                3: "prot1_pep4",
                4: "prot1_pep5",
                5: "prot1_pep6",
                6: "prot2_pep1",
                7: "prot2_pep2",
                8: "prot2_pep3",
                9: "prot2_pep4",
                10: "prot2_pep5",
                11: "prot2_pep6",
                12: "prot2_pep7",
                13: "prot3_pep1",
                14: "prot3_pep2",
                15: "prot3_pep3",
                16: "prot3_pep4",
                17: "prot3_pep5",
                18: "prot3_pep6",
            },
        }
    )

    annotations = {
        "identifier": {
            0: "prot1",
            1: "prot1",
            2: "prot1",
            3: "prot1",
            4: "prot1",
            5: "prot2",
            6: "prot2",
            7: "prot2",
            8: "prot2",
            9: "prot2",
            10: "prot3",
            11: "prot3",
            12: "prot3",
            13: "prot3",
            14: "prot3",
            15: "prot3",
        },
        "source": {
            0: "Gene Ontology (biological process)",
            1: "Gene Ontology (cellular component)",
            2: "Catalytic activity",
            3: "Gene Ontology (biological process)",
            4: "Gene Ontology (biological process)",
            5: "Gene Ontology (biological process)",
            6: "Gene Ontology (cellular component)",
            7: "Catalytic activity",
            8: "Gene Ontology (biological process)",
            9: "Gene Ontology (biological process)",
            10: "Gene Ontology (biological process)",
            11: "Gene Ontology (cellular component)",
            12: "pH dependence",
            13: "Activity regulation",
            14: "Catalytic activity",
            15: "Gene Ontology (biological process)",
        },
        "annotation": {
            0: "pathway1",
            1: "pathway2",
            2: "pathway3",
            3: "pathway4",
            4: "pathway5",
            5: "pathway6",
            6: "pathway7",
            7: "pathway8",
            8: "pathway9",
            9: "pathway10",
            10: "pathway8",
            11: "pathway12",
            12: "pathway13",
            13: "pathway14",
            14: "pathway15",
            15: "pathway1",
        },
    }

    annotations = pd.DataFrame.from_dict(annotations)

    annotations_extended = reg_df[["leading_protein", "identifier"]].join(
        annotations.set_index("identifier"), on="leading_protein"
    )
    annotations_extended

    actual = ea.run_up_down_regulation_enrichment(
        regulation_data=reg_df,
        annotation=annotations_extended,
        identifier="identifier",
        pval_col="padj",
        min_detected_in_set=1,  # ! default is 2, so more conservative
        lfc_cutoff=0.1,  # ! the default is 1
    ).reset_index(drop=True)
    expected = pd.DataFrame.from_dict(
        {
            "terms": [
                "pathway2",
                "pathway3",
                "pathway4",
                "pathway5",
                "pathway8",
                "pathway1",
                "pathway10",
                "pathway6",
                "pathway7",
                "pathway9",
                "pathway12",
                "pathway13",
                "pathway14",
                "pathway15",
                "pathway1",
                "pathway8",
            ],
            "identifiers": [
                "prot1_pep3,prot1_pep4,prot1_pep5",
                "prot1_pep3,prot1_pep4,prot1_pep5",
                "prot1_pep3,prot1_pep4,prot1_pep5",
                "prot1_pep3,prot1_pep4,prot1_pep5",
                "prot2_pep4",
                "prot1_pep3,prot1_pep4,prot1_pep5",
                "prot2_pep4",
                "prot2_pep4",
                "prot2_pep4",
                "prot2_pep4",
                "prot3_pep5,prot3_pep6",
                "prot3_pep5,prot3_pep6",
                "prot3_pep5,prot3_pep6",
                "prot3_pep5,prot3_pep6",
                "prot3_pep5,prot3_pep6",
                "prot3_pep5,prot3_pep6",
            ],
            "foreground": [3, 3, 3, 3, 1, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2],
            "background": [3, 3, 3, 3, 12, 9, 6, 6, 6, 6, 4, 4, 4, 4, 10, 11],
            "foreground_pop": [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2],
            "background_pop": [
                19,
                19,
                19,
                19,
                19,
                19,
                19,
                19,
                19,
                19,
                19,
                19,
                19,
                19,
                19,
                19,
            ],
            "pvalue": [
                0.07094943240454077,
                0.07094943240454077,
                0.07094943240454077,
                0.07094943240454077,
                0.07094943240454077,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.08771929824561403,
                0.08771929824561403,
                0.08771929824561403,
                0.08771929824561403,
                0.5087719298245614,
                1.0,
            ],
            "padj": [
                0.14189886480908154,
                0.14189886480908154,
                0.14189886480908154,
                0.14189886480908154,
                0.14189886480908154,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.13157894736842105,
                0.13157894736842105,
                0.13157894736842105,
                0.13157894736842105,
                0.6105263157894737,
                1.0,
            ],
            "rejected": [
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
            ],
            "direction": [
                "upregulated",
                "upregulated",
                "upregulated",
                "upregulated",
                "upregulated",
                "upregulated",
                "upregulated",
                "upregulated",
                "upregulated",
                "upregulated",
                "downregulated",
                "downregulated",
                "downregulated",
                "downregulated",
                "downregulated",
                "downregulated",
            ],
            "comparison": [
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
                "timepoint1~timepoint2",
            ],
        }
    )

    actual.equals(expected)


def test_run_regulation_enrichment_with_duplicates():
    """Integration test for run_regulation_enrichment. Indirectly tests
    run_enrichment from enrichment_analysis module.

    should throw an error if unique identifier is enforced.
    """
    annotation = {
        "annotation": [
            "path1",
            "path1",
            "path1",
            "path2",
            "path2",  # e.g. protein 5 has two peptides associated
            "path2",  # e.g. protein 5 has two peptides associated
            "path3",
            "path3",
        ],
        "identifier": [
            "protein1",
            "protein2",
            "protein3",
            "protein1",
            "protein5",  # duplicated pep for protein identifier
            "protein5",  # duplicated pep for protein identifier
            "protein6",
            "protein9",
        ],
        "source": ["GO", "GO", "GO", "GO_P", "GO_P", "GO_P", "GO_P", "GO_P"],
    }
    annotation = pd.DataFrame(annotation)
    regulation_res = {
        "identifier": [
            "protein1",
            "protein2",
            "protein3",
            "protein4",
            "protein5",
            "protein5",
            "protein6",
        ],
        "rejected": [True, True, False, False, True, True, True],
    }
    regulation_res = pd.DataFrame(regulation_res)

    with pytest.raises(ValueError):
        _ = ea.run_regulation_enrichment(
            regulation_data=regulation_res,
            annotation=annotation,
            min_detected_in_set=1,
        )
