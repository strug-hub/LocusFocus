import pytest
import pandas as pd

from app.utils.helpers import validate_chromosome, adjust_snp_column


def test_validate_chromosome():
    assert validate_chromosome(chr=1, prefix=None)
    assert validate_chromosome(chr="chr1", prefix="chr")
    assert validate_chromosome(chr=23, prefix=None, x_y_numeric=True)
    assert validate_chromosome(chr="chX", prefix="ch", x_y_numeric=False)

    with pytest.raises(ValueError):
        validate_chromosome(chr=1, prefix="ch")
    with pytest.raises(ValueError):
        validate_chromosome(chr="chr1", prefix="None")
    with pytest.raises(ValueError):
        validate_chromosome(chr=23, prefix=None, x_y_numeric=False)
    with pytest.raises(ValueError):
        validate_chromosome(chr="chX", prefix="ch", x_y_numeric=True)


def test_adjust_snp_column_rsid():
    """Verify that rsIDs are ignored"""
    snps_df = pd.DataFrame(
        {
            "SNP": ["rs12345", "rs12346", "rs12347"],
            "CHROM": ["1", "2", "3"],
            "POS": ["100", "200", "300"],
            "REF": ["A", "C", "G"],
            "ALT": ["T", "G", "C"],
        }
    )

    expected = snps_df.copy()
    actual = adjust_snp_column(snps_df, target_build="hg19")

    assert expected.equals(actual)


def test_adjust_snp_column_chr_pos_ref_alt_build():
    """Verify that chr_pos_ref_alt_build format is maintained"""
    snps_df = pd.DataFrame(
        {
            "SNP": ["1_100_A_T_b37", "2_200_C_G_b37", "3_300_G_C_b37"],
            "CHROM": ["1", "2", "3"],
            "POS": ["100", "200", "300"],
            "REF": ["A", "C", "G"],
            "ALT": ["T", "G", "C"],
        }
    )

    expected = pd.DataFrame(
        {
            "SNP": ["1_100_A_T_b38", "2_200_C_G_b38", "3_300_G_C_b38"],
            "CHROM": ["1", "2", "3"],
            "POS": ["100", "200", "300"],
            "REF": ["A", "C", "G"],
            "ALT": ["T", "G", "C"],
        }
    )
    actual = adjust_snp_column(snps_df, target_build="hg38")

    assert expected.equals(actual)


def test_adjust_snp_column_bad_input():
    """Verify that an error is raised for invalid input"""
    snps_df = pd.DataFrame(
        {
            "SNP": ["1_100_A_T_b37", "2_200_C_G_b37", "3_300_G_C_b37"],
            "CHROM": ["1", "2", "3"],
            "POS": ["100", "200", "300"],
            "REF": ["A", "C", "G"],
            "ALT": ["T", "G", "C"],
        }
    )

    with pytest.raises(ValueError):
        adjust_snp_column(snps_df, target_build="19")
    with pytest.raises(ValueError):
        adjust_snp_column(snps_df, target_build="b37")


def test_adjust_snp_column_mixture():
    """Verify that a mixture of rsIDs and chr_pos_ref_alt_build format is maintained"""
    snps_df = pd.DataFrame(
        {
            "SNP": ["rs12345", "2_200_C_G_b37", "3_300_G_C_b37", "4_400_T_A_b37"],
            "CHROM": ["1", "2", "3", "4"],
            "POS": ["100", "200", "300", "400"],
            "REF": ["A", "C", "G", "T"],
            "ALT": ["T", "G", "C", "A"],
        }
    )

    expected = pd.DataFrame(
        {
            "SNP": ["rs12345", "2_200_C_G_b38", "3_300_G_C_b38", "4_400_T_A_b38"],
            "CHROM": ["1", "2", "3", "4"],
            "POS": ["100", "200", "300", "400"],
            "REF": ["A", "C", "G", "T"],
            "ALT": ["T", "G", "C", "A"],
        }
    )
    actual = adjust_snp_column(snps_df, target_build="hg38")

    assert expected.equals(actual)


def test_adjust_snp_column_ignore_alleles():
    """Verify that alleles are ignored"""
    snps_df = pd.DataFrame(
        {
            "SNP": ["rs12345", "2_200_C_G_b37", "3_300_G_C_b37", "4_400_T_A_b37"],
            "CHROM": ["1", "2", "3", "4"],
            "POS": ["100", "200", "300", "400"],
        }
    )

    expected = pd.DataFrame(
        {
            "SNP": ["rs12345", "2_200_C_G_b38", "3_300_G_C_b38", "4_400_T_A_b38"],
            "CHROM": ["1", "2", "3", "4"],
            "POS": ["100", "200", "300", "400"],
        }
    )

    actual = adjust_snp_column(snps_df, target_build="hg38", ignore_alleles=True)

    assert expected.equals(actual)
