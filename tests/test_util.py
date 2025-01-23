import pytest

from requests.exceptions import HTTPError

from app.utils.util import format_ensembl_variant
from app.utils.apis.ensembl import fetch_ensembl_variant_info


def test_fetch_variant_info():
    # grch37
    res = fetch_ensembl_variant_info("grch37", "chr12", 120881837)
    assert len(res) == 1
    assert res[0]["id"].startswith("rs")
    assert isinstance(res[0]["alleles"], list)

    # grch38
    res = fetch_ensembl_variant_info("grch38", "chr12", 120881839)
    assert len(res) == 1
    assert res[0]["id"].startswith("rs")
    assert isinstance(res[0]["alleles"], list)

    # bad input (caller should catch HTTPError)
    with pytest.raises(HTTPError):
        res = fetch_ensembl_variant_info("grch38", "chr12", -27)


def test_format_variant():
    variant = {
        "alleles": ["A", "G"],
        "end": 120881839,
        "assembly_name": "GRCh38",
        "start": 120881839,
        "source": "dbSNP",
        "id": "rs1873295566",
        "feature_type": "variation",
        "consequence_type": "intron_variant",
        "clinical_significance": [],
        "seq_region_name": "12",
        "strand": 1,
    }

    formatted = format_ensembl_variant(variant, "abc")

    assert formatted == "12_120881839_A_G_abc"

    biallelic_variant = {
        "alleles": ["A", "G", "T"],
        "end": 120881839,
        "assembly_name": "GRCh38",
        "start": 120881839,
        "source": "dbSNP",
        "id": "rs1873295566",
        "feature_type": "variation",
        "consequence_type": "intron_variant",
        "clinical_significance": [],
        "seq_region_name": "12",
        "strand": 1,
    }

    assert format_ensembl_variant(biallelic_variant, "abc") == "12_120881839_A_G,T_abc"


def test_standardizeSNPs(app):
    with app.app_context():
        from app.routes import standardizeSNPs

        varlist = ["chr1_16896_GT_G_b38"]
        snps = standardizeSNPs(varlist, "chr1:16896-16896", "hg38")
        assert len(snps) == 1
        assert snps[0] == "1_16896_GT_G_b38"

        # biallelic
        varlist = ["chr12_120881839_A_G,T_b38"]
        snps = standardizeSNPs(varlist, "chr12:120881839-120881839", "hg38")
        assert len(snps) == 1
        assert snps[0] == "12_120881839_A_G,T_b38"

        # rsid (not in gtex)
        varlist = ["rs1873295566"]
        snps = standardizeSNPs(varlist, "chr12:120881839-120881839", "hg38")
        assert len(snps) == 1
        assert snps[0] == "."

        # rsid (in gtex) (todo, one of these needs to work with db before porting over)
        # problem is my local variant db is evidently misnamed (has hg37 but labeled 38)
        # anyway, it seems 7 is old anyway
        # in any case, might want to use ensemble rest API to convert from gh19 to 38
        # but this is failing b/c w/ those that start w/ rs, it uses the lookup, which has
        # the wrong column name for rs id
        # however, with the gtex API we can query on rsid and get the b37 and 38 ids...
        # yeah, so can replace that lookup, then grab the right ID based on build
        # and this removes it!
        varlist = ["rs533210981"]
        snps = standardizeSNPs(varlist, "chr1:16911-16911", "hg38")
        assert len(snps) == 1
        assert snps[0] == "1_16911_CA_C_b37"

