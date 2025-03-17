from random import randint

from app import mongo
from app.utils.apis.gtex import get_variants


def test_can_fetch_v10_variants_from_region_string():
    """Sanity check for v10 variant fetch"""
    region_string = "chr11:0-200000"
    chr, pos = region_string.split(":")
    start, end = pos.split("-")
    results = get_variants(
        dataset_id="gtex_v10", start=int(start), end=int(end), chromosome=chr
    )

    assert len(results.data) > 0
    assert results.data[0].snp_id.startswith("rs")


def test_can_fetch_v8_variants_from_region_string(app):
    region_string = "chr11:0-200000"
    chr, pos = region_string.split(":")
    start, end = pos.split("-")
    with app.app_context():
        client = mongo.cx
        db = client.GTEx_V8
        collection = db["variant_table"]
        variants_query = collection.find(
            {
                "$and": [
                    {"chr": int(chr.replace("chr", ""))},
                    {"variant_pos": {"$gte": int(start), "$lte": int(end)}},
                ]
            }
        )
        variants_list_mongo = sorted(
            list(variants_query), key=lambda x: x["variant_pos"]
        )

    variants_api = get_variants(
        dataset_id="gtex_v8", start=int(start), end=int(end), chromosome=chr
    )

    variants_list_api = sorted(variants_api.data, key=lambda x: x.pos)

    assert len(variants_list_api) == len(variants_list_mongo)

    idx = randint(0, len(variants_list_api) - 1)

    assert (
        variants_list_api[idx].variant_id.replace("chr", "")
        == variants_list_mongo[idx]["variant_id"]
    )
