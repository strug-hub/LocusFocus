from typing import Dict, List

import requests


def fetch_variant_info(
    build: str, chr: str, location: int
) -> List[Dict]:
    """Fetch variant info from ensembl API

    :param build: The genome build
    :type build: str
    :param chr: The chromosome
    :type chr: str
    :param location: The location
    :type location: int
    :return: List of variants
    :rtype: List[Dict]

    ## Example Return value

    ```json
        [
            {
                "strand": 1,
                "feature_type": "variation",
                "id": "rs1873295566",
                "source": "dbSNP",
                "consequence_type": "intron_variant",
                "clinical_significance": [],
                "start": 120881839,
                "end": 120881839,
                "seq_region_name": "12",
                "alleles": [
                    "A",
                    "G"
                ],
                "assembly_name": "GRCh38"
            }
        ]
    ```
    """
    subd = "grch37." if build in ["hg19", "grch37"] else ""

    chr = f"chr{chr}" if not chr.startswith("chr") else chr

    query = f"https://{subd}rest.ensembl.org/overlap/region/homo_sapiens/{chr}:{location}-{location}?feature=variation"

    results = requests.get(query, headers={"accept":"application/json"})

    results.raise_for_status()

    return results.json()
