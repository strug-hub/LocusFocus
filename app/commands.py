from os import path
from pathlib import Path
import shutil
import tarfile
import zipfile

import click
from flask import Flask, current_app
import requests
import tqdm


def register_cli(app: Flask):
    @app.cli.command(name="download-smr-mqtl", help="Download mQTL files")
    @click.option(
        "--lite",
        is_flag=True,
        help="Download 'lite' version of the McRae et al. mQTL data (only SNPs with P < 1e-5 are included; 241 MB)",
    )
    def download_smr_mqtl(lite):
        # https://yanglab.westlake.edu.cn/software/smr/#mQTLsummarydata
        file_list = [
            # Whole blood mQTL data set used in Hannon et al. (2018 AJHG).(121MB)
            # Saved as US_mQTLS_SMR_format
            "https://yanglab.westlake.edu.cn/data/SMR/US_mQTLS_SMR_format.zip",
            # 42MB
            "https://yanglab.westlake.edu.cn/data/SMR/Hannon_Blood_dataset1.zip",
            # 25MB
            "https://yanglab.westlake.edu.cn/data/SMR/Hannon_Blood_dataset2.zip",
            # https://yanglab.westlake.edu.cn/data/SMR/Hannon_FetalBrain.zip (4.8MB)
            "https://yanglab.westlake.edu.cn/data/SMR/Hannon_FetalBrain.zip",
            # mQTL summary data from a meta-analysis of samples of East Asian ancestry. (2.5GB)
            # no particular tissue? saved as EAS
            "https://yanglab.westlake.edu.cn/data/SMR/EAS.tar.gz",
            # mQTL summary data from a meta-analysis of samples of European ancestry. (3.7GB)
            # no particular tissue? saved as EUR
            "https://yanglab.westlake.edu.cn/data/SMR/EUR.tar.gz",
            # Brain-mMeta mQTL summary data (Qi et al. 2018 Nat Commun) in SMR binary (BESD) format: Brain-mMeta.tar.gz (893 MB)
            # brain (from meta-analysis)
            "https://yanglab.westlake.edu.cn/data/SMR/Brain-mMeta.tar.gz",
        ] + (
            # Lite version of the McRae et al. mQTL data (only SNPs with P < 1e-5 are included; 241 MB)
            # peripheral blood
            ["https://yanglab.westlake.edu.cn/data/SMR/LBC_BSGS_meta_lite.tar.gz"]
            if lite
            else [
                # McRae et al. mQTL summary data (7.5 GB)
                "https://yanglab.westlake.edu.cn/data/SMR/LBC_BSGS_meta.tar.gz",
            ]
        )

        data_dir = path.join(path.dirname(path.dirname(__file__)), "data", "smr_mqtl")
        Path(data_dir).mkdir(exist_ok=True, parents=True)

        for file_url in file_list:
            filename = path.basename(file_url)
            _, ext = path.splitext(filename)
            tmp_save_path = path.join(data_dir, filename)
            with requests.get(file_url, stream=True) as r:
                current_app.logger.info(f"Downloading {filename}...")
                if r.status_code != 200:
                    r.raise_for_status()
                    raise RuntimeError(
                        f"Request to {file_url} returned status code {r.status_code}"
                    )
                file_size = int(r.headers.get("Content-Length", 0))
                desc = "(Unknown total file size)" if file_size == 0 else ""
                with tqdm.tqdm.wrapattr(
                    r.raw, "read", total=file_size, desc=desc
                ) as r_raw:
                    with open(tmp_save_path, "wb") as fd:
                        shutil.copyfileobj(r_raw, fd)

            if ext == ".zip":
                with zipfile.ZipFile(tmp_save_path) as zf:
                    zf.extractall(data_dir)
                    Path(tmp_save_path).unlink()
            elif ext == ".gz":
                with tarfile.open(tmp_save_path) as tf:
                    tf.extractall(data_dir, filter="data")
                    Path(tmp_save_path).unlink()
                    tf.close()
