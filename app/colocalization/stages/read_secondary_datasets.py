from bs4 import BeautifulSoup as bs
import numpy as np

import app.utils.htmltableparser as htmltableparser

from app.colocalization.payload import SessionPayload
from app.utils import get_file_with_ext
from app.pipeline.pipeline_stage import PipelineStage
from app.utils.errors import InvalidUsage


class ReadSecondaryDatasetsStage(PipelineStage):
    """
    Stage responsible for reading user-uploaded secondary dataset files
    and storing them in the session payload, if they exist.
    """

    def name(self):
        return "read-secondary-datasets"

    def description(self) -> str:
        return "Read secondary datasets (if any)"

    def invoke(self, payload: SessionPayload) -> SessionPayload:
        secondary_datasets = self._read_dataset_file(payload)
        payload.secondary_datasets = secondary_datasets

        return payload

    def _read_dataset_file(self, payload: SessionPayload):
        """
        Check if the secondary dataset file exists, and then parse the HTML for dataset tables.
        """

        html_filepath = get_file_with_ext(payload.uploaded_files, ["html"])

        if html_filepath is None or html_filepath == "":
            return None

        secondary_datasets = dict()

        with open(html_filepath, encoding="utf-8", errors="replace") as f:
            html = f.read()
            if (
                (not html.startswith("<h3>"))
                and (not html.startswith("<html>"))
                and (
                    not html.startswith("<table>")
                    and (not html.startswith("<!DOCTYPE html>"))
                )
            ):
                raise InvalidUsage(
                    "Secondary dataset(s) provided are not formatted correctly. Please use the merge_and_convert_to_html.py script for formatting.",
                    status_code=410,
                )
        soup = bs(html, "lxml")
        table_titles = soup.find_all("h3")
        table_titles = [x.text for x in table_titles]
        tables = soup.find_all("table")
        hp = htmltableparser.HTMLTableParser()
        for i in np.arange(len(tables)):
            try:
                table = hp.parse_html_table(tables[i])
                secondary_datasets[table_titles[i]] = table.fillna(-1).to_dict(
                    orient="records"
                )
            except Exception:
                secondary_datasets[table_titles[i]] = []

        return secondary_datasets
