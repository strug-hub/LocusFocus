function buildParamsTable(data, sessionid) {
  var tableselect = d3.select("#params-table");

  // Column headers
  header = tableselect.append("thead").append("tr");
  header.append("th").attr("class", "th-sm").text("Field");
  header.append("th").attr("class", "th-sm").text("Value");

  chrom = data["chrom"];
  startbp = data["startbp"];
  endbp = data["endbp"];
  SS_start = data["SS_region"][0];
  SS_end = data["SS_region"][1];

  // Table body:
  let tbody = tableselect.append("tbody");
  // field-value pairs in order displayed in table
  let table_data = [
    ["Session ID", sessionid],
    ["Lead SNP", data["lead_snp"]],
    ["Chromosome", chrom],
    ["Start position", startbp],
    ["End position", endbp],
    ["Build", data["coordinate"]],
    ["Infer variants", data["inferVariant"]],
    [`Number of SNPs in ${chrom}:${startbp}-${endbp}`, data["snps"].length],
    ["LD Population", data["ld_populations"]],
    ["GTEx version", data["gtex_version"]],
    ["Number of GTEx tissues selected", data["gtex_tissues"].length],
    [
      "Number of GTEx genes selected",
      data["gtex_genes"] !== undefined ? data["gtex_genes"].length : undefined,
    ],
    ["SS region", data["SS_region"]],
    [`Number of SNPs in ${chrom}:${SS_start}-${SS_end}`, data["num_SS_snps"]],
    ["First stage -log10(SS P-value) threshold", data["set_based_p"]],
    ["Many SNPs not matching GTEx SNPs", data["snp_warning"]],
    ["SNPs matching threshold level", data["thresh"]],
    ["Number of SNPs matching with GTEx", data["numGTExMatches"]],
    [
      "Number of user-provided secondary datasets",
      data["secondary_dataset_titles"].length,
    ],
    ["Run COLOC2", data["runcoloc2"]],
  ];

  for (const [key, value] of table_data) {
    if (value === undefined) continue;
    let row = tbody.append("tr");
    row.append("td").text(key);
    row.append("td").text(value);
  }

  // Add DataTables functionality:
  paramsTable = $(document).ready(function () {
    $("#params-table").DataTable({
      dom: "Bfrtipl",
      buttons: [
        "copy",
        {
          extend: "csv",
          filename: "parameters_table",
        },
        {
          extend: "excel",
          filename: "parameters_table",
          messageTop: "Input parameters",
        },
      ],
    });
  });
}
