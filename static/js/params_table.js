const START = 0;
const END = 1;

function buildParamsTable(data, sessionid, type = "default") {
  var tableselect = d3.select("#params-table");

  // Column headers
  header = tableselect.append("thead").append("tr");
  header.append("th").attr("class", "th-sm").text("Field");
  header.append("th").attr("class", "th-sm").text("Value");

  // Table body:
  let tbody = tableselect.append("tbody");
  // field-value pairs in order displayed in table
  let table_data = [];
  if (type === "default") {
    table_data = [
      ["Session ID", sessionid],
      ["Lead SNP", data["lead_snp"]],
      ["Chromosome", data["chrom"]],
      ["Start position", data["startbp"]],
      ["End position", data["endbp"]],
      ["Build", data["coordinate"]],
      ["Infer variants", data["inferVariant"]],
      [
        `Number of SNPs in ${data["chrom"]}:${data["startbp"]}-${data["endbp"]}`,
        data["snps"].length,
      ],
      ["LD Population", data["ld_populations"]],
      ["GTEx version", data["gtex_version"]],
      ["Number of GTEx tissues selected", data["gtex_tissues"].length],
      [
        "Number of GTEx genes selected",
        data["gtex_genes"] !== undefined
          ? data["gtex_genes"].length
          : undefined,
      ],
      ["SS region", data["SS_region"]],
      [
        `Number of SNPs in ${data["chrom"]}:${data["SS_region"][START]}-${data["SS_region"][END]}`,
        data["num_SS_snps"],
      ],
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
  } else if (type === "set-based-test") {
    table_data = [
      ["Session ID", sessionid],
      ["Chromosome", data["chrom"]],
      ["Start position", data["startbp"]],
      ["End position", data["endbp"]],
      ["Build", data["coordinate"]],
      [
        `Number of SNPs in ${data["chrom"]}:${data["startbp"]}-${data["endbp"]}`,
        data["snps"].length,
      ],
      ["LD Population", data["ld_populations"]],
      ["First stage -log10(SS P-value) threshold", data["set_based_p"]],
      [
        "Number of user-provided secondary datasets",
        data["secondary_dataset_titles"].length,
      ],
    ];
  } else {
    // Shouldn't get here
    throw Error(`Unexpected params table type: ${type}`);
  }

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
