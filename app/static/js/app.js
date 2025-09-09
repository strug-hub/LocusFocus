const startingChr = "1";
const startingPos = "205,500,000";
const endingPos = "206,000,000";
const genomicWindowLimit = 2e6;
var submitButton = d3.select("#submit-btn");
var errorDiv = d3.select("#error-messages");
var theTable = d3.select("#variants-table");
var gtexTissuesMsgDiv = d3.select("#tissue-select");
var coordinate = "hg19"; // default
var gtex_version = "v7"; // default
var gtexurl = `/gtex/${gtex_version}/tissues_list`;
var markerColDiv = d3.select("#snp");
var variantInputsDiv = d3.select("#variantInputs");
var chromColDiv = d3.select("#chrom");
var posColDiv = d3.select("#pos");
var refColDiv = d3.select("#ref");
var altColDiv = d3.select("#alt");
var statsDiv = d3.select("#statsDiv");
var statsDiv2 = d3.select("#statsDiv2");

const formFields = {
  locus: {
    selector: "#locus",
    type: "input",
  },
  snpCol: {
    selector: "#snp-col",
    type: "input",
  },
  ssLocus: {
    selector: "#SSlocus",
    type: "input",
  },
  leanSnp: {
    selector: "#leadsnp",
    type: "input",
  },
  setBasedP: {
    selector: "#setbasedP",
    type: "input",
  },
  coordinate: {
    selector: "#coordinate",
    type: "input",
  },
  chrom: {
    selector: "#chrom input",
    type: "input",
  },
  pos: {
    selector: "#pos input",
    type: "input",
  },
  ref: {
    selector: "#ref input",
    type: "input",
  },
  alt: {
    selector: "#alt input",
    type: "input",
  },
  ldPop: {
    selector: "#LD-populations",
    type: "input",
  },
  getxVersion: {
    selector: "#GTEx-version",
    type: "select",
  },
  gtexTissues: {
    selector: "#GTEx-tissues",
    type: "multiselect",
  },
  regionGenes: {
    selector: "#region-genes",
    type: "multiselect",
  },
};

function storeFormArgs(e) {
  const vals = Object.entries(formFields)
    .map(([k, v]) => {
      if (v.type === "input") {
        return { [k]: d3.select(v.selector).node().value };
      } else if (v.type === "select") {
        return { [k]: d3.select(v.selector).node().value };
      } else if (v.type === "multiselect") {
        return { [k]: $(v.selector).val() };
      }
    })
    .reduce((acc, curr) => ({
      ...acc,
      ...curr,
    }));

  localStorage.setItem("formVals", JSON.stringify(vals));
}

const addLoadingIndicator = () => {
  $("#loading-modal").modal({ closeClass: "icon-remove" });
};

const closeActiveModal = () => $.modal.close();

const setErrorModalOpen = (message) => {
  d3.select("#error-modal p.message").text(message);
  $("#error-modal").modal();
};

const setPendingModalOpen = () => {
  $("#job-pending-modal").modal({ closeClass: "icon-remove" });
};

const setRunningModalOpen = (stage_index, stage_count) => {
  $("#job-running-modal").modal({ closeClass: "icon-remove" });
  d3.select("#job-running-modal h4.stage").text(
    `Currently on stage ${stage_index} of ${stage_count}`
  );
};

d3.select(".inputs-and-upload-form form").on("submit", storeFormArgs);

$(".inputs-and-upload-form form").on("submit", async (e) => {
  e.preventDefault();

  addLoadingIndicator();

  try {
    const res = await fetch("/", {
      method: "POST",
      body: new FormData(e.currentTarget),
    });

    if (!res.status.toString().startsWith("2")) {
      throw res;
    }

    const content = await res.json();

    if (!content.queued) {
      window.location = `${window.location}/session_id/${content.session_id}`;
    } else {
      handleJobStatus(`/job/status/${content.session_id}`, content.session_id);
    }
  } catch (e) {
    closeActiveModal();
    console.error(e);
    if (e?.status.toString().startsWith("4")) {
      const error = await e.json();
      const message = error.message;
      setErrorModalOpen(message);
    } else {
      setErrorModalOpen(
        "The job failed due to an unexpected error, please try again later."
      );
    }
  }
});

/**
 * Handle checking the status of a job and
 * updating the progress bar on the waiting page.
 */
async function handleJobStatus(jobStatusURL, sessionId) {
  const statusResponse = await fetch(jobStatusURL);
  const statusData = await statusResponse.json();
  const jobStatus = statusData.status;
  if (jobStatus == "PENDING") {
    closeActiveModal();
    setPendingModalOpen();
    await new Promise((resolve) =>
      setTimeout(() => resolve(handleJobStatus(jobStatusURL, sessionId)), 10000)
    );
  }

  if (jobStatus == "RUNNING") {
    const stage_index = statusData.stage_index + 1;
    const stage_count = statusData.stage_count;
    closeActiveModal();
    setRunningModalOpen(stage_index, stage_count);
    await new Promise((resolve) =>
      setTimeout(() => resolve(handleJobStatus(jobStatusURL, sessionId)), 500)
    );
  }

  if (jobStatus == "FAILURE") {
    setErrorModalOpen(
      `Error: ${statusData.error_title}; Details: ${statusData.error_message}`
    );
  }

  if (jobStatus == "SUCCESS") {
    const { redirect_url } = statusData;
    window.location = redirect_url;
  }
}

var locText = d3.select("#locusText").empty()
  ? null
  : d3.select("#locusText").text();
d3.select("#locusText").text(
  `${locText} (max: ${genomicWindowLimit / 1e6} Mbp):`
);
d3.select("#locus").attr("value", `${startingChr}:${startingPos}-${endingPos}`);

// FUNCTIONS

function setLoading(loading) {
  d3.select("#loading").style("display", loading ? "block" : "none");
}

// GTEx version selection change
function gtexVersionChange(newVersion) {
  gtex_version = newVersion.toLowerCase();
  gtexTissuesMsgDiv.text(`Select GTEx (${gtex_version.toUpperCase()}) Tissues`);
}

function updateGtexVersionsAllowed(newCoordinate) {
  if (newCoordinate.toLowerCase() === "hg19") {
    $("#gtex-v7").prop("disabled", false);
    $("#gtex-v8").prop("disabled", true);

    $("#gtex-v7").prop("selected", true);
    $("#gtex-v8").prop("selected", false);
    gtex_version = "v7";
  } else {
    $("#gtex-v7").prop("disabled", true);
    $("#gtex-v8").prop("disabled", false);

    $("#gtex-v7").prop("selected", false);
    $("#gtex-v8").prop("selected", true);
    gtex_version = "v8";
  }
  gtexTissuesMsgDiv.text(`Select GTEx (${gtex_version.toUpperCase()}) Tissues`);
}

// Coordinate system selection change
function coordinateChange(newCoordinate) {
  $("#LD-populations").multiselect("destroy");
  $("#GTEx-tissues").multiselect("destroy");
  $("#region-genes").multiselect("destroy");
  $("#GTEx-version").multiselect("destroy");
  d3.select("#locus").property(
    "value",
    `${startingChr}:${startingPos}-${endingPos}`
  );
  if (!["hg38", "hg19"].includes(newCoordinate)) {
    alert("Invalid coordinate system. Please select hg38 or hg19 coordinates.");
    return;
  }
  setLoading(true);
  if (newCoordinate === "hg38") {
    gtex_version = "v8";
    gtexurl = `/gtex/${gtex_version}/tissues_list`;
    coordinate = "hg38";
    d3.select("#genes-select").text(
      "Select Genes (enter coordinates above to populate)"
    );
  } else if (newCoordinate.toLowerCase() == "hg19") {
    gtex_version = "v7";
    gtexurl = `/gtex/${gtex_version}/tissues_list`;
    coordinate = "hg19";
    d3.select("#genes-select").text(
      "Select Genes (enter coordinates above to populate)"
    );
  }
  gtexTissuesMsgDiv.text(`Select GTEx (${gtex_version.toUpperCase()}) Tissues`);
  updateGtexVersionsAllowed(newCoordinate);
  d3.select("#region-genes").text("");
  init();
}

function askBetaInput(betaColDiv) {
  betaColDiv.html("");
  betaColDiv
    .append("label")
    .attr("for", "betaColname")
    .attr("data-html", "true")
    .attr("data-toggle", "tooltip")
    .attr(
      "title",
      "Enter the header text corresponding to the beta column in your txt/tsv file (primary dataset)"
    )
    .text("Beta Column Name:");
  betaColDiv
    .append("input")
    .attr("class", "form-control")
    .attr("name", "beta-col")
    .attr("type", "text")
    .attr("value", "BETA")
    .attr("onfocus", "this.value=''")
    .attr("data-toggle", "tooltip")
    .attr("data-html", "true")
    .attr(
      "title",
      "Enter the header text corresponding to the beta column in your txt/tsv file (primary dataset)"
    );
}
function askStdErrInput(stderrColDiv) {
  stderrColDiv.html("");
  stderrColDiv
    .append("label")
    .attr("for", "stderrColname")
    .attr("data-html", "true")
    .attr("data-toggle", "tooltip")
    .attr(
      "title",
      "Enter the header text corresponding to the standard error column in your txt/tsv file (primary dataset)"
    )
    .text("Standard Error Column Name:");
  stderrColDiv
    .append("input")
    .attr("class", "form-control")
    .attr("name", "stderr-col")
    .attr("type", "text")
    .attr("value", "SE")
    .attr("onfocus", "this.value=''")
    .attr("data-toggle", "tooltip")
    .attr("data-html", "true")
    .attr(
      "title",
      "Enter the header text corresponding to the standard error column in your txt/tsv file (primary dataset)"
    );
}
function askNumSamplesInput(numSamplesDiv) {
  numSamplesDiv.html("");
  numSamplesDiv
    .append("label")
    .attr("for", "numSamples")
    .attr("data-html", "true")
    .attr("data-toggle", "tooltip")
    .attr(
      "title",
      "Enter the header text corresponding to the number of samples column in your txt/tsv file (primary dataset)"
    )
    .text("Number of Samples Column Name:");
  numSamplesDiv
    .append("input")
    .attr("class", "form-control")
    .attr("name", "numsamples-col")
    .attr("id", "numsamples-col")
    .attr("type", "text")
    .attr("value", "N")
    .attr("onfocus", "this.value=''")
    .attr("data-toggle", "tooltip")
    .attr("data-html", "true")
    .attr(
      "title",
      "Enter the header text corresponding to the number of samples column in your txt/tsv file (primary dataset)"
    );
}
function askPvalueInput(pvalueColDiv) {
  pvalueColDiv.html("");
  pvalueColDiv
    .append("label")
    .attr("for", "p-value")
    .attr("data-html", "true")
    .attr("data-toggle", "tooltip")
    .attr(
      "title",
      "Header text corresponding to the p-value column in your txt/tsv file (primary dataset)"
    )
    .text("P-value Column Name:");
  pvalueColDiv
    .append("input")
    .attr("class", "form-control")
    .attr("name", "pval-col")
    .attr("type", "text")
    .attr("value", "P")
    .attr("onfocus", "this.value=''")
    .attr("data-toggle", "tooltip")
    .attr("data-html", "true")
    .attr(
      "title",
      "Enter the header text corresponding to the p-value column in your txt/tsv file (primary dataset)"
    );
}
function askMafInput(mafColDiv) {
  mafColDiv.html("");
  mafColDiv
    .append("label")
    .attr("for", "maf")
    .attr("data-html", "true")
    .attr("data-toggle", "tooltip")
    .attr(
      "title",
      "Header text corresponding to the MAF column in your txt/tsv file (primary dataset)"
    )
    .text("MAF Column Name:");
  mafColDiv
    .append("input")
    .attr("class", "form-control")
    .attr("name", "maf-col")
    .attr("type", "text")
    .attr("value", "MAF")
    .attr("onfocus", "this.value=''")
    .attr("data-toggle", "tooltip")
    .attr("data-html", "true")
    .attr(
      "title",
      "Enter the header text corresponding to the MAF column in your txt/tsv file (primary dataset)"
    );
}
function askNumCasesInput(studytype) {
  var numCasesDiv = d3.select("#numCases");
  numCasesDiv.html("");
  if (studytype === "cc") {
    numCasesDiv
      .append("label")
      .attr("for", "numcases")
      .attr("data-html", "true")
      .attr("data-toggle", "tooltip")
      .attr("title", "Enter the number of cases in the study")
      .text("Number of Cases:");
    numCasesDiv
      .append("input")
      .attr("class", "form-control")
      .attr("name", "numcases")
      .attr("id", "numcases")
      .attr("type", "text")
      .attr("value", 100)
      .attr("onfocus", "this.value=''")
      .attr("data-toggle", "tooltip")
      .attr("data-html", "true")
      .attr("title", "Enter the number of cases in the study")
      .attr("onchange", "checkNumSamplesInput(this.value)");
    numCasesDiv
      .append("div")
      .attr("class", "input-error")
      .attr("id", "numSamplesError-message");
  }
}
function askStudytypeInput(studytypeDiv) {
  studytypeDiv.html("");
  studytypeDiv
    .append("label")
    .attr("for", "studytype")
    .attr("data-toggle", "tooltip")
    .attr(
      "title",
      "Select whether the phenotype is quantitative, or whether  it is a case-control design"
    )
    .text("Select Study Type:");
  studytypeDivSelect = studytypeDiv
    .append("p")
    .append("select")
    .attr("id", "studytype")
    .attr("name", "studytype")
    .attr("onchange", "askNumCasesInput(this.value)");
  studytypeDivSelect
    .append("option")
    .text("Quantitative")
    .property("value", "quant");
  studytypeDivSelect
    .append("option")
    .text("Case-control")
    .property("value", "cc");
}
function askColocInputs() {
  statsDiv.html("");
  statsDiv2.html("");
  var betaColDiv = statsDiv
    .append("div")
    .attr("class", "col-md-3")
    .attr("id", "beta");
  var stderrColDiv = statsDiv
    .append("div")
    .attr("class", "col-md-3")
    .attr("id", "stderr");
  var numSamplesDiv = statsDiv
    .append("div")
    .attr("class", "col-md-3")
    .attr("id", "numSamples");
  var pvalueColDiv = statsDiv
    .append("div")
    .attr("class", "col-md-3")
    .attr("id", "p_value");
  var mafColDiv = statsDiv2
    .append("div")
    .attr("class", "col-md-3")
    .attr("id", "maf");
  var studytypeDiv = statsDiv2
    .append("div")
    .attr("class", "col-md-3")
    .attr("id", "studytype");
  statsDiv2.append("div").attr("class", "col-md-3").attr("id", "numCases");
  askBetaInput(betaColDiv);
  askStdErrInput(stderrColDiv);
  askNumSamplesInput(numSamplesDiv);
  askPvalueInput(pvalueColDiv);
  askMafInput(mafColDiv);
  askStudytypeInput(studytypeDiv);
}

function inferVariant(snpbox) {
  if (snpbox.checked) {
    $("#variantInputs").hide();
    ["#chrom-col", "#pos-col", "#ref-col", "#alt-col"].forEach((id) =>
      $(id).prop("disabled", true)
    );
  } else {
    $("#variantInputs").show();
    ["#chrom-col", "#pos-col", "#ref-col", "#alt-col"].forEach((id) =>
      $(id).prop("disabled", false)
    );
  }
  // re-initialize popover and tooltip
  $(function () {
    $('[data-toggle="popover"]').popover();
  });
  $(document).ready(function () {
    $('[data-toggle="tooltip"]').tooltip();
  });
}

function addColoc2Inputs(colocinput) {
  if (colocinput.checked) {
    askColocInputs();
  } else {
    statsDiv.html("");
    statsDiv2.html("");
    var pvalueColDiv = statsDiv
      .append("div")
      .attr("class", "col-md-3")
      .attr("id", "p_value");
    askPvalueInput(pvalueColDiv);
  }
  // re-initialize popover and tooltip
  $(function () {
    $('[data-toggle="popover"]').popover();
  });
  $(document).ready(function () {
    $('[data-toggle="tooltip"]').tooltip();
  });
}

function loadGenes(build, region) {
  var chrom = parseInt(
    region.split(":")[0].toLowerCase().replace("chr", "").replace("x", "23")
  );
  var startbp = parseInt(
    region.split(":")[1].split("-")[0].replaceAll(",", "")
  );
  var endbp = parseInt(region.split(":")[1].split("-")[1].replaceAll(",", ""));
  var genesdiv = d3.select("#region-genes");
  var genesMsgDiv = d3.select("#genes-select");
  d3.json(`/genenames/${build}/${chrom}/${startbp}/${endbp}`).then(
    (response) => {
      genesdiv.text("");
      var genenames = response.map((k) => k);
      if (genenames.length === 0) {
        genesMsgDiv.text(`No Genes Found in ${region} (${build})`);
      } else {
        genesMsgDiv.text(`Select Genes Found in ${region} (${build})`);
      }
      for (var i = 0; i < genenames.length; i++) {
        genesdiv
          .append("option")
          .attr("value", genenames[i])
          .text(genenames[i]);
      }
      $("#region-genes").multiselect("destroy");
      $(document).ready(function () {
        $("#region-genes").multiselect({
          enableFiltering: true,
          includeSelectAllOption: true,
          maxHeight: 400,
          buttonWidth: "400px",
          checkboxName: function (option) {
            return "multiselect[]";
          },
        });
      });
    }
  );
}

// Checks the database connection, and updates error messages appropriately
async function checkDB() {
  const response = await (await fetch("/dbstatus")).json();
  const errorSpan = errorDiv.select("#database-status-error-message");
  if (response.status !== "ok") {
    if (errorSpan.empty()) {
      errorDiv
        .append("span")
        .attr("id", "database-status-error-message")
        .text("Database connection has been lost. Please try again later.");
    }
  } else {
    if (!errorSpan.empty()) {
      errorSpan.remove();
    }
  }
}

function handleLDPopulations() {
  var lddiv = d3.select("#LD-populations");
  if (coordinate.toLowerCase() === "hg19") {
    popcodes = ["EUR", "AFR", "AMR", "ASN"];
    popdesc = [
      "1000 Genomes 2012 EUR",
      "1000 Genomes 2012 AFR",
      "1000 Genomes 2012 AMR",
      "1000 Genomes 2012 ASN",
    ];
    lddiv.text("");
    for (var i = 0; i < popcodes.length; i++) {
      lddiv.append("option").text(popdesc[i]).property("value", popcodes[i]);
    }
    lddiv.property("selectedIndex", 0); // EUR selected by default
  } else if (coordinate.toLowerCase() === "hg38") {
    popcodes = ["EUR", "AFR", "AMR", "EAS", "SAS", "NFE"];
    popdesc = [
      "1000 Genomes 2018 EUR",
      "1000 Genomes 2018 AFR",
      "1000 Genomes 2018 AMR",
      "1000 Genomes 2018 EAS",
      "1000 Genomes 2018 SAS",
      "1000 Genomes 2018 NFE",
    ];
    lddiv.text("");
    for (var i = 0; i < popcodes.length; i++) {
      lddiv.append("option").text(popdesc[i]).property("value", popcodes[i]);
    }
    lddiv.property("selectedIndex", 0); // EUR selected by default
  }
}

function updateChrXWarning() {
  const regionText = $("#locus").val().toLowerCase();
  const selectedAssembly = $("#coordinate").val();
  const uploadedFileNames = $.map(
    $("#file-upload").prop("files"),
    (f) => f.name
  );

  console.log(
    `[updateChrXWarning]
    regionText: ${regionText}
    selectedAssembly: ${selectedAssembly}
    uploadedFileNames: ${uploadedFileNames}`
  );

  const hideWarning = !(
    ["x", "23", "chrx"].some((start) => regionText.startsWith(start)) && // x chromosome
    selectedAssembly.toLowerCase() === "hg38" && // hg38
    uploadedFileNames.every((name) => !name.toLowerCase().endsWith(".ld"))
  ); // no LD provided

  $("#chrX-warning").prop("hidden", hideWarning);
}

function init() {
  String.prototype.replaceAll = function (search, replacement) {
    var target = this;
    return target.split(search).join(replacement);
  };

  // Check the database connection
  checkDB();
  setInterval(checkDB, 5 * 60_000);

  // // askSNPInput(markerColDiv);
  // askChromInput(chromColDiv);
  // askPosInput(posColDiv);
  // askRefInput(refColDiv);
  // askAltInput(altColDiv);

  statsDiv.html("");
  var pvalueColDiv = statsDiv
    .append("div")
    .attr("class", "col-md-3")
    .attr("id", "p_value");
  askPvalueInput(pvalueColDiv);

  // Build LD population selections depending on coordinate build chosen:
  handleLDPopulations();

  // Auto-complete for genes field:
  // d3.json(`/genenames/${coordinate}`).then(response => {
  //     $( "#gencodeID" ).autocomplete({
  //         source: response
  //     });
  // });

  loadGenes(coordinate, "1:205,500,000-206,000,000");

  d3.json(gtexurl).then((response) => {
    var gtex_tissues = response.map((k) => k);

    // Build GTEx tissues multiselect dropdown
    var gtexdiv = d3.select("#GTEx-tissues");
    gtexdiv.text("");
    for (var i = 0; i < gtex_tissues.length; i++) {
      gtexdiv
        .append("option")
        .attr("value", gtex_tissues[i])
        .text(gtex_tissues[i].replaceAll("_", " "));
    }

    setLoading(false);

    // Multi-Select Initialization
    $(document).ready(function () {
      //     $('#LD-populations').multiselect({
      //         enableClickableOptGroups: true,
      //         maxHeight: 400,
      //         buttonWidth: '400px',
      //         checkboxName: function(option) {
      //             return 'multiselect[]';
      //         }
      //     });
      $("#coordinate").multiselect({
        buttonWidth: "200px",
        checkboxName: function (option) {
          return "multiselect[]";
        },
      });
      $("#marker-format").multiselect({
        buttonWidth: "400px",
        checkboxName: function (option) {
          return "multiselect[]";
        },
      });
      $("#LD-populations").multiselect({
        buttonWidth: "400px",
        checkboxName: function (option) {
          return "multiselect[]";
        },
      });
      $("#GTEx-tissues").multiselect({
        enableCaseInsensitiveFiltering: true,
        includeSelectAllOption: true,
        maxHeight: 400,
        buttonWidth: "400px",
        checkboxName: function (option) {
          return "multiselect[]";
        },
      });
      $("#region-genes").multiselect({
        enableCaseInsensitiveFiltering: true,
        includeSelectAllOption: true,
        maxHeight: 400,
        buttonWidth: "400px",
        checkboxName: function (option) {
          return "multiselect[]";
        },
      });
      $("#GTEx-version").multiselect({
        buttonWidth: "400px",
        checkboxName: function (option) {
          return "multiselect[]";
        },
      });

      const inputs = localStorage.getItem("formVals");

      if (inputs) {
        const formVals = JSON.parse(inputs);

        Object.entries(formVals).forEach(([k, v]) => {
          if (formFields[k].type === "input") {
            d3.select(formFields[k].selector).node().value = v;
          } else if (formFields[k].type === "select") {
            d3.select(formFields[k].selector).node().value = v;
          } else if (formFields[k].type === "multiselect") {
            if (v) {
              $(formFields[k].selector).val(v);
              $(formFields[k].selector).multiselect("refresh");
            }
          }
        });
      }
    });
  });

  // initialize popover and tooltip
  $(function () {
    $('[data-toggle="popover"]').popover();
  });
  $(document).ready(function () {
    $('[data-toggle="tooltip"]').tooltip({
      delay: { show: 500, hide: 100 },
    });
  });

  // add warning listeners
  ["#locus", "#coordinate", "#file-upload"].forEach((id) =>
    $(id).on("change", updateChrXWarning)
  );
}

init();
