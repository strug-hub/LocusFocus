var genesfile = $("head").data("genesfile");
var sessionfile = $("head").data("session_file");
var SSPvalues_file = $("head").data("sspvalues_file");
var coloc2_file = $("head").data("coloc2_file");
var sessionid = $("head").data("sessionid");
var metadata_file = $("head").data("metadata_file");
var staticEndpoint = $("head").data("static_endpoint");
var transpose = false;

// initialize popover, tooltip, callbacks
$(document).ready(function () {
  $('[data-toggle="tooltip"]').tooltip({
    delay: { show: 500, hide: 100 },
  });
  $('[data-toggle="popover"]').popover();
  $("#selGene").on("change", (e) => optionChanged(e.currentTarget.value));
  $("#pval-filter").on("change", (e) => plot_fullgwas(e.currentTarget));
});

if (staticEndpoint.endsWith("/")) {
  staticEndpoint = staticEndpoint.slice(0, -1);
}

function getResourceUrl(resource) {
  return `${staticEndpoint}/${resource}`;
}

function showSection(name, hidden = false) {
  let rowId = "#" + name + "-row";
  let dividerId = "#" + name + "-divider";
  $(rowId).prop("hidden", hidden);
  $(dividerId).prop("hidden", hidden);
}

function init() {
  d3.select("#plot_message")
    .append("i")
    .classed("fa", true)
    .classed("fa-info-circle", true)
    .text(`Please wait while loading colocalization plot`);
  d3.select("#sessionid")
    .append("a")
    .attr("href", `/session_id/${sessionid}`)
    .text(sessionid);

  let filePromises = [
    d3.json(getResourceUrl(sessionfile)),
    d3.json(getResourceUrl(metadata_file)),
  ];

  Promise.allSettled(filePromises).then((results) => {
    let [sessionData, metadataResponseJson] = results.map((v) => {
      return v.status === "fulfilled" ? v.value : null;
    });

    $("#loading-spinner").remove();
    if (
      metadataResponseJson !== null &&
      metadataResponseJson["type"] === "set-based-test"
    ) {
      // set based test case
      // show elements used for set-based test results
      prepareSBTResults(sessionData, sessionid);
    } else {
      // fetch other files for default tests
      Promise.allSettled([
        d3.json(getResourceUrl(genesfile)),
        d3.json(getResourceUrl(SSPvalues_file)),
        d3.json(getResourceUrl(coloc2_file)),
      ]).then((SSfiles) => {
        let [genesData, SSResponseJson, coloc2ResponseJson] = SSfiles.map((v) =>
          v.status === "fulfilled" ? v.value : null
        );
        prepareResults(
          sessionData,
          sessionid,
          genesData,
          SSResponseJson,
          coloc2ResponseJson
        );
      });
    }
  });
}

function populate_genes(genesdata, response) {
  var selector = d3.select("#selGene");
  var genenames = genesdata.map((gene) => gene["name"]);
  var geneindex = 0;
  Object.keys(genenames).forEach((gene, i) => {
    if (genenames[i] === response["gene"]) {
      geneindex = i;
    }
    selector
      .append("option")
      .text(genenames[i])
      .property("value", genenames[i]);
  });
  selector.property("selectedIndex", geneindex);
}

function optionChanged(newgene) {
  var newurl = `/update/${sessionid}/${newgene}`;
  d3.select("#plot").text("");
  // d3.select('#plot_message').append("p").attr("class","text-warning").text(`Please wait while re-drawing eQTLs for gene ${newgene}`);
  d3.select("#plot_message")
    .append("i")
    .classed("fa", true)
    .classed("fa-info-circle", true)
    .text(`Please wait while re-drawing eQTLs for gene ${newgene}`);
  // d3.json("{{ url_for('update_colocalizing_gene', session_id=sessionid, newgene=newgene) }}").then(response => {
  d3.json(newurl).then((newresponse) => {
    d3.json(getResourceUrl(genesfile)).then((genesResponse) => {
      var newdata = newresponse;
      var genesdata = genesResponse;
      plot_gwas(newdata, genesdata);
      d3.select("#plot_message").text("");
      // plot_heatmap(SSResponse.Genes, SSResponse.Tissues, SSResponse.SSPvalues);
      // buildTable(SSResponse.Genes, SSResponse.Tissues, SSResponse.SSPvalues);
    });
  });
}

function transposeTable() {
  transpose = !transpose;
  d3.json(getResourceUrl(SSPvalues_file)).then((SSResponse) => {
    buildTable(
      SSResponse.Genes,
      SSResponse.Tissues,
      SSResponse.SSPvalues,
      transpose
    );
    buildNTable(
      SSResponse.Genes,
      SSResponse.Tissues,
      SSResponse.Num_SNPs_Used_for_SS,
      transpose
    );
  });
}

function prepareResults(
  sessionData,
  sessionid,
  genesData,
  SSResponseJson,
  coloc2ResponseJson
) {
  // All code for displaying results for simple sum, coloc2, etc.
  [
    "params-table",
    "plot",
    "heatmap",
    "variants-table",
    "numSSsnpsUsed-table",
    "secondary-table",
    "SSguidance-table",
    "coloc2-table",
    "set-based-test-table",
  ].forEach((name) => showSection(name));

  if (sessionData !== null) {
    if (sessionData["snp_warning"] === true) {
      d3.select("#snpWarning")
        .append("div")
        .classed("alert alert-danger", true)
        .append("p")
        .append("center")
        .append("strong")
        .html(
          "Warning: The number of SNPs that match GTEx variants is lower than 80%.<br> Please ensure variant names given are accurate, including REF and ALT columns."
        );
    }
    buildParamsTable(sessionData, sessionid, (type = "default"));
    if (genesData !== null) {
      populate_genes(genesData, sessionData);
      plot_gwas(sessionData, genesData);
    }
    if (sessionData["liftover_warning"]?.length > 0) {
      d3.select("#liftover-warning")
        .property("hidden", false)
        .append("strong")
        .text(sessionData["liftover_warning"]);
    }
  }
  if (SSResponseJson !== null) {
    if (SSResponseJson["Tissues"].length === 0) {
      d3.selectAll(".gtex_heatmap_message")
        .append("p")
        .attr("class", "text-warning")
        .text(`No GTEx tissues were selected`);
    }
    if (SSResponseJson["Tissues"].length > 0) {
      plot_heatmap(
        SSResponseJson.Genes,
        SSResponseJson.Tissues,
        SSResponseJson.SSPvalues,
        SSResponseJson.SSPvalues_secondary
      );
      buildTable(
        SSResponseJson.Genes,
        SSResponseJson.Tissues,
        SSResponseJson.SSPvalues
      );
      buildNTable(
        SSResponseJson.Genes,
        SSResponseJson.Tissues,
        SSResponseJson.Num_SNPs_Used_for_SS
      );
    }
    if (SSResponseJson["Secondary_dataset_titles"].length === 0) {
      d3.select(".secondary_datasets_message")
        .append("p")
        .attr("class", "text-warning")
        .text(`No secondary datasets uploaded`);
    }
    if (SSResponseJson["Secondary_dataset_titles"].length > 0) {
      list_secondary_SSPvalues(
        SSResponseJson.Secondary_dataset_titles,
        SSResponseJson.SSPvalues_secondary,
        SSResponseJson.Num_SNPs_Used_for_SS_secondary
      );
    }
    buildSSguidanceTable(
      SSResponseJson.Genes,
      SSResponseJson.Tissues,
      SSResponseJson.SSPvalues,
      SSResponseJson.SSPvalues_secondary
    );
  }
  if (coloc2ResponseJson !== null) {
    if (coloc2ResponseJson["ProbeID"].length === 0) {
      d3.select("#coloc2-message")
        .append("p")
        .attr("class", "text-warning")
        .text(`Insufficient data for COLOC2 calculations`);
    }
    if (coloc2ResponseJson["ProbeID"].length > 0) {
      buildColoc2Table(coloc2ResponseJson);
    }
  }
  d3.select("#plot_message").text("");
  buildFirstStageTable(sessionData);
}

function prepareSBTResults(sessionData, sessionid) {
  // All code for displaying set-based test results
  ["params-table", "set-based-test-table"].forEach((name) => showSection(name));
  buildParamsTable(sessionData, sessionid, (type = "set-based-test"));
  buildFirstStageTable(sessionData);
}

init();

// If unchecked p-value filter, redraw all points:
function plot_fullgwas(pval_filter_box) {
  d3.select("#plot").text("");
  d3.select("#plot_message")
    .append("i")
    .classed("fa", true)
    .classed("fa-info-circle", true)
    .text(`Please wait while re-drawing`);
  d3.json(getResourceUrl(sessionfile)).then((response) => {
    d3.json(getResourceUrl(genesfile)).then((genesResponse) => {
      var data = response;
      var genesdata = genesResponse;
      if (pval_filter_box.checked === false) {
        // d3.select('#plot_message').append("p").attr("class","text-warning").text(`Please wait while re-drawing with all GWAS points (slow)`);
        d3.select("#plot_message")
          .append("i")
          .classed("fa", true)
          .classed("fa-info-circle", true)
          .text(`Please wait while re-drawing with all GWAS points (slow)`);
        plot_gwas(
          data,
          genesdata,
          (eqtl_smoothing_window_size = -1),
          (percent_occupied_by_one_char = 0.02),
          (inputHeight = 720),
          (inputWidth = 1080),
          (font_size = 14),
          (legend_offset = 0.1),
          (pval_filter = false)
        );
      } else {
        plot_gwas(data, genesdata);
      }
      d3.select("#plot_message").text("");
    });
  });
}

// Listen for redraw click events:
var coloc_plot_redraw_btn = d3.select("#colocPlot-redraw-btn");
coloc_plot_redraw_btn.on("click", function () {
  d3.event.preventDefault();

  // Get input values:
  var coloc_plot_width = d3.select("#colocPlotWidth").property("value");
  var coloc_plot_height = d3.select("#colocPlotHeight").property("value");
  var eqtl_window_size = d3.select("#eqtlWindow").property("value");
  var percent_occupied_by_one_char = d3.select("#pctOneChar").property("value");
  var coloc_plot_fontsize = d3.select("#colocPlotFontSize").property("value");
  var legendOffset = d3.select("#legendOffset").property("value");
  var selectedgene = d3.select("#selGene").property("value");
  var url = `/update/${sessionid}/${selectedgene}`;

  if (!coloc_plot_width) {
    coloc_plot_width = 1080;
  }
  if (!coloc_plot_height) {
    coloc_plot_height = 720;
  }
  if (!eqtl_window_size) {
    eqtl_window_size = -1;
  }
  if (!percent_occupied_by_one_char) {
    percent_occupied_by_one_char = 0.1;
  }
  if (!coloc_plot_fontsize) {
    coloc_plot_fontsize = 14;
  }
  if (!legendOffset) {
    legendOffset = 0.1;
  }

  // Clear current colocalization plot
  d3.select("#plot").text("");

  // Redraw colocalization plot with updated parameters
  d3.json(url).then((response) => {
    d3.json(getResourceUrl(genesfile)).then((genesResponse) => {
      var data = response;
      var genesdata = genesResponse;
      plot_gwas(
        data,
        genesdata,
        (eqtl_smoothing_window_size = +eqtl_window_size),
        (percent_occupied_by_one_char = +percent_occupied_by_one_char),
        (inputHeight = +coloc_plot_height),
        (inputWidth = +coloc_plot_width),
        (font_size = +coloc_plot_fontsize),
        (legend_offset = +legendOffset)
      );
    });
  });
});

// Listen for svg plot requests for colocalization plot:
var coloc_plot_svg_btn = d3.select("#colocPlot-svg-btn");
coloc_plot_svg_btn.on("click", function () {
  d3.event.preventDefault();

  colocPlot_container = d3
    .select("#plot")
    .select("div.plot-container")
    .select("div.svg-container");

  // Get input values:
  var coloc_plot_width = d3.select("#colocPlotWidth").property("value");
  var coloc_plot_height = d3.select("#colocPlotHeight").property("value");

  if (!coloc_plot_width) {
    coloc_plot_width = 1080;
  }
  if (!coloc_plot_height) {
    coloc_plot_height = 720;
  }

  saveSvg(
    colocPlot_container,
    "colocalization-plot.svg",
    coloc_plot_width,
    coloc_plot_height
  );
});

function saveSvg(svgDiv, name, w = 1080, h = 720) {
  var svgEls = svgDiv.selectAll(".main-svg");
  var head = `<svg class="main-svg" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="${w}" height="${h}" style="background: rgb(255, 255, 255);">`;

  var svgData = head;
  svgEls.each(function (d, i) {
    svgEl = d3.select(this);
    svgData += svgEl.html();
  });
  svgData += "</svg>";
  svgData = svgData.replace('"<', '"&lt;');
  svgData = svgData.replace("<br>", "\\r\\n");
  var preface = '<?xml version="1.0" standalone="no"?>\r\n';
  var svgBlob = new Blob([preface, svgData], {
    type: "image/svg+xml;charset=utf-8",
  });
  saveAs(svgBlob, name);
}

// Listen for redraw click events:
var heatmap_redraw_btn = d3.select("#heatmap-redraw-btn");
heatmap_redraw_btn.on("click", function () {
  d3.event.preventDefault();

  // Get width and height
  var widthinput = +d3.select("#heatmapWidth").property("value");
  var heightinput = +d3.select("#heatmapHeight").property("value");
  var heatmap_fontsize = +d3.select("#heatmapFontsize").property("value");

  if (!widthinput) {
    widthinput = 1080;
  }
  if (!heightinput) {
    heightinput = 720;
  }
  if (!heatmap_fontsize) {
    heatmap_fontsize = 12;
  }

  // Clear current heatmap
  d3.select("#heatmap").text("");

  // Redraw heatmap with updated width and height
  d3.json(getResourceUrl(SSPvalues_file)).then((SSResponse) => {
    plot_heatmap(
      SSResponse.Genes,
      SSResponse.Tissues,
      SSResponse.SSPvalues,
      SSResponse.SSPvalues_secondary,
      widthinput,
      heightinput,
      heatmap_fontsize
    );
  });
});

// Listen for svg plot requests for heatmap plot:
var heatmap_svg_btn = d3.select("#heatmap-svg-btn");
heatmap_svg_btn.on("click", function () {
  d3.event.preventDefault();

  heatmap_container = d3
    .select("#heatmap")
    .select("div.plot-container")
    .select("div.svg-container");

  // Get input values:
  var widthinput = +d3.select("#heatmapWidth").property("value");
  var heightinput = +d3.select("#heatmapHeight").property("value");

  if (!widthinput) {
    widthinput = 1080;
  }
  if (!heightinput) {
    heightinput = 720;
  }

  saveSvg(heatmap_container, "SS-heatmap.svg", widthinput, heightinput);
});
