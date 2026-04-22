function smoothing(x, y, xrange, windowCount) {
  if (windowCount <= 0) {
    // default window count
    windowCount = ((xrange[1] - xrange[0]) / 1000000) * 150;
  }

  if (!x.length) return [[], []];

  const windowWidth = (xrange[1] - xrange[0]) / windowCount;

  const smoothX = [];
  const smoothY = [];

  let left = 0;
  let right = 0;
  let curr = xrange[0];

  const n = x.length;

  while (curr < xrange[1]) {
    const windowEnd = curr + windowWidth;

    // Advance right pointer to include points in window
    while (right < n && x[right] <= windowEnd) {
      right++;
    }

    // Advance left pointer to exclude points before window
    while (left < n && x[left] < curr) {
      left++;
    }

    // Find max y in [left, right)
    let maxY = -Infinity;
    let maxIdx = -1;

    for (let i = left; i < right; i++) {
      if (y[i] > maxY) {
        maxY = y[i];
        maxIdx = i;
      }
    }

    if (maxIdx !== -1) {
      smoothX.push(x[maxIdx]);
      smoothY.push(Math.max(0, maxY));
    }

    curr += windowWidth + 1;
  }

  return [smoothX, smoothY];
}


function withDefaults(opts) {
  return {
    eqtlSmoothingWindow: -1,
    pvalFilter: true,
    markerSize: 10,
    leadMarkerScale: 1.5,
    fontSize: 14,
    legendOffset: 0.1,
    width: 1080,
    height: 720,
    percentCharWidth: 0.02,
    ...opts
  };
}


function prepareGWASData(data, cfg) {
  let {
    positions,
    pvalues,
    ld_values,
    snps,
    lead_snp,
    startbp,
    endbp,
    chrom,
    coordinate
  } = data;

  if (cfg.pvalFilter) {
    const kept = pvalues
      .map((p, i) => (p < 0.1 ? i : -1))
      .filter(i => i !== -1);

    positions = kept.map(i => positions[i]);
    pvalues = kept.map(i => pvalues[i]);
    ld_values = kept.map(i => ld_values[i]);
    snps = kept.map(i => snps[i]);
  }

  const log10p = pvalues.map(p => -Math.log10(p));
  const regionSize = endbp - startbp;

  return {
    positions,
    pvalues,
    log10p,
    ld_values,
    snps,
    lead_snp,
    leadIndex: snps.indexOf(lead_snp),
    startbp,
    endbp,
    chrom: chrom === 23 ? "X" : chrom,
    build: coordinate,
    regionSize,
    yMax: d3.max(log10p)
  };
}


const LD_BINS = [
  { name: "No LD Info", test: v => v < 0, color: "#7f7f7f" },
  { name: "< 0.2", test: v => Math.abs(v) < 0.2, color: "#1f77b4" },
  { name: "0.2–0.4", test: v => Math.abs(v) >= 0.2 && Math.abs(v) < 0.4, color: "#17becf" },
  { name: "0.4–0.6", test: v => Math.abs(v) >= 0.4 && Math.abs(v) < 0.6, color: "#bcbd22" },
  { name: "0.6–0.8", test: v => Math.abs(v) >= 0.6 && Math.abs(v) < 0.8, color: "#ff7f0e" },
  { name: "≥ 0.8", test: v => Math.abs(v) >= 0.8, color: "#d62728" }
];


function makeLDTraces(gwas, cfg) {
  const traces = [];

  const leadIndex = gwas.leadIndex;
  const used = new Set();

  // Never allow lead SNP to be picked up by LD bins
  if (leadIndex !== -1) {
    used.add(leadIndex);
  }

  const LD_BINS = [
    {
      name: "< 0.2",
      test: v => v >= 0 && Math.abs(v) < 0.2,
      color: "#1f77b4"
    },
    {
      name: "0.2–0.4",
      test: v => Math.abs(v) >= 0.2 && Math.abs(v) < 0.4,
      color: "#17becf"
    },
    {
      name: "0.4–0.6",
      test: v => Math.abs(v) >= 0.4 && Math.abs(v) < 0.6,
      color: "#bcbd22"
    },
    {
      name: "0.6–0.8",
      test: v => Math.abs(v) >= 0.6 && Math.abs(v) < 0.8,
      color: "#ff7f0e"
    },
    {
      name: "≥ 0.8",
      test: v => Math.abs(v) >= 0.8,
      color: "#d62728"
    }
  ];

  // Build LD bin traces
  LD_BINS.forEach(bin => {
    const indices = [];

    for (let i = 0; i < gwas.ld_values.length; i++) {
      const ld = gwas.ld_values[i];

      // Critical guard:
      // Do NOT plot any negative LD values (error sentinels)
      if (ld < 0) continue;
      if (used.has(i)) continue;

      if (bin.test(ld)) {
        indices.push(i);
        used.add(i);
      }
    }

    if (indices.length > 0) {
      traces.push({
        x: indices.map(i => gwas.positions[i]),
        y: indices.map(i => gwas.log10p[i]),
        text: indices.map(i => gwas.snps[i]),
        name: bin.name,
        type: "scatter",
        mode: "markers",
        marker: {
          size: cfg.markerSize,
          color: bin.color
        }
      });
    }
  });

  // Lead SNP trace (always plotted)
  if (leadIndex !== -1) {
    traces.push({
      x: [gwas.positions[leadIndex]],
      y: [gwas.log10p[leadIndex]],
      text: [gwas.lead_snp],
      name: "Lead SNP",
      type: "scatter",
      mode: "markers",
      marker: {
        size: cfg.markerSize * cfg.leadMarkerScale,
        color: "#9467bd"
      }
    });
  }

  return traces;
}


function prepareEQTLData(data, names, gwas, cfg, cols = {}) {
  const result = {};

  names.forEach(name => {
    const x = [];
    const y = [];
    const snps = [];

    data[name].forEach(row => {
      x.push(+row[cols.pos || "variant_pos" || "pos"] ?? +row.pos);
      y.push(-Math.log10(+row[cols.pval || "pval"]));
      if (cols.snp || row.rs_id) snps.push(row[cols.snp || "rs_id"]);
    });

    result[name] = {
      raw: { x, y, snps },
      smooth: smoothing(x, y, [gwas.startbp, gwas.endbp], cfg.eqtlSmoothingWindow)
    };
  });

  return result;
}


function prepareSecondaryData(data, gwas, cfg) {
  const titles = data.secondary_dataset_titles;
  const colnames = data.secondary_dataset_colnames;

  const positionCol = colnames[1];
  const snpCol = colnames[2];
  const pvalCol = colnames[3];

  const result = {};

  titles.forEach(title => {
    const x = [];
    const y = [];
    const snps = [];

    const records = data[title] || [];

    records.forEach(row => {
      const pos = +row[positionCol];
      const pval = +row[pvalCol];

      if (Number.isFinite(pos) && Number.isFinite(pval) && pval > 0) {
        x.push(pos);
        y.push(-Math.log10(pval));
        snps.push(row[snpCol]);
      }
    });

    const smooth =
      cfg.eqtlSmoothingWindow > 0
        ? smoothing(
            x,
            y,
            [gwas.startbp, gwas.endbp],
            cfg.eqtlSmoothingWindow
          )
        : [[], []];

    result[title] = {
      raw: {
        x,
        y,
        snps
      },
      smooth
    };
  });

  return result;
}


function makeEQTLTraces(datasets) {
  return Object.entries(datasets).map(([name, d]) => ({
    x: d.smooth[0],
    y: d.smooth[1],
    name,
    mode: "lines",
    yaxis: "y2"
  }));
}


function makeEQTLMarkerTraces(datasets) {
  return Object.entries(datasets).map(([name, d]) => ({
    x: d.raw.x,
    y: d.raw.y,
    text: d.raw.snps,
    name,
    mode: "markers",
    visible: "legendonly",
    marker: { size: 8, opacity: 0.3 },
    yaxis: "y2"
  }));
}


function layoutGenes(genesData, gwas, cfg) {
  // Defensive copy + sort by start
  const genes = [...genesData].sort((a, b) => a.txStart - b.txStart);

  /* Row assignment (no overlap) */

  function overlaps(a, b) {
    return (b[0] <= a[1] && b[1] >= a[0]);
  }

  const rows = [];
  genes.forEach((gene, i) => {
    const bin = [gene.txStart, gene.txEnd];
    let rowIndex = 0;

    while (true) {
      if (!rows[rowIndex]) {
        rows[rowIndex] = [bin];
        gene.geneRow = rowIndex + 1;
        break;
      }
      if (!rows[rowIndex].some(r => overlaps(r, bin))) {
        rows[rowIndex].push(bin);
        gene.geneRow = rowIndex + 1;
        break;
      }
      rowIndex++;
    }
  });

  /* Vertical geometry */

  const logRange = gwas.yMax - d3.min(gwas.log10p);
  const geneAreaPct = 0.15;
  const maxRows = 7;

  const geneAreaHeight = Math.min(
    geneAreaPct * logRange * rows.length,
    geneAreaPct * logRange * maxRows
  );

  const rowHeight = geneAreaHeight / rows.length;
  const textHeight = rowHeight * 0.15;
  const geneMargin = rowHeight * 0.15;
  const exonHeight = rowHeight - 2 * (geneMargin + textHeight);
  const intronHeight = exonHeight * 0.4;

  /* Shapes (exons & introns) */

  const shapes = [];

  genes.forEach(gene => {
    const baseY =
      -(gene.geneRow * rowHeight) + textHeight + geneMargin;

    // Intron
    shapes.push({
      type: "rect",
      xref: "x",
      yref: "y",
      x0: gene.txStart,
      x1: gene.txEnd,
      y0: baseY + (exonHeight - intronHeight) / 2,
      y1: baseY + (exonHeight - intronHeight) / 2 + intronHeight,
      line: { color: "rgb(55,128,191)", width: 1 },
      fillcolor: "rgb(55,128,191)"
    });

    // Exons
    for (let i = 0; i < gene.exonStarts.length; i++) {
      shapes.push({
        type: "rect",
        xref: "x",
        yref: "y",
        x0: gene.exonStarts[i],
        x1: gene.exonEnds[i],
        y0: baseY,
        y1: baseY + exonHeight,
        line: { color: "rgb(55,128,191)", width: 1 },
        fillcolor: "rgb(55,128,191)"
      });
    }
  });

  /* Gene name placement (collision‑aware) */

  const regionSize = gwas.regionSize;
  const extraX = 0.05 * regionSize;
  const fullY = gwas.yMax + geneAreaHeight;

  function rectOverlap(a, b) {
    return (
      a[0] <= b[2] &&
      a[2] >= b[0] &&
      a[1] <= b[3] &&
      a[3] >= b[1]
    );
  }

  const placedRects = [];
  const finalX = [];
  const finalY = [];
  const finalText = [];

  genes.forEach(gene => {
    const midX = (gene.txStart + gene.txEnd) / 2;
    const baseY =
      -(gene.geneRow * rowHeight) + textHeight + geneMargin;

    const nameWidth =
      (gene.name.length / 2) * cfg.percentCharWidth;

    const xNorm =
      (midX - (gwas.startbp - extraX)) /
      (regionSize + 2 * extraX);

    const x0 = Math.max(xNorm - nameWidth, 0);
    const x1 = xNorm + nameWidth;

    // Bottom position
    let y0 = (geneAreaHeight - -baseY - cfg.fontSize * 0.5) / fullY;
    let y1 = (geneAreaHeight - -baseY) / fullY;

    let rect = [x0, y0, x1, y1];
    let yText = baseY;

    // Try top if overlapping
    if (placedRects.some(r => rectOverlap(rect, r))) {
      y0 =
        (geneAreaHeight - -baseY - cfg.fontSize * 0.5 + exonHeight) /
        fullY;
      y1 =
        (geneAreaHeight - -baseY + exonHeight) / fullY;
      rect = [x0, y0, x1, y1];
      yText = baseY + exonHeight;

      if (placedRects.some(r => rectOverlap(rect, r))) {
        return; // hide label
      }
    }

    placedRects.push(rect);
    finalX.push(midX);
    finalY.push(yText);
    finalText.push(gene.name);
  });

  const nameTrace = {
    x: finalX,
    y: finalY,
    text: finalText,
    mode: "markers+text",
    type: "scatter",
    marker: { opacity: 0 },
    textposition: "bottom",
    font: { style: "italic" },
    showlegend: false,
    yaxis: "y"
  };

  return {
    shapes,
    nameTrace,
    height: geneAreaHeight,
    rows: rows.length
  };
}


function getSecondaryYMax(...datasets) {
  let ymax = 1;

  datasets.forEach(group => {
    Object.values(group).forEach(dataset => {
      const y = dataset.smooth?.[1];
      if (Array.isArray(y) && y.length > 0) {
        const m = d3.max(y);
        if (m > ymax) ymax = m;
      }
    });
  });

  return ymax;
}


function makeLayout(gwas, genes, gtex, secondary, cfg, gtexYMax) {
  const extraX = 0.05 * gwas.regionSize;
  const extraY = 0.05 * gwas.yMax;

  // Primary axis (GWAS)
  const yMin1 = -genes.height;
  const yMax1 = gwas.yMax + extraY;

  // Fractional position of zero on primary axis
  const zeroFrac = (0 - yMin1) / (yMax1 - yMin1);

  // Secondary axis (GTEx / secondary datasets)
  const yMax2 = gtexYMax + extraY * (gtexYMax / gwas.yMax);
  const yMin2 = -zeroFrac * yMax2 / (1 - zeroFrac);

  return {
    xaxis: {
      range: [gwas.startbp - extraX, gwas.endbp + extraX],
      title: `Chromosome ${gwas.chrom} (${gwas.build})`
    },
    yaxis: {
      range: [yMin1, yMax1],
      title: "-log10(p)"
    },
    yaxis2: {
      range: [yMin2, yMax2],
      overlaying: "y",
      side: "right",
      title: "Secondary datasets -log10(p)",
      showgrid: false
    },
    shapes: genes.shapes,
    legend: {
      x: 1 + cfg.legendOffset,
      y: 1,
      font: { size: cfg.fontSize }
    },
    width: cfg.width,
    height: cfg.height,
    hovermode: "closest",
    font: { size: cfg.fontSize }
  };
}


function makeGeneNameTraces(genes) {
  if (!genes || !genes.nameTrace) {
    return [];
  }

  return [genes.nameTrace];
}


function plotGWAS(
  data,
  genesData,
  options = {}
) {
  const cfg = withDefaults(options);

  const gwas = prepareGWASData(data, cfg);
  const genes = layoutGenes(genesData, gwas, cfg);

  const gtex = prepareEQTLData(data, data.gtex_tissues, gwas, cfg);
  const secondary = prepareSecondaryData(data, gwas, cfg);

  const traces = [
    ...makeLDTraces(gwas, cfg),
    ...makeEQTLTraces(gtex),
    ...makeEQTLMarkerTraces(gtex),
    ...makeEQTLTraces(secondary),
    ...makeEQTLMarkerTraces(secondary),
    ...makeGeneNameTraces(genes)
  ];

  const gtexYMax = getSecondaryYMax(gtex, secondary);
  const layout = makeLayout(gwas, genes, gtex, secondary, cfg, gtexYMax);

  Plotly.newPlot("plot", traces, layout);
}

function plot_gwas(data, genesdata,
  eqtl_smoothing_window_size = -1,
  percent_occupied_by_one_char = 0.020,
  inputHeight = 720,
  inputWidth = 1080,
  font_size = 14,
  legend_offset = 0.1,
  pval_filter = true) {

  var positions = data.positions;
  var pvalues = data.pvalues;
  var ld_values = data.ld_values;
  var chrom = data.chrom;
  if (chrom === 23) chrom = "X";
  var startbp = data.startbp;
  var endbp = data.endbp;
  var snps = data.snps;
  var lead_snp = data.lead_snp;
  var gtex_tissues = data.gtex_tissues;
  var SS_start = +data.SS_region[0];
  var SS_end = +data.SS_region[1];
  var build = data.coordinate;

  var secondary_dataset_titles = data.secondary_dataset_titles;
  var secondary_dataset_position_colname = data.secondary_dataset_colnames[1];
  var secondary_dataset_snp_colname = data.secondary_dataset_colnames[2];
  var secondary_dataset_pval_colname = data.secondary_dataset_colnames[3];

  if (pval_filter) {
    pindices = [];
    pvalues.map((p, i) => {
      if (p < 0.1) pindices.push(i);
    });
    pAtIndices = pindices.map(i => pvalues[i]);
    positionsAtIndices = pindices.map(i => positions[i]);
    ld_valuesAtIndices = pindices.map(i => ld_values[i]);
    snpsAtIndices = pindices.map(i => snps[i]);
    pvalues = pAtIndices;
    positions = positionsAtIndices;
    ld_values = ld_valuesAtIndices;
    snps = snpsAtIndices;
  }
  var log10pvalues = pvalues.map(p => -Math.log10(p));
  var no_ld_info_snps = [], no_ld_info_snps_color = "#7f7f7f"; // grey
  var ld_lt_20_group = [], ld_lt_20_group_color = "#1f77b4"; // very low ld points (dark blue)
  var ld_20_group = [], ld_20_group_color = "#17becf"; // low ld points (light blue)
  var ld_40_group = [], ld_40_group_color = "#bcbd22"; // green
  var ld_60_group = [], ld_60_group_color = "#ff7f0e"; // orange
  var ld_80_group = [], ld_80_group_color = "#d62728"; // red
  var lead_snp_index = 0, lead_snp_color = "#9467bd"; // purple
  var markersize = 10;
  var lead_markersize = 10 * 1.5;
  var ld_colors = [];
  var regionsize = endbp - startbp;
  var log10pvalue_range = d3.max(pvalues.map(p => -Math.log10(p))) - d3.min(pvalues.map(p => -Math.log10(p)));
  var extra_x_range = 0.05 * regionsize;
  var extra_y_range = 0.05 * log10pvalue_range;
  var eqtl_window_multiplier = 150;
  if (eqtl_smoothing_window_size === -1) {
    eqtl_smoothing_window_size = (regionsize / 1000000) * eqtl_window_multiplier;
  }
  var font_height = 0.5;

  // Helper functions:

  function getYmax(gtex_traces, secondary_traces) {
    ymax = 1;
    for (var i = 0; i < gtex_tissues.length; i++) {
      currymax = d3.max(gtex_traces[gtex_tissues[i]][1]);
      if (currymax > ymax) ymax = currymax;
    }
    for (var i = 0; i < secondary_dataset_titles.length; i++) {
      currymax = d3.max(secondary_traces[secondary_dataset_titles[i]][1]);
      if (currymax > ymax) ymax = currymax;
    }
    return ymax;
  }

  // Add the row number each gene will be plotted into:
  genesdata.sort(function (a, b) {
    return a.txStart - b.txStart;
  });

  function overlap(bin1, bin2) {
    return (bin2[0] >= bin1[0] && bin2[0] <= bin1[1]) || (bin2[1] >= bin1[0] && bin2[1] <= bin1[1]);
  }

  function checkSpace(geneRows, rowNum, bin) {
    occupied = false;
    for (var i = 0; i < geneRows[rowNum].length; i++) {
      if (overlap(geneRows[rowNum][i], bin)) {
        occupied = true;
        return occupied;
      }
    }
    return occupied;
  }

  var firstBin = [genesdata[0].txStart, genesdata[0].txEnd];
  var geneRows = [];
  geneRows[0] = [firstBin];
  genesdata[0]['geneRow'] = 1;

  for (var i = 1; i < genesdata.length; i++) {
    currRow = 0;
    geneBin = [genesdata[i].txStart, genesdata[i].txEnd];
    occupied = checkSpace(geneRows, currRow, geneBin);
    while (occupied && currRow < geneRows.length) {
      currRow++;
      if (currRow < geneRows.length) occupied = checkSpace(geneRows, currRow, geneBin);
    }
    if (currRow === geneRows.length) {
      geneRows[currRow] = [geneBin];
      genesdata[i]['geneRow'] = currRow + 1;
    }
    else {
      geneRows[currRow].push(geneBin);
      genesdata[i]['geneRow'] = currRow + 1;
    }
  }

  var gene_area_percentage = 0.15;
  var max_num_gene_rows = 7;
  var gene_area_height = d3.min([gene_area_percentage * log10pvalue_range * geneRows.length, gene_area_percentage * max_num_gene_rows * log10pvalue_range]);
  var row_height = gene_area_height / geneRows.length;
  var text_height = row_height * 0.15
  var gene_margin = row_height * 0.15;
  var exon_height = row_height - (2 * (gene_margin + text_height));
  var intron_height = exon_height * 0.4;
  var rectangle_shapes = [];
  var annotations_x = [];
  var annotations_y = [];
  var annotations_text = [];
  for (var i = 0; i < genesdata.length; i++) {
    // build intron rectangle shapes for each gene:
    var rectangle_shape = {
      type: 'rect',
      xref: 'x',
      yref: 'y',
      x0: genesdata[i]['txStart'],
      y0: -(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin + ((exon_height - intron_height) / 2),
      x1: genesdata[i]['txEnd'],
      y1: -(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin + ((exon_height - intron_height) / 2) + intron_height,
      line: {
        color: 'rgb(55, 128, 191)',
        width: 1
      },
      fillcolor: 'rgba(55, 128, 191, 1)'
    }
    rectangle_shapes.push(rectangle_shape);
    annotations_x.push(genesdata[i]['txStart']);
    annotations_x.push((genesdata[i]['txStart'] + genesdata[i]['txEnd']) / 2);
    annotations_x.push(genesdata[i]['txEnd']);
    var y = -(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin + ((exon_height - intron_height) / 2) + intron_height / 2;
    annotations_y.push(y);
    annotations_y.push(y);
    annotations_y.push(y);
    annotations_text.push(genesdata[i]['name']);
    annotations_text.push(genesdata[i]['name']);
    annotations_text.push(genesdata[i]['name']);
    for (var j = 0; j < genesdata[i]['exonStarts'].length; j++) {
      // build exon rectangle shapes for current gene
      var rectangle_shape = {
        type: 'rect',
        xref: 'x',
        yref: 'y',
        x0: genesdata[i]['exonStarts'][j],
        y0: -(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin,
        x1: genesdata[i]['exonEnds'][j],
        y1: -(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin + exon_height,
        line: {
          color: 'rgb(55, 128, 191)',
          width: 1
        },
        fillcolor: 'rgba(55, 128, 191, 1)'
      }
      rectangle_shapes.push(rectangle_shape);
    }
  }

  // Smooth out each GTEx tissue's association results for plotting as lines:
  gtex_line_traces = {};
  gtex_positions = {};
  gtex_log10_pvalues = {};
  gtex_snps = {};
  for (var i = 0; i < gtex_tissues.length; i++) {
    gtex_positions[gtex_tissues[i]] = [];
    gtex_log10_pvalues[gtex_tissues[i]] = [];
    gtex_snps[gtex_tissues[i]] = [];
    data[gtex_tissues[i]].forEach(eqtl => {
      Object.keys(eqtl).forEach(k => {
        if (k === 'variant_pos' || k === 'pos') {
          gtex_positions[gtex_tissues[i]].push(+eqtl[k]);
        }
        else if (k === 'pval') {
          gtex_log10_pvalues[gtex_tissues[i]].push(-Math.log10(+eqtl[k]));
        }
        else if (k === 'rs_id') {
          gtex_snps[gtex_tissues[i]].push(eqtl[k]);
        }
      });
    });
    gtex_line_traces[gtex_tissues[i]] = smoothing(gtex_positions[gtex_tissues[i]], gtex_log10_pvalues[gtex_tissues[i]],
      [startbp, endbp], eqtl_smoothing_window_size);
    // console.log(data['Pancreas']);
  }

  secondary_line_traces = {};
  secondary_positions = {};
  secondary_log10_pvalues = {};
  secondary_snps = {};
  for (var i = 0; i < secondary_dataset_titles.length; i++) {
    secondary_positions[secondary_dataset_titles[i]] = [];
    secondary_log10_pvalues[secondary_dataset_titles[i]] = [];
    secondary_snps[secondary_dataset_titles[i]] = [];
    data[secondary_dataset_titles[i]].forEach(marker => {
      Object.keys(marker).forEach(k => {
        if (k === secondary_dataset_position_colname) {
          secondary_positions[secondary_dataset_titles[i]].push(+marker[k]);
        }
        else if (k === secondary_dataset_pval_colname) {
          secondary_log10_pvalues[secondary_dataset_titles[i]].push(-Math.log10(+marker[k]));
        }
        else if (k === secondary_dataset_snp_colname) {
          secondary_snps[secondary_dataset_titles[i]].push(marker[k]);
        }
      });
    });
    secondary_line_traces[secondary_dataset_titles[i]] = smoothing(secondary_positions[secondary_dataset_titles[i]], secondary_log10_pvalues[secondary_dataset_titles[i]],
      [startbp, endbp], eqtl_smoothing_window_size);
  }

  // Assign each SNP to an LD group:
  for (i = 0; i < ld_values.length; i++) {
    if (snps[i] === lead_snp) {
      lead_snp_index = i;
      ld_colors[i] = lead_snp_color;
    }
    else if (ld_values[i] === -1 || ld_values[i] < 0) {
      no_ld_info_snps.push(i);
      ld_colors[i] = no_ld_info_snps_color;
    }
    else if (Math.abs(ld_values[i]) < 0.2) {
      ld_lt_20_group.push(i);
      ld_colors[i] = ld_lt_20_group_color;
    }
    else if (Math.abs(ld_values[i]) >= 0.2 && Math.abs(ld_values[i]) < 0.4) {
      ld_20_group.push(i);
      ld_colors[i] = ld_20_group_color;
    }
    else if (Math.abs(ld_values[i]) >= 0.4 && Math.abs(ld_values[i]) < 0.6) {
      ld_40_group.push(i);
      ld_colors[i] = ld_40_group_color;
    }
    else if (Math.abs(ld_values[i]) >= 0.6 && Math.abs(ld_values[i]) < 0.8) {
      ld_60_group.push(i);
      ld_colors[i] = ld_60_group_color;
    }
    else if (Math.abs(ld_values[i]) >= 0.8) {
      ld_80_group.push(i);
      ld_colors[i] = ld_80_group_color;
    }
    else if (snps[i] == lead_snp) {
      lead_snp_index = i;
      ld_colors[i] = lead_snp_color;
    }
    else {
      no_ld_info_snps.push(i);
      ld_colors[i] = no_ld_info_snps_color;
    }
  }

  // plot the 7 LD groups
  var no_ld_trace = {
    x: no_ld_info_snps.map(i => positions[i]),
    y: no_ld_info_snps.map(i => log10pvalues[i]),
    name: 'No LD Info',
    mode: 'markers',
    type: 'scatter',
    text: no_ld_info_snps.map(i => snps[i]),
    marker: {
      size: markersize,
      color: no_ld_info_snps_color
    }
  };

  var ld_lt_20_trace = {
    x: ld_lt_20_group.map(i => positions[i]),
    y: ld_lt_20_group.map(i => log10pvalues[i]),
    name: '&lt; 0.2',
    mode: 'markers',
    type: 'scatter',
    text: ld_lt_20_group.map(i => snps[i]),
    marker: {
      size: markersize,
      color: ld_lt_20_group_color
    }
  };

  var ld_20_trace = {
    x: ld_20_group.map(i => positions[i]),
    y: ld_20_group.map(i => log10pvalues[i]),
    name: '0.2',
    mode: 'markers',
    type: 'scatter',
    text: ld_20_group.map(i => snps[i]),
    marker: {
      size: markersize,
      color: ld_20_group_color
    }
  };

  var ld_40_trace = {
    x: ld_40_group.map(i => positions[i]),
    y: ld_40_group.map(i => log10pvalues[i]),
    name: '0.4',
    mode: 'markers',
    type: 'scatter',
    text: ld_40_group.map(i => snps[i]),
    marker: {
      size: markersize,
      color: ld_40_group_color
    }
  };

  var ld_60_trace = {
    x: ld_60_group.map(i => positions[i]),
    y: ld_60_group.map(i => log10pvalues[i]),
    name: '0.6',
    mode: 'markers',
    type: 'scatter',
    text: ld_60_group.map(i => snps[i]),
    marker: {
      size: markersize,
      color: ld_60_group_color
    }
  };

  var ld_80_trace = {
    x: ld_80_group.map(i => positions[i]),
    y: ld_80_group.map(i => log10pvalues[i]),
    name: '0.8',
    mode: 'markers',
    type: 'scatter',
    text: ld_80_group.map(i => snps[i]),
    marker: {
      size: markersize,
      color: ld_80_group_color
    }
  };

  var lead_snp_trace = {
    x: [positions[lead_snp_index]],
    y: [log10pvalues[lead_snp_index]],
    name: 'Lead SNP',
    mode: 'markers',
    type: 'scatter',
    text: lead_snp,
    marker: {
      size: lead_markersize,
      color: lead_snp_color
    },
    yaxis: 'y1'
  };

  all_traces = [no_ld_trace, ld_lt_20_trace, ld_20_trace, ld_40_trace, ld_60_trace, ld_80_trace, lead_snp_trace];


  // Plot the GTEx lines (gtex_line_traces):
  for (var i = 0; i < gtex_tissues.length; i++) {
    var gtex_tissue_trace = {
      x: gtex_line_traces[gtex_tissues[i]][0],
      y: gtex_line_traces[gtex_tissues[i]][1],
      name: gtex_tissues[i],
      mode: 'lines',
      xaxis: 'x1',
      yaxis: 'y2'
    };
    all_traces.push(gtex_tissue_trace);
  }

  for (var i = 0; i < gtex_tissues.length; i++) {
    var gtex_tissue_trace = {
      x: gtex_positions[gtex_tissues[i]],
      y: gtex_log10_pvalues[gtex_tissues[i]],
      name: gtex_tissues[i],
      mode: 'markers',
      type: 'scatter',
      text: gtex_snps[gtex_tissues[i]],
      marker: {
        size: markersize,
        opacity: 0.3
      },
      xaxis: 'x1',
      yaxis: 'y2',
      visible: 'legendonly'
    };
    all_traces.push(gtex_tissue_trace);
  }

  // Plot secondary dataset lines (secondary_line_traces):
  for (var i = 0; i < secondary_dataset_titles.length; i++) {
    var secondary_trace = {
      x: secondary_line_traces[secondary_dataset_titles[i]][0],
      y: secondary_line_traces[secondary_dataset_titles[i]][1],
      name: secondary_dataset_titles[i],
      mode: 'lines',
      xaxis: 'x1',
      yaxis: 'y2'
    };
    all_traces.push(secondary_trace);
  }

  for (var i = 0; i < secondary_dataset_titles.length; i++) {
    var secondary_trace = {
      x: secondary_positions[secondary_dataset_titles[i]],
      y: secondary_log10_pvalues[secondary_dataset_titles[i]],
      name: secondary_dataset_titles[i],
      mode: 'markers',
      type: 'scatter',
      text: secondary_snps[secondary_dataset_titles[i]],
      marker: {
        size: markersize,
        opacity: 0.3
      },
      xaxis: 'x1',
      yaxis: 'y2',
      visible: 'legendonly'
    };
    all_traces.push(secondary_trace);
  }

  var genenames_trace = {
    x: annotations_x,
    y: annotations_y,
    text: annotations_text,
    type: 'scatter',
    mode: 'markers',
    marker: {
      opacity: 0
    },
    yaxis: 'y1',
    showlegend: false,
    name: 'Gene name'
  }
  all_traces.push(genenames_trace);

  var gwas_ymax = d3.max(log10pvalues);
  var gtex_ymax = getYmax(gtex_line_traces, secondary_line_traces);

  function rect_overlap(rect1, rect2) {
    var overlap = false;
    // If one rectangle is on the left side of other
    rectA_x1 = rect1[0];
    rectA_y1 = rect1[1];
    rectA_x2 = rect1[2];
    rectA_y2 = rect1[3];
    rectB_x1 = rect2[0];
    rectB_y1 = rect2[1];
    rectB_x2 = rect2[2];
    rectB_y2 = rect2[3];

    // from https://stackoverflow.com/questions/306316/determine-if-two-rectangles-overlap-each-other
    if (rectA_x1 <= rectB_x2 && rectA_x2 >= rectB_x1 && rectA_y1 <= rectB_y2 && rectA_y2 >= rectB_y1) {
      return (true);
    }
    return (overlap);
  }

  function curr_rect_overlaps(rect_bin, rect_bins) {
    overlap = false;
    for (var i = 0; i < rect_bins.length; i++) {
      var curr_rect_bin = rect_bins[i];
      if (rect_overlap(rect_bin, curr_rect_bin)) {
        return (true);
      }
    }
    return (overlap);
  }

  var full_y_range = gwas_ymax + extra_y_range + gene_area_height;
  i = 0;
  midx = (genesdata[i]['txStart'] + genesdata[i]['txEnd']) / 2;
  // midy = -(genesdata[i]['geneRow'] * row_height) + row_height/2;
  midy = -(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin;
  thegenename = genesdata[i]['name'];
  var xrefloc = (midx - (startbp - extra_x_range)) / (regionsize + 2 * extra_x_range);
  var x0 = d3.max([xrefloc - ((thegenename.length / 2) * percent_occupied_by_one_char), 0]);
  var x1 = xrefloc + ((thegenename.length / 2) * percent_occupied_by_one_char);
  var tempx2 = (genesdata[i]['txStart'] + genesdata[i]['txEnd']) / 2;
  var y0 = (gene_area_height - (-1 * midy) - font_height) / full_y_range;
  var y1 = (gene_area_height - (-1 * midy)) / full_y_range;
  var first_rect_bin = [x0, y0, x1, y1];

  var rect_bins = [];
  rect_bins = [first_rect_bin];
  var annotations_x2 = [];
  var annotations_y2 = [];
  var annotations_text2 = [];
  var locations = [];
  annotations_x2.push(tempx2);
  annotations_y2.push(-(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin);
  annotations_text2.push(genesdata[i]['name']);
  locations.push('bottom');

  for (var i = 1; i < genesdata.length; i++) {
    midx = (genesdata[i]['txStart'] + genesdata[i]['txEnd']) / 2;
    midy = -(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin;
    thegenename = genesdata[i]['name'];
    var xrefloc = (midx - (startbp - extra_x_range)) / (regionsize + 2 * extra_x_range);
    var x0 = d3.max([xrefloc - ((thegenename.length / 2) * percent_occupied_by_one_char), 0]);
    var x1 = xrefloc + ((thegenename.length / 2) * percent_occupied_by_one_char);
    var tempx2 = (genesdata[i]['txStart'] + genesdata[i]['txEnd']) / 2;
    var y0 = (gene_area_height - (-1 * midy) - font_height) / full_y_range;
    var y1 = (gene_area_height - (-1 * midy)) / full_y_range;
    var curr_rect_bin = [x0, y0, x1, y1];
    if (curr_rect_overlaps(curr_rect_bin, rect_bins)) {
      // try the top area of the gene
      y0 = (gene_area_height - (-1 * midy) - font_height + exon_height) / full_y_range;
      y1 = (gene_area_height - (-1 * midy) + exon_height) / full_y_range;
      curr_rect_bin = [x0, y0, x1, y1];
      if (curr_rect_overlaps(curr_rect_bin, rect_bins)) {
        rect_bins.push([-1, -1, -1, -1]); // don't output the genename text then
        annotations_x2.push(-1);
        annotations_y2.push(-1);
        annotations_text2.push(genesdata[i]['name']);
        locations.push('hidden');
      } else {
        rect_bins.push(curr_rect_bin); // put gene name text at the top of the gene
        annotations_x2.push(tempx2);
        annotations_y2.push(-(genesdata[i]['geneRow'] * row_height) + row_height);
        annotations_text2.push(genesdata[i]['name']);
        locations.push('top');
      }
    } else {
      // default to putting the gene name text at the bottom of the gene
      rect_bins.push(curr_rect_bin);
      annotations_x2.push(tempx2);
      annotations_y2.push(-(genesdata[i]['geneRow'] * row_height) + text_height + gene_margin);
      annotations_text2.push(genesdata[i]['name']);
      locations.push('bottom');
    }

    var final_x = [];
    var final_y = [];
    var final_text = [];
    for (var i = 0; i < locations.length; i++) {
      if (locations[i] !== 'hidden') {
        final_x.push(annotations_x2[i]);
        final_y.push(annotations_y2[i]);
        final_text.push(annotations_text2[i]);
      }
    }
    var genenames_trace2 = {
      x: final_x,
      y: final_y,
      text: final_text,
      type: 'scatter',
      mode: 'markers+text',
      marker: {
        opacity: 0
      },
      yaxis: 'y1',
      showlegend: false,
      name: 'Gene name',
      textposition: 'bottom'
    }
  }
  var genenames_trace2 = {
    x: final_x,
    y: final_y,
    text: final_text,
    type: 'scatter',
    mode: 'markers+text',
    marker: {
      opacity: 0
    },
    yaxis: 'y1',
    showlegend: false,
    name: 'Gene name',
    textposition: 'bottom',
    font: { style: 'italic' }
  }
  all_traces.push(genenames_trace2);

  // Shade the Simple Sum Region
  var SS_shade_shape = {
    type: 'rect',
    xref: 'x',
    yref: 'y',
    x0: SS_start,
    y0: 0,
    x1: SS_end,
    y1: gwas_ymax,
    fillcolor: "#d3d3d3",
    opacity: 0.3,
    line: { width: 0 }
  }
  rectangle_shapes.push(SS_shade_shape);


  var layout = {
    xaxis: {
      range: [startbp - extra_x_range, endbp + extra_x_range],
      zeroline: false,
      title: { text: `Chromosome ${chrom} (${build})` }
    },
    yaxis: {
      range: [0 - gene_area_height, gwas_ymax + extra_y_range],
      title: { text: `GWAS -log10(p-value)` }
    },
    yaxis2: {
      range: [0 - gene_area_height * (gtex_ymax / gwas_ymax), gtex_ymax + extra_y_range * (gtex_ymax / gwas_ymax)],
      overlaying: 'y',
      anchor: 'x',
      side: 'right',
      showgrid: false,
      title: 'Secondary datasets -log10(p-value)'
    },
    height: inputHeight,
    width: inputWidth,
    showlegend: true,
    legend: {
      x: 1 + legend_offset,
      y: 1,
      font: { size: font_size }
    },
    zeroline: true,
    hovermode: "closest",
    shapes: rectangle_shapes,
    font: { size: font_size }
  };

  var img_svg = d3.select("#svg-try");
  Plotly.newPlot('plot', all_traces, layout);

}
