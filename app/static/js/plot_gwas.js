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


function prepareXQTLData(data, names, gwas, cfg, cols = {}) {
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

  const gtex = prepareXQTLData(data, data.gtex_tissues, gwas, cfg);
  const mqtl = prepareXQTLData(data?.["xqtl"] || {}, data?.["xqtl_names"] || [], gwas, cfg);
  const secondary = prepareSecondaryData(data, gwas, cfg);

  const traces = [
    ...makeLDTraces(gwas, cfg),
    ...makeEQTLTraces(gtex),
    ...makeEQTLMarkerTraces(gtex),
    ...makeEQTLTraces(secondary),
    ...makeEQTLMarkerTraces(secondary),
    ...makeEQTLTraces(mqtl),
    ...makeEQTLMarkerTraces(mqtl),
    ...makeGeneNameTraces(genes)
  ];

  const gtexYMax = getSecondaryYMax(gtex, secondary);
  const layout = makeLayout(gwas, genes, gtex, secondary, cfg, gtexYMax);

  Plotly.newPlot("plot", traces, layout);
}
