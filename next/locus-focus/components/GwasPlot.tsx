"use client";

import dynamic from "next/dynamic";

import { useEffect, useMemo, useState } from "react";
import {
  Checkbox,
  FormControlLabel,
  Grid,
  MenuItem,
  TextField,
} from "@mui/material";
import { Layout, PlotData, Shape } from "plotly.js";
import { GeneRecord, SessionFile } from "@/lib/ts/types";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false });

interface GwasGenesData extends GeneRecord {
  geneRow?: number;
}

const rectOverlap = (rect1: number[], rect2: number[]) => {
  const overlap = false;
  // If one rectangle is on the left side of other
  const rectAX1 = rect1[0];
  const rectAY1 = rect1[1];
  const rectAX2 = rect1[2];
  const rectAY2 = rect1[3];
  const rectBX1 = rect2[0];
  const rectBY1 = rect2[1];
  const rectBX2 = rect2[2];
  const rectBY2 = rect2[3];

  // from https://stackoverflow.com/questions/306316/determine-if-two-rectangles-overlap-each-other
  if (
    rectAX1 <= rectBX2 &&
    rectAX2 >= rectBX1 &&
    rectAY1 <= rectBY2 &&
    rectAY2 >= rectBY1
  ) {
    return true;
  }
  return overlap;
};

const getLayoutsAndTraces = (
  data: SessionFile,
  genesdata: GwasGenesData[],
  eqtlSmoothingWindowSize: number,
  percentOccupiedByOneChar: number,
  inputHeight: number,
  inputWidth: number,
  fontSize: number,
  legendOffset: number,
  pvalFilter: boolean
) => {
  let positions = data.positions;
  let pvalues = data.pvalues;
  let ldValues = data.ld_values;
  let chrom = data.chrom.toString();
  if (chrom === "23") chrom = "X";
  const startbp = data.startbp;
  const endbp = data.endbp;
  let snps = data.snps;
  const leadSnp = data.lead_snp;
  const gtexTissues = data.gtex_tissues;
  const ssStart = +data.SS_region[0];
  const ssEnd = +data.SS_region[1];
  const build = data.coordinate;
  const secondaryDatasetTitles = data.secondary_dataset_titles;
  const secondaryDatasetPositionColname = data.secondary_dataset_colnames[1];
  const secondaryDatasetSnpColname = data.secondary_dataset_colnames[2];
  const secondaryDatasetPvalColname = data.secondary_dataset_colnames[3];
  if (pvalFilter) {
    const pindices: number[] = [];
    pvalues.forEach((p, i) => {
      if (p < 0.1) pindices.push(i);
    });
    pvalues = pindices.map((i) => pvalues[i]);
    positions = pindices.map((i) => positions[i]);
    ldValues = pindices.map((i) => ldValues[i]);
    snps = pindices.map((i) => snps[i]);
  }
  const log10pvalues = pvalues.map((p) => -Math.log10(p));
  const noLdInfoSnps = [],
    noLdInfoSnpsColor = "#7f7f7f"; // grey
  const ldLt_20Group = [],
    ldLt_20GroupColor = "#1f77b4"; // very low ld points (dark blue)
  const ld_20Group = [],
    ld_20GroupColor = "#17becf"; // low ld points (light blue)
  const ld_40Group = [],
    ld_40GroupColor = "#bcbd22"; // green
  const ld_60Group = [],
    ld_60GroupColor = "#ff7f0e"; // orange
  const ld_80Group = [],
    ld_80GroupColor = "#d62728"; // red
  const leadSnpIndex = 0,
    leadSnpColor = "#9467bd"; // purple
  const markersize = 10;
  const leadMarkersize = 10 * 1.5;
  const ldColors = [];
  const regionsize = endbp - startbp;
  const log10pvalueRange =
    Math.max(...pvalues.map((p) => -Math.log10(p))) -
    Math.min(...pvalues.map((p) => -Math.log10(p)));
  const extraXRange = 0.05 * regionsize;
  const extraYRange = 0.05 * log10pvalueRange;
  const eqtlWindowMultiplier = 150;
  if (eqtlSmoothingWindowSize === -1) {
    eqtlSmoothingWindowSize = (regionsize / 1000000) * eqtlWindowMultiplier;
  }

  const fontHeight = 0.5;

  // Helper functions:
  function smoothing(
    x: number[],
    y: number[],
    xrange: number[],
    windowSize: number
  ) {
    const windowPartition = windowSize;
    const windowing = (xrange[1] - xrange[0]) / windowPartition;

    let curr = xrange[0];
    const smoothCurveX = [];
    const smoothCurveY = [];
    let indices: number[] = [];
    x.map((v, i) => {
      if (v >= curr && v <= curr + windowing) indices.push(i);
    });

    // console.log(x);
    while (indices.length == 0 && curr < xrange[1]) {
      // console.log(curr);
      curr = curr + windowing + 1;
      x.map((v, i) => {
        if (v >= curr && v <= curr + windowing) indices.push(i);
      });
    }

    while (curr < xrange[1]) {
      const yAtIndices = indices.map((i) => y[i]);
      let ymaxAtIndices = Math.max(...yAtIndices);
      if (ymaxAtIndices === -1 || ymaxAtIndices < 0) ymaxAtIndices = 0;
      const desiredXindex: number[] = [];
      indices.forEach((i) => {
        if (y[i] === ymaxAtIndices) desiredXindex.push(i);
      });
      smoothCurveX.push(x[desiredXindex[0]]);
      smoothCurveY.push(ymaxAtIndices);
      curr = curr + windowing + 1;
      indices = [];
      x.map((v, i) => {
        if (v >= curr && v <= curr + windowing) indices.push(i);
      });
      while (indices.length == 0 && curr < xrange[1]) {
        curr = curr + windowing + 1;
        x.map((v, i) => {
          if (v >= curr && v <= curr + windowing) indices.push(i);
        });
      }
    }
    return [smoothCurveX, smoothCurveY];
  }

  function getYmax(
    // gtexTraces: Record<string, Partial<PlotData>>,
    gtexTraces: any,
    secondaryTraces: any
  ) {
    let ymax = 1;
    for (let i = 0; i < gtexTissues.length; i++) {
      const currymax = Math.max(...gtexTraces[gtexTissues[i]][1]);
      if (currymax > ymax) ymax = currymax;
    }
    for (let i = 0; i < secondaryDatasetTitles.length; i++) {
      const currymax = Math.max(
        ...secondaryTraces[secondaryDatasetTitles[i]][1]
      );
      if (currymax > ymax) ymax = currymax;
    }
    return ymax;
  }

  // Add the row number each gene will be plotted into:
  genesdata.sort((a: GeneRecord, b: GeneRecord) => a.txStart - b.txStart);

  function overlap(bin1: number[], bin2: number[]) {
    return (
      (bin2[0] >= bin1[0] && bin2[0] <= bin1[1]) ||
      (bin2[1] >= bin1[0] && bin2[1] <= bin1[1])
    );
  }

  function checkSpace(geneRows: number[][][], rowNum: number, bin: number[]) {
    let occupied = false;
    for (let i = 0; i < geneRows[rowNum].length; i++) {
      if (overlap(geneRows[rowNum][i], bin)) {
        occupied = true;
        return occupied;
      }
    }
    return occupied;
  }

  const firstBin = [genesdata[0].txStart, genesdata[0].txEnd];
  const geneRows = [];
  geneRows[0] = [firstBin];
  genesdata[0].geneRow = 1;

  for (let i = 1; i < genesdata.length; i++) {
    let currRow = 0;
    const geneBin = [genesdata[i].txStart, genesdata[i].txEnd];
    let occupied = checkSpace(geneRows, currRow, geneBin);
    while (occupied && currRow < geneRows.length) {
      currRow++;
      if (currRow < geneRows.length)
        occupied = checkSpace(geneRows, currRow, geneBin);
    }
    if (currRow === geneRows.length) {
      geneRows[currRow] = [geneBin];
      genesdata[i].geneRow = currRow + 1;
    } else {
      geneRows[currRow].push(geneBin);
      genesdata[i].geneRow = currRow + 1;
    }
  }

  const geneAreaPercentage = 0.25;
  const maxNumGeneRows = 7;
  const geneAreaHeight = Math.min(
    geneAreaPercentage * log10pvalueRange * geneRows.length,
    geneAreaPercentage * maxNumGeneRows * log10pvalueRange
  );
  const rowHeight = geneAreaHeight / geneRows.length;
  const textHeight = rowHeight * 0.15;
  const geneMargin = rowHeight * 0.15;
  const exonHeight = rowHeight - 2 * (geneMargin + textHeight);
  const intronHeight = exonHeight * 0.4;

  const rectangleShapes: Partial<Shape>[] = [];
  const annotationsX = [];
  const annotationsY = [];
  const annotationsText = [];
  for (let i = 0; i < genesdata.length; i++) {
    // build intron rectangle shapes for each gene:
    const rectangleShape: Partial<Shape> = {
      type: "rect",
      xref: "x",
      yref: "y",
      x0: genesdata[i].txStart,
      y0:
        -(genesdata[i].geneRow! * rowHeight) +
        textHeight +
        geneMargin +
        (exonHeight - intronHeight) / 2,
      x1: genesdata[i].txEnd,
      y1:
        -(genesdata[i].geneRow! * rowHeight) +
        textHeight +
        geneMargin +
        (exonHeight - intronHeight) / 2 +
        intronHeight,
      line: {
        color: "rgb(55, 128, 191)",
        width: 1,
      },
      fillcolor: "rgba(55, 128, 191, 1)",
    };
    rectangleShapes.push(rectangleShape);
    annotationsX.push(genesdata[i].txStart);
    annotationsX.push((genesdata[i].txStart + genesdata[i].txEnd) / 2);
    annotationsX.push(genesdata[i].txEnd);
    const y =
      -(genesdata[i].geneRow! * rowHeight) +
      textHeight +
      geneMargin +
      (exonHeight - intronHeight) / 2 +
      intronHeight / 2;
    annotationsY.push(y);
    annotationsY.push(y);
    annotationsY.push(y);
    annotationsText.push(genesdata[i].name);
    annotationsText.push(genesdata[i].name);
    annotationsText.push(genesdata[i].name);
    for (let j = 0; j < genesdata[i].exonStarts.length; j++) {
      // build exon rectangle shapes for current gene
      const rectangleShape: Partial<Shape> = {
        type: "rect",
        xref: "x",
        yref: "y",
        x0: genesdata[i].exonStarts[j],
        y0: -(genesdata[i].geneRow! * rowHeight) + textHeight + geneMargin,
        x1: genesdata[i].exonEnds[j],
        y1:
          -(genesdata[i].geneRow! * rowHeight) +
          textHeight +
          geneMargin +
          exonHeight,
        line: {
          color: "rgb(55, 128, 191)",
          width: 1,
        },
        fillcolor: "rgba(55, 128, 191, 1)",
      };
      rectangleShapes.push(rectangleShape);
    }
  }

  // Smooth out each GTEx tissue's association results for plotting as lines:
  const gtexLineTraces: Record<string, number[][]> = {};
  const gtexPositions: Record<string, number[]> = {};
  const gtexLog10Pvalues: Record<string, number[]> = {};
  const gtexSnps: Record<string, number[]> = {};
  for (let i = 0; i < gtexTissues.length; i++) {
    gtexPositions[gtexTissues[i]] = [];
    gtexLog10Pvalues[gtexTissues[i]] = [];
    gtexSnps[gtexTissues[i]] = [];
    (
      data[gtexTissues[i] as keyof SessionFile]! as unknown as Record<
        string,
        any
      >[]
    ).forEach((eqtl) => {
      Object.keys(eqtl).forEach((k) => {
        if (k === "constiantPos") {
          gtexPositions[gtexTissues[i]].push(+eqtl[k]);
        } else if (k === "pval") {
          gtexLog10Pvalues[gtexTissues[i]].push(-Math.log10(+eqtl[k]));
        } else if (k === "rsId") {
          gtexSnps[gtexTissues[i]].push(eqtl[k]);
        }
      });
    });
    gtexLineTraces[gtexTissues[i]] = smoothing(
      gtexPositions[gtexTissues[i]],
      gtexLog10Pvalues[gtexTissues[i]],
      [startbp, endbp],
      eqtlSmoothingWindowSize
    );
  }

  const secondaryLineTraces: Record<string, number[][]> = {};
  const secondaryPositions: Record<string, number[]> = {};
  const secondaryLog10Pvalues: Record<string, number[]> = {};
  const secondarySnps: Record<string, number[]> = {};
  for (let i = 0; i < secondaryDatasetTitles.length; i++) {
    secondaryPositions[secondaryDatasetTitles[i]] = [];
    secondaryLog10Pvalues[secondaryDatasetTitles[i]] = [];
    secondarySnps[secondaryDatasetTitles[i]] = [];
    (
      data[gtexTissues[i] as keyof SessionFile]! as unknown as Record<
        string,
        number
      >[]
    ).forEach((marker) => {
      Object.keys(marker).forEach((k) => {
        if (k === secondaryDatasetPositionColname) {
          secondaryPositions[secondaryDatasetTitles[i]].push(+marker[k]);
        } else if (k === secondaryDatasetPvalColname) {
          secondaryLog10Pvalues[secondaryDatasetTitles[i]].push(
            -Math.log10(+marker[k])
          );
        } else if (k === secondaryDatasetSnpColname) {
          secondarySnps[secondaryDatasetTitles[i]].push(marker[k]);
        }
      });
    });
    secondaryLineTraces[secondaryDatasetTitles[i]] = smoothing(
      secondaryPositions[secondaryDatasetTitles[i]],
      secondaryLog10Pvalues[secondaryDatasetTitles[i]],
      [startbp, endbp],
      eqtlSmoothingWindowSize
    );
  }

  // Assign each SNP to an LD group:
  for (let i = 0; i < ldValues.length; i++) {
    if (snps[i] === leadSnp) {
      ldColors[i] = leadSnpColor;
    } else if (ldValues[i] === -1 || ldValues[i] < 0) {
      noLdInfoSnps.push(i);
      ldColors[i] = noLdInfoSnpsColor;
    } else if (Math.abs(ldValues[i]) < 0.2) {
      ldLt_20Group.push(i);
      ldColors[i] = ldLt_20GroupColor;
    } else if (Math.abs(ldValues[i]) >= 0.2 && Math.abs(ldValues[i]) < 0.4) {
      ld_20Group.push(i);
      ldColors[i] = ld_20GroupColor;
    } else if (Math.abs(ldValues[i]) >= 0.4 && Math.abs(ldValues[i]) < 0.6) {
      ld_40Group.push(i);
      ldColors[i] = ld_40GroupColor;
    } else if (Math.abs(ldValues[i]) >= 0.6 && Math.abs(ldValues[i]) < 0.8) {
      ld_60Group.push(i);
      ldColors[i] = ld_60GroupColor;
    } else if (Math.abs(ldValues[i]) >= 0.8) {
      ld_80Group.push(i);
      ldColors[i] = ld_80GroupColor;
    } else if (snps[i] == leadSnp) {
      ldColors[i] = leadSnpColor;
    } else {
      noLdInfoSnps.push(i);
      ldColors[i] = noLdInfoSnpsColor;
    }
  }

  // plot the 7 LD groups
  const noLdTrace: Partial<PlotData> = {
    x: noLdInfoSnps.map((i) => positions[i]),
    y: noLdInfoSnps.map((i) => log10pvalues[i]),
    name: "No LD Info",
    mode: "markers",
    type: "scatter",
    text: noLdInfoSnps.map((i) => snps[i]),
    marker: {
      size: markersize,
      color: noLdInfoSnpsColor,
    },
  };

  const ldLt_20Trace: Partial<PlotData> = {
    x: ldLt_20Group.map((i) => positions[i]),
    y: ldLt_20Group.map((i) => log10pvalues[i]),
    name: "&lt; 0.2",
    mode: "markers",
    type: "scatter",
    text: ldLt_20Group.map((i) => snps[i]),
    marker: {
      size: markersize,
      color: ldLt_20GroupColor,
    },
  };

  const ld_20Trace: Partial<PlotData> = {
    x: ld_20Group.map((i) => positions[i]),
    y: ld_20Group.map((i) => log10pvalues[i]),
    name: "0.2",
    mode: "markers",
    type: "scatter",
    text: ld_20Group.map((i) => snps[i]),
    marker: {
      size: markersize,
      color: ld_20GroupColor,
    },
  };

  const ld_40Trace: Partial<PlotData> = {
    x: ld_40Group.map((i) => positions[i]),
    y: ld_40Group.map((i) => log10pvalues[i]),
    name: "0.4",
    mode: "markers",
    type: "scatter",
    text: ld_40Group.map((i) => snps[i]),
    marker: {
      size: markersize,
      color: ld_40GroupColor,
    },
  };

  const ld_60Trace: Partial<PlotData> = {
    x: ld_60Group.map((i) => positions[i]),
    y: ld_60Group.map((i) => log10pvalues[i]),
    name: "0.6",
    mode: "markers",
    type: "scatter",
    text: ld_60Group.map((i) => snps[i]),
    marker: {
      size: markersize,
      color: ld_60GroupColor,
    },
  };

  const ld_80Trace: Partial<PlotData> = {
    x: ld_80Group.map((i) => positions[i]),
    y: ld_80Group.map((i) => log10pvalues[i]),
    name: "0.8",
    mode: "markers",
    type: "scatter",
    text: ld_80Group.map((i) => snps[i]),
    marker: {
      size: markersize,
      color: ld_80GroupColor,
    },
  };

  const leadSnpTrace: Partial<PlotData> = {
    x: [positions[leadSnpIndex]],
    y: [log10pvalues[leadSnpIndex]],
    name: "Lead SNP",
    mode: "markers",
    type: "scatter",
    text: leadSnp,
    marker: {
      size: leadMarkersize,
      color: leadSnpColor,
    },
    yaxis: "y1",
  };

  const allTraces: Partial<PlotData>[] = [
    noLdTrace,
    ldLt_20Trace,
    ld_20Trace,
    ld_40Trace,
    ld_60Trace,
    ld_80Trace,
    leadSnpTrace,
  ];

  // Plot the GTEx lines (gtexLineTraces):
  for (let i = 0; i < gtexTissues.length; i++) {
    const gtexTissueTrace: Partial<PlotData> = {
      x: gtexLineTraces[gtexTissues[i]][0],
      y: gtexLineTraces[gtexTissues[i]][1],
      name: gtexTissues[i],
      mode: "lines",
      xaxis: "x1",
      yaxis: "y2",
    };
    allTraces.push(gtexTissueTrace);
  }

  for (let i = 0; i < gtexTissues.length; i++) {
    const gtexTissueTrace: Partial<PlotData> = {
      x: gtexPositions[gtexTissues[i]],
      y: gtexLog10Pvalues[gtexTissues[i]],
      name: gtexTissues[i],
      mode: "markers",
      type: "scatter",
      text: gtexSnps[gtexTissues[i]].toString(),
      marker: {
        size: markersize,
        opacity: 0.3,
      },
      xaxis: "x1",
      yaxis: "y2",
      visible: "legendonly",
    };
    allTraces.push(gtexTissueTrace);
  }

  // Plot secondary dataset lines (secondaryLineTraces):
  for (let i = 0; i < secondaryDatasetTitles.length; i++) {
    const secondaryTrace: Partial<PlotData> = {
      x: secondaryLineTraces[secondaryDatasetTitles[i]][0],
      y: secondaryLineTraces[secondaryDatasetTitles[i]][1],
      name: secondaryDatasetTitles[i],
      mode: "lines",
      xaxis: "x1",
      yaxis: "y2",
    };
    allTraces.push(secondaryTrace);
  }

  for (let i = 0; i < secondaryDatasetTitles.length; i++) {
    const secondaryTrace: Partial<PlotData> = {
      x: secondaryPositions[secondaryDatasetTitles[i]],
      y: secondaryLog10Pvalues[secondaryDatasetTitles[i]],
      name: secondaryDatasetTitles[i],
      mode: "markers",
      type: "scatter",
      text: secondarySnps[secondaryDatasetTitles[i]].toString(),
      marker: {
        size: markersize,
        opacity: 0.3,
      },
      xaxis: "x1",
      yaxis: "y2",
      visible: "legendonly",
    };
    allTraces.push(secondaryTrace);
  }

  const genenamesTrace: Partial<PlotData> = {
    x: annotationsX,
    y: annotationsY,
    text: annotationsText,
    type: "scatter",
    mode: "markers",
    marker: {
      opacity: 0,
    },
    yaxis: "y1",
    showlegend: false,
    name: "Gene name",
  };
  allTraces.push(genenamesTrace);

  const gwasYmax = Math.max(...log10pvalues);
  const gtexYmax = getYmax(gtexLineTraces, secondaryLineTraces);

  function currRectOverlaps(rectBin: number[], rectBins: number[][]) {
    for (let i = 0; i < rectBins.length; i++) {
      const currRectBin = rectBins[i];
      if (rectOverlap(rectBin, currRectBin)) {
        return true;
      }
    }
    return false;
  }

  const fullYRange = gwasYmax + extraYRange + geneAreaHeight;
  const i = 0;
  const midx = (genesdata[i].txStart + genesdata[i].txEnd) / 2;
  const midy = -(genesdata[i].geneRow * rowHeight) + textHeight + geneMargin;
  const thegenename = genesdata[i].name;
  const xrefloc =
    (midx - (startbp - extraXRange)) / (regionsize + 2 * extraXRange);
  const x0 = Math.max(
    xrefloc - (thegenename.length / 2) * percentOccupiedByOneChar,
    0
  );
  const x1 = xrefloc + (thegenename.length / 2) * percentOccupiedByOneChar;
  const tempx2 = (genesdata[i].txStart + genesdata[i].txEnd) / 2;

  const y0 = (geneAreaHeight - -1 * midy - fontHeight) / fullYRange;
  const y1 = (geneAreaHeight - -1 * midy) / fullYRange;
  const firstRectBin = [x0, y0, x1, y1];

  const rectBins = [firstRectBin];
  const annotationsX2 = [];
  const annotationsY2 = [];
  const annotationsText2 = [];
  const locations = [];
  annotationsX2.push(tempx2);
  annotationsY2.push(
    -(genesdata[i].geneRow * rowHeight) + textHeight + geneMargin
  );
  annotationsText2.push(genesdata[i].name);
  locations.push("bottom");

  const layouts = [];
  for (let i = 1; i < genesdata.length; i++) {
    const finalX = [];
    const finalY = [];
    const finalText = [];
    const midx = (genesdata[i].txStart + genesdata[i].txEnd) / 2;
    const midy = -(genesdata[i].geneRow! * rowHeight) + textHeight + geneMargin;
    const thegenename = genesdata[i].name;
    const xrefloc =
      (midx - (startbp - extraXRange)) / (regionsize + 2 * extraXRange);
    const x0 = Math.max(
      xrefloc - (thegenename.length / 2) * percentOccupiedByOneChar,
      0
    );
    const x1 = xrefloc + (thegenename.length / 2) * percentOccupiedByOneChar;
    const tempx2 = (genesdata[i].txStart + genesdata[i].txEnd) / 2;

    const y0 = (geneAreaHeight - -1 * midy - fontHeight) / fullYRange;
    const y1 = (geneAreaHeight - -1 * midy) / fullYRange;
    const currRectBin: number[] = [x0, y0, x1, y1];

    if (currRectOverlaps(currRectBin, rectBins)) {
      // try the top area of the gene
      const y0 =
        (geneAreaHeight - -1 * midy - fontHeight + exonHeight) / fullYRange;
      const y1 = (geneAreaHeight - -1 * midy + exonHeight) / fullYRange;
      const currRectBin = [x0, y0, x1, y1];
      if (currRectOverlaps(currRectBin, rectBins)) {
        rectBins.push([-1, -1, -1, -1]); // don't output the genename text then
        annotationsX2.push(-1);
        annotationsY2.push(-1);
        annotationsText2.push(genesdata[i].name);
        locations.push("hidden");
      } else {
        rectBins.push(currRectBin); // put gene name text at the top of the gene
        annotationsX2.push(tempx2);
        annotationsY2.push(-(genesdata[i].geneRow! * rowHeight) + rowHeight);
        annotationsText2.push(genesdata[i].name);
        locations.push("top");
      }
    } else {
      // default to putting the gene name text at the bottom of the gene
      rectBins.push(currRectBin);
      annotationsX2.push(tempx2);
      annotationsY2.push(
        -(genesdata[i].geneRow! * rowHeight) + textHeight + geneMargin
      );
      annotationsText2.push(genesdata[i].name);
      locations.push("bottom");
    }

    for (let i = 0; i < locations.length; i++) {
      if (locations[i] !== "hidden") {
        finalX.push(annotationsX2[i]);
        finalY.push(annotationsY2[i]);
        finalText.push(annotationsText2[i]);
      }
    }

    if (i === genesdata.length - 1) {
      const genenamesTrace2: Partial<PlotData> = {
        x: finalX,
        y: finalY,
        text: finalText,
        type: "scatter",
        mode: "text+markers",
        marker: {
          opacity: 0,
        },
        yaxis: "y1",
        showlegend: false,
        name: "Gene name",
        textposition: "bottom center",
        textfont: { style: "italic" },
      };
      allTraces.push(genenamesTrace2);
    }

    // Shade the Simple Sum Region
    const ssShadeShape: Partial<Shape> = {
      type: "rect",
      xref: "x",
      yref: "y",
      x0: ssStart,
      y0: 0,
      x1: ssEnd,
      y1: gwasYmax,
      fillcolor: "#d3d3d3",
      opacity: 0.05,
      line: { width: 0 },
    };
    rectangleShapes.push(ssShadeShape);

    const layout: Partial<Layout> = {
      xaxis: {
        range: [startbp - extraXRange, endbp + extraXRange],
        zeroline: false,
        title: { text: `Chromosome ${chrom} (${build})` },
      },
      yaxis: {
        range: [0 - geneAreaHeight, gwasYmax + extraYRange],
        title: { text: `GWAS -log10(p-value)` },
      },
      yaxis2: {
        range: [
          0 - geneAreaHeight * (gtexYmax / gwasYmax),
          gtexYmax + extraYRange * (gtexYmax / gwasYmax),
        ],
        overlaying: "y",
        anchor: "x",
        side: "right",
        showgrid: false,
        title: { text: "Secondary datasets -log10(p-value)" },
      },
      height: inputHeight,
      width: inputWidth!,
      showlegend: true,
      legend: {
        x: 1 + legendOffset,
        y: 1,
        font: { size: fontSize },
      },
      hovermode: "closest",
      shapes: rectangleShapes,
      font: { size: fontSize },
    };

    layouts.push(layout);
  }

  return { layouts, allTraces };
};

const GwasPlot: React.FC<{
  data: SessionFile;
  genesdata: GwasGenesData[];
  sessionId: string;
  eqtlSmoothingWindowSize?: number;
  percentOccupiedByOneChar?: number;
  inputHeight?: number;
  inputWidth?: number;
  fontSize?: number;
  legendOffset?: number;
}> = ({
  data,
  genesdata,
  sessionId,
  eqtlSmoothingWindowSize,
  percentOccupiedByOneChar,
  inputHeight,
  inputWidth,
  fontSize,
  legendOffset,
}) => {
  const [pvalueFilter, setPvalueFilter] = useState(true);

  const geneNames = genesdata
    .map((g) => g.name)
    .filter(Boolean)
    .sort();

  const [selectedGene, setSelectedGene] = useState(geneNames[0]);

  useEffect(() => {
    fetch(
      `${process.env.NEXT_PUBLIC_BROWSER_API_HOST}/update/${sessionId}/${selectedGene}`
    ).then((r) =>
      r
        .json()
        .then((r) =>
          getLayoutsAndTraces(
            r,
            genesdata,
            eqtlSmoothingWindowSize || -1,
            percentOccupiedByOneChar || 0.02,
            inputHeight || 720,
            inputWidth || 1080,
            fontSize || 14,
            legendOffset || 0.1,
            pvalueFilter
          )
        )
    );
  }, [selectedGene, data.sessionid]);

  const { layouts, allTraces: traces } = useMemo(
    () =>
      getLayoutsAndTraces(
        data,
        genesdata,
        eqtlSmoothingWindowSize || -1,
        percentOccupiedByOneChar || 0.02,
        inputHeight || 720,
        inputWidth || 1080,
        fontSize || 14,
        legendOffset || 0.1,
        pvalueFilter
      ),
    [data, genesdata, pvalueFilter]
  );

  return (
    <Grid spacing={2} container direction="column">
      <Grid>
        <TextField
          select
          label="Select a gene"
          value={selectedGene}
          onChange={(e) => setSelectedGene(e.target.value)}
        >
          {geneNames.map((g) => (
            <MenuItem key={g} value={g}>
              {g}
            </MenuItem>
          ))}
        </TextField>
      </Grid>
      <Grid>
        <FormControlLabel
          control={
            <Checkbox
              checked={pvalueFilter}
              onChange={(e) => setPvalueFilter(e.currentTarget.checked)}
            />
          }
          label="Uncheck to draw all GWAS points below -log10P < 1 (slow)"
        />
      </Grid>
      <Grid>
        {/* It appears all the layouts are the same? */}
        <Plot layout={layouts[0]} data={traces} />
      </Grid>
    </Grid>
  );
};

export default GwasPlot;
