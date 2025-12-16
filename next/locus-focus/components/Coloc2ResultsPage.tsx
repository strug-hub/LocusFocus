"use client";

import dynamic from "next/dynamic";

import { useEffect, useMemo, useState } from "react";
import { json as fetchJson } from "d3-fetch";
import { Button, Grid, Link, Paper, Typography } from "@mui/material";
import { GridColDef } from "@mui/x-data-grid";
import { useRouter, useSearchParams } from "next/navigation";
import { PlotData } from "plotly.js";
import { Coloc2AnalysisResults, SessionFile } from "@/lib/ts/types";
import { GwasPlot, InfoGrid, LFDataTable, LoadingOverlay } from "@/components";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false });

const getParamsRows = (
  data: SessionFile,
  type: "default" | "set-based-test"
) => {
  const snpCountHeader = `Number of SNPs in ${data.chrom}:${data.startbp}-${data.endbp}`;
  const ssSnpCountHeader = `Number of SNPs in ${data.chrom}:${data.SS_region[0]}-${data.SS_region[1]}`;

  let rows = [];

  if (type === "default") {
    rows = [
      { field: data.sessionid, label: "Session ID" },
      { field: data.lead_snp, label: "Lead SNP" },
      { field: data.chrom, label: "Chromosome" },
      { field: data.startbp, label: "Start position" },
      { field: data.endbp, label: "End position" },
      { field: data.coordinate, label: "Build" },
      { field: data.inferVariant, label: "Infer variants" },
      {
        field: data.snps.length,
        label: snpCountHeader,
      },
      { field: data.ld_populations, label: "LD Population" },
      { field: data.gtex_version, label: "GTEx version" },
      {
        field: data.gtex_tissues.length,
        label: "Number of GTEx tissues selected",
      },
      {
        field: data.gtex_genes.length,
        label: "Number of GTEx genes selected",
      },
      {
        field: data.SS_region,
        label: "SS region",
      },
      {
        field: data.num_SS_snps,
        label: ssSnpCountHeader,
      },
      {
        field: data.set_based_p,
        label: "First stage -log10(SS P-value) threshold",
      },
      { field: data.snp_warning, label: "Many SNPs not matching GTEx SNPs" },
      { field: data.thresh, label: "SNPs matching threshold level" },
      {
        field: data.numGTExMatches,
        label: "Number of SNPs matching with GTEx",
      },
      {
        field: data.secondary_dataset_titles.length,
        label: "Number of user-provided secondary datasetsx",
      },
      { field: data.runcoloc2, label: "Run COLOC2" },
    ];
  } else if (type === "set-based-test") {
    rows = [
      { field: data.sessionid, label: "Session ID" },
      { field: data.coordinate, label: "Build" },
      { field: data.ld_populations, label: "LD Population" },
      {
        field: data.first_stages.length,
        label: "Total set-based tests performed",
      },
      {
        field: data?.regions?.length || 0,
        label: "Number of regions",
      },
      {
        field: data.multiple_tests ? "Yes" : "No",
        label: "Multiple tests?",
      },
    ].concat(
      !!data.multiple_tests
        ? (data?.regions || []).map((r, i) => ({
            label: `SNPs used in region '${r}'`,
            field: data.snps_used_in_test
              ? data.snps_used_in_test[i].length
              : 0,
          }))
        : [
            {
              label: "Number of SNPs used",
              field: data.snps_used_in_test
                ? data.snps_used_in_test[0].length
                : 0,
            },
          ]
    );
  } else {
    // Shouldn't get here
    throw Error(`Unexpected params table type: ${type}`);
  }

  return rows.map((r, id) => ({ id, ...r }));
};

const ColocResultsPage: React.FC = () => {
  const [loading, setLoading] = useState(false);
  const [location, setLocation] = useState("false");
  const [coloc2Results, setColoc2Results] = useState<Coloc2AnalysisResults>();
  const [ssTableType, setSSTableType] = useState<"byGene" | "byTissue">(
    "byGene"
  );

  useEffect(() => {
    setLocation(window.location.href);
  }, []);

  const params = useSearchParams();

  const router = useRouter();

  const firstStageTableData = useMemo(() => {
    if (coloc2Results) {
      const tableconfig: {
        data: Record<string, any>[];
        columns: GridColDef[];
      } = {
        data: [],
        columns: [
          {
            field: "firstStageP",
            headerName: "Set-based test P-value",
            flex: 1,
          },
          {
            field: "positionLength",
            headerName: "Number of SNPs used",
            flex: 1,
          },
        ],
      };
      const sessionData = coloc2Results.sessionfile;
      const positions = sessionData.snps_used_in_test;
      if (sessionData.hasOwnProperty("regions") && !!positions) {
        // set-based only
        const numTests = sessionData.first_stages.length;
        const multipleTestsRequested = sessionData.multiple_tests;
        if (multipleTestsRequested) {
          // multiple tests
          tableconfig.data = sessionData.regions!.map((regiontext, i) => ({
            id: i,
            region: regiontext,
            firstStageP: sessionData.first_stage_Pvalues[i],
            positionLength: positions[i].length,
          }));

          tableconfig.columns = [
            {
              field: "region",
              headerName: "Test Region",
              flex: 1,
            },
          ];
        } else {
          // should be 1 test
          if (numTests > 1)
            console.error(
              `'${numTests}' test detected despite not requesting multiple tests.`
            );
          tableconfig.data = [
            {
              id: 1,
              title: sessionData.dataset_title,
              firstStageP: sessionData.first_stage_Pvalues[0],
              positionLength: positions[0].length,
            },
          ];

          tableconfig.columns.unshift({
            field: "title",
            headerName: "Dataset Description",
            flex: 1,
          });
        }
      } else {
        tableconfig.data = sessionData.secondary_dataset_titles.map(
          (title, i) => ({
            id: i,
            title,
            firstStageP: sessionData.first_stage_Pvalues[i],
          })
        );
        tableconfig.columns.unshift({
          field: "title",
          headerName: "Dataset Description",
          flex: 1,
        });
      }
      return tableconfig;
    }
  }, [coloc2Results]);

  const secondaryTableData = useMemo(() => {
    if (coloc2Results) {
      const sspVals = coloc2Results.SSPvalues_file;
      return sspVals.Secondary_dataset_titles.map((description, i) => ({
        id: i,
        description,
        values: sspVals.SSPvalues_secondary[i],
        count: sspVals.Num_SNPs_Used_for_SS_secondary[i],
      }));
    }
  }, [coloc2Results]);

  const ssGuidanceTableData = useMemo(() => {
    if (coloc2Results) {
      const {
        Genes: genes,
        Tissues: tissues,
        SSPvalues: ssp,
        SSPvalues_secondary: ssp2,
      } = coloc2Results.SSPvalues_file;

      let numGTEx = 0;
      const numSecondary = ssp2.length;
      let numNoeQTL = 0;
      let numFirstStage = 0;
      let numTested = 0;
      let numFailed = 0;
      for (let i = 0; i < tissues.length; i++) {
        // for each tissue
        for (let j = 0; j < genes.length; j++) {
          // for each gene
          numGTEx += 1;
          if (ssp[i][j] == -1) {
            numNoeQTL += 1;
          } else if (ssp[i][j] == -2) {
            numFirstStage += 1;
          } else if (ssp[i][j] == -3) {
            numFailed += 1;
          } else if (ssp[i][j] > 0) {
            numTested += 1;
          } else {
            numFailed += 1;
          }
        }
      }
      for (let i = 0; i < ssp2.length; i++) {
        if (ssp2[i] == -1) {
          numNoeQTL += 1;
        } else if (ssp2[i] == -2) {
          numFirstStage += 1;
        } else if (ssp2[i] == -3) {
          numFailed += 1;
        } else if (ssp2[i] > 0) {
          numTested += 1;
        } else {
          numFailed += 1;
        }
      }

      const suggested_SSP =
        numTested > 0
          ? -Math.log10(0.05 / numTested)
          : "N/A (No datasets were tested successfully)";

      return [
        {
          id: 1,
          field: "Total number of secondary datasets (including GTEx)",
          value: numGTEx + numSecondary,
        },
        {
          id: 2,
          field: "Total number of secondary datasets (including GTEx)",
          value: numGTEx + numSecondary,
        },

        { id: 3, field: "Total number of GTEx datasets", value: numGTEx },
        {
          id: 4,
          field: "Total number of user-uploaded secondary datasets",
          value: numSecondary,
        },
        {
          id: 5,
          field: "Number of datasets with no eQTL data (-1)",
          value: numNoeQTL,
        },
        {
          id: 6,
          field: "Number of datasets not passing first stage (-2)",
          value: numFirstStage,
        },
        {
          id: 7,
          field: "Number of datasets with computation error (-3)",
          value: numFailed,
        },
        {
          id: 8,
          field: "Number of datasets tested for colocalization",
          value: numTested,
        },
        {
          id: 9,
          field:
            "Suggested Simple Sum colocalization threshold at alpha 0.05 (-log10P)",
          value: suggested_SSP,
        },
      ];
    }
  }, [coloc2Results]);

  const ssTableData = useMemo(() => {
    if (coloc2Results) {
      const byTissue = [];
      const bySnpCount: Record<string, string | number>[] = [];
      const byGene = coloc2Results.SSPvalues_file.Genes.map<
        Record<string, number | string>
      >((Gene, id) => ({
        id,
        Gene,
      }));
      for (let i = 0; i < coloc2Results.SSPvalues_file.Tissues.length; i++) {
        const tissueRow: Record<string, string | number> = {
          Tissue: coloc2Results.SSPvalues_file.Tissues[i],
          id: i,
        };
        const snpCountRow = { ...tissueRow };
        for (let j = 0; j < coloc2Results.SSPvalues_file.Genes.length; j++) {
          tissueRow[coloc2Results.SSPvalues_file.Genes[j]] =
            coloc2Results.SSPvalues_file.SSPvalues[i][j];
          snpCountRow[coloc2Results.SSPvalues_file.Genes[j]] =
            coloc2Results.SSPvalues_file.Num_SNPs_Used_for_SS[i][j];
          byGene[j][coloc2Results.SSPvalues_file.Tissues[i]] =
            coloc2Results.SSPvalues_file.SSPvalues[i][j];
        }
        bySnpCount.push(snpCountRow);
        byTissue.push(tissueRow);
      }
      return { byTissue, byGene, bySnpCount };
    }
  }, [coloc2Results]);

  const heatMapInput = useMemo(() => {
    if (coloc2Results) {
      const { gtex_genes, gtex_tissues } = coloc2Results.sessionfile;
      const { SSPvalues, SSPvalues_secondary } = coloc2Results.SSPvalues_file;

      return getHeatmapInput(
        gtex_genes,
        gtex_tissues,
        SSPvalues,
        SSPvalues_secondary
      );
    }
  }, [coloc2Results]);

  const sessionId = useMemo(() => {
    return params.get("sessionId");
  }, [params]);

  const fetchSessionData = async (sessionId: string) => {
    setLoading(true);
    try {
      const sessionData = await fetch(
        `${process.env.NEXT_PUBLIC_BROWSER_API_HOST}/session_id/${sessionId}`
      );

      const endpoints: Record<string, string> = await sessionData.json();
      const dataFetches = Object.entries(endpoints)
        .filter(([k]) => k !== "sessionid")
        .map(async ([key, endpoint]) => {
          const data = await fetchJson(
            `${process.env.NEXT_PUBLIC_BROWSER_API_HOST}/static/${endpoint}`
          );
          return [key, data];
        }) as Promise<[keyof Coloc2AnalysisResults, any]>[];

      const allData = await Promise.all(dataFetches);
      const results = {} as Coloc2AnalysisResults;
      allData.forEach(([key, data]) => {
        results[key] = data;
      });
      setColoc2Results(results);
    } catch (e) {
      console.error(e);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    if (!sessionId) {
      console.error("No session ID!");
      return router.push("/");
    } else {
      fetchSessionData(sessionId);
    }
  }, [router, sessionId]);

  return (
    <Grid container direction="column" marginTop={3} spacing={2}>
      <Grid container direction="row" justifyContent="center">
        <SessionBox sessionId={params.get("sessionId")!} location={location} />
      </Grid>
      <Grid>
        <Typography variant="h5">Params Table</Typography>
      </Grid>
      <Grid width="100%" justifyContent="center" container direction="row">
        <Grid flexGrow={1} sx={{ maxWidth: "950px" }}>
          {!!coloc2Results && (
            <LFDataTable
              columns={[
                { field: "label", headerName: "Field", flex: 1 },
                { field: "field", headerName: "Value", flex: 1 },
              ]}
              rows={getParamsRows(coloc2Results.sessionfile, "default")}
            />
          )}
        </Grid>
      </Grid>
      <Grid container direction="column">
        <Grid>
          <Typography variant="h5">Colocalization Plot</Typography>
        </Grid>
        {!!coloc2Results && !!sessionId && (
          <Grid>
            <GwasPlot
              sessionId={sessionId}
              data={coloc2Results.sessionfile}
              genesdata={coloc2Results.genesfile}
            />
          </Grid>
        )}
      </Grid>
      <Grid container direction="column">
        <Grid>
          <Typography variant="h5">
            Heatmap of Simple Sum Colocalization P-values for GTEx Tissues
            Selected
          </Typography>
        </Grid>
        {!!heatMapInput && (
          <Grid>
            <Plot
              data={heatMapInput.heatmapData}
              layout={heatMapInput.heatmapLayout}
            />
          </Grid>
        )}
      </Grid>
      <Grid container direction="column" spacing={2}>
        <Grid>
          <Typography variant="h5">
            Simple Sum -log10 P-values Table for GTEx Genes/Tissues Selected
          </Typography>
        </Grid>
        <Grid>
          <Button
            onClick={() =>
              setSSTableType(ssTableType === "byGene" ? "byTissue" : "byGene")
            }
            variant="contained"
          >
            Transpose
          </Button>
        </Grid>
        {!!ssTableData && !!ssTableType && (
          <Grid>
            <LFDataTable
              columns={Object.keys(ssTableData[ssTableType][0])
                .filter((k) => k !== "id")
                .map((field) => ({
                  field,
                  flex: 1,
                }))}
              rows={ssTableData[ssTableType]}
            />
          </Grid>
        )}
        <Grid>
          <InfoGrid>
            <Typography>
              -1 values correspond to gene-tissue pairs with no eQTL data
              (likely due to little or no expression).
            </Typography>
            <Typography>
              -2 values correspond to gene-tissue pairs that did not pass the
              Bonferroni-corrected first stage testing for signficance among the
              secondary datasets chosen
            </Typography>
            <Typography>
              -3 values correspond to gene-tissue pairs where the Simple Sum
              P-value computation failed, likely due to insufficient SNPs.
            </Typography>
            <Typography>
              Please note that we leave up to the user to determine the
              threshold of significance among the datasets that passed the
              first-stage significance test. For example, if a user selected 3
              tissues and 4 genes for testing, and 3 other secondary datasets (a
              total of 3 × 4 + 3 = 15 tests) and among these, 6 datasets were
              tested for colocalization, then one would conservatively choose to
              consider a Bonferroni-corrected p-value threshold of 0.05 / 6 =
              8.3 × 10-3. Thus, Simple Sum colocalization tests above this
              threshold would be considered as significant.
            </Typography>
          </InfoGrid>
        </Grid>
        <Grid container direction="column">
          <Grid>
            <Typography variant="h5">
              Number of SNPs Used for Simple Sum Calculation for GTEx
              Genes/Tissues Selected
            </Typography>
          </Grid>
          <Grid>
            {!!ssTableData && (
              <LFDataTable
                columns={Object.keys(ssTableData.bySnpCount[0])
                  .filter((k) => k !== "id")
                  .map((field) => ({
                    field,
                    flex: 1,
                  }))}
                rows={ssTableData.bySnpCount}
              />
            )}
          </Grid>
        </Grid>
        <Grid container direction="column" spacing={2}>
          <Grid>
            <Typography variant="h5">
              Simple Sum -log10 P-values Table for Secondary Datasets Uploaded
            </Typography>
          </Grid>
          <Grid>
            {!!secondaryTableData && (
              <LFDataTable
                columns={[
                  {
                    field: "description",
                    headerName: "Dataset Description",
                    flex: 1,
                  },
                  {
                    field: "value",
                    headerName: "Simple Sum -log10P",
                    flex: 1,
                  },
                  {
                    field: "count",
                    headerName: "Number of SNPs Used in SS Calculation",
                    flex: 1,
                  },
                ]}
                rows={secondaryTableData}
              />
            )}
          </Grid>
          <InfoGrid>
            <Typography>
              -1 values correspond to gene-tissue pairs with no eQTL data
              (likely due to little or no expression)
            </Typography>
            <Typography>
              -2 values correspond to gene-tissue pairs that did not pass the
              Bonferroni-corrected first stage testing for signficance among the
              secondary datasets chosen
            </Typography>
            -3 values correspond to gene-tissue pairs where the Simple Sum
            P-value computation failed, likely due to insufficient SNPs
            <Typography>
              Please note that we leave up to the user to determine the
              threshold of significance among the datasets that passed the
              first-stage significance test. For example, if a user selected 3
              tissues and 4 genes for testing, and 3 other secondary datasets (a
              total of 3 x 4 + 3 = 15 tests) and among these, 6 datasets were
              tested for colocalization, then one would conservatively choose to
              consider a Bonferroni-corrected p-value threshold of 0.05 / 6 =
              8.3 x 10-3. Thus, Simple Sum colocalization tests above this
              threshold would be considered as significant.
            </Typography>
          </InfoGrid>
        </Grid>
      </Grid>
      <Grid container direction="column" spacing={2}>
        <Grid>
          <Typography variant="h5">Simple Sum Guidance Summary</Typography>
        </Grid>
        {!!ssGuidanceTableData && (
          <LFDataTable
            columns={[
              {
                field: "field",
                headerName: "Field",
                flex: 1,
              },
              {
                field: "value",
                headerName: "Value",
                flex: 1,
              },
            ]}
            rows={ssGuidanceTableData}
          />
        )}
      </Grid>
      <Grid container direction="column" spacing={2}>
        <Grid>
          <Typography variant="h5">
            COLOC2 Posterior Probability Results Table
          </Typography>
        </Grid>
        <Grid>
          {coloc2Results?.coloc2_file?.ProbeID?.length ? (
            <LFDataTable
              columns={[
                {
                  field: "ProbeID",
                  headerName: "Probe ID",
                },
                {
                  field: "PPH4abf",
                  headerName: "PP.H4.abf",
                },
              ]}
              rows={coloc2Results?.coloc2_file.ProbeID.map((p, i) => ({
                p,
                ph: coloc2Results?.coloc2_file.PPH4abf[i],
              }))}
            />
          ) : (
            <Typography color="info">
              Insufficient data for COLOC2 calculations
            </Typography>
          )}
        </Grid>
      </Grid>
      <Grid container direction="column" spacing={2}>
        <Grid>
          <Typography variant="h5">
            Set-based Significance Test -log10 P-values Table for Datasets
            Uploaded
          </Typography>
        </Grid>
        <Grid>
          {!!firstStageTableData && (
            <LFDataTable
              columns={firstStageTableData.columns}
              rows={firstStageTableData.data}
            />
          )}
        </Grid>
      </Grid>
      <LoadingOverlay open={loading} />
    </Grid>
  );
};

export default ColocResultsPage;

const normalize = (value: number, min: number, max: number) => {
  return Math.abs((value - min) / (max - min));
};

const getHeatmapInput = (
  genes: string[],
  tissues: string[],
  SSPvalues: number[][],
  SSPvalues_secondary: number[],
  imageWidth = 1080,
  imageHeight = 720,
  fontSize = 12
) => {
  // remember that these are -log10P-values
  // want the different negative p-value statuses to have white/black/grey colors
  // can't really do it that way b/c of the different cases that may occur
  // prior checks ensure that we have at least one positive -log10 SS Pvalue
  let pmax = 0;
  const newSSPvalues = SSPvalues.map((v) => [...v]);
  let numDatasets = 0;
  for (let i = 0; i < tissues.length; i++) {
    for (let j = 0; j < genes.length; j++) {
      newSSPvalues[i][j] = newSSPvalues[i][j];
      if ([-2, -3].includes(newSSPvalues[i][j])) {
        newSSPvalues[i][j] = -1;
      } else {
        if (newSSPvalues[i][j] > 0) {
          numDatasets += 1;
        }
        if (newSSPvalues[i][j] > pmax) {
          pmax = newSSPvalues[i][j];
        }
      }
    }
  }
  SSPvalues_secondary.forEach((p) => {
    if (p > 0) {
      numDatasets += 1;
    }
  });

  const suggestedThreshold =
    numDatasets > 0 ? -Math.log10(0.05 / numDatasets) : 1.3;

  const colorCutoffs = (
    [
      [-1, "rgb(105,105,105)"], // gray
      [0, "rgb(105,105,105)"],

      // Next trying to follow LD colors

      // dark blue to bright red HSV gradient in 5 steps:
      [0, "rgb(0, 0, 128)"], // dark blue (not significant)
      [suggestedThreshold, "rgb(0, 0, 128)"], // dark blue

      [suggestedThreshold, "rgb(0, 147, 142)"], // teal
      [suggestedThreshold + 1, "rgb(0, 147, 142)"],

      [suggestedThreshold + 1, "rgb(10, 166, 0)"], // green
      [suggestedThreshold + 2.5, "rgb(10, 166, 0)"],

      [suggestedThreshold + 2.5, "rgb(185, 167, 0)"], // yellow/orange
      [suggestedThreshold + 5, "rgb(185, 167, 0)"],

      [suggestedThreshold + 5, "rgb(204, 0, 24)"], // dark red
    ] as [number, string][]
  ).filter(([threshold]) => threshold >= -1 && threshold <= pmax); // remove colors that are out of bounds of data

  const normColorCutoffs = colorCutoffs.map(([threshold, color]) => [
    normalize(threshold, -1, pmax),
    color,
  ]);

  normColorCutoffs.push([1, colorCutoffs.slice(-1)[0][1]]);

  const tickParams = colorCutoffs.map(
    ([threshold]) => Math.round((threshold + Number.EPSILON) * 1e3) / 1e3
  );

  const data: Partial<PlotData>[] = [
    {
      z: newSSPvalues,
      x: genes.map((gene) => `${gene}`),
      y: tissues,
      //colorscale: 'Portland',
      type: "heatmap",
      name: "-log10(Simple Sum P-value)",
      hovertemplate:
        "Gene: %{x}" +
        "<br>Tissue: %{y}<br>" +
        `-log10(SS P-value): %{z: null}`,
      colorbar: {
        title: { text: "-log10(Simple<br>Sum P-value)" },
        // dtick0: 0,
        // dtick: 1,
        // autotick: false,
        tickmode: "array",
        tickvals: tickParams,
        ticktext: tickParams,
      },
      colorscale: normColorCutoffs as [number, string][],
    },
  ];

  const layout = {
    annotations: [],
    margin: {
      r: 50,
      t: 50,
      b: 125,
      l: 300,
    },
    autosize: false,
    width: imageWidth,
    height: imageHeight,
    font: { size: fontSize },
    xaxis: {
      //uncomment to set a max range of visible genes and set `dragmode=pan` to pan
      //range: [0, Math.min(genes.length, 10)],
      dtick: 1,
    },
  };

  return { heatmapLayout: layout, heatmapData: data };
};

interface SessionBoxProps {
  location: string;
  sessionId: string;
}

const SessionBox: React.FC<SessionBoxProps> = ({ sessionId, location }) => (
  <Grid
    component={Paper}
    maxWidth="400px"
    container
    direction="column"
    paddingBottom={1}
    sx={{
      ["& p, span"]: {
        padding: 2,
      },
    }}
  >
    <Grid sx={(theme) => ({ backgroundColor: theme.palette.grey[100] })}>
      <Typography>Session ID</Typography>
    </Grid>
    <Grid>
      <Link href={location}>
        <Typography>{sessionId}</Typography>
      </Link>
    </Grid>
    <Grid>
      <Typography variant="caption">
        Save the above string for your records to load or share your plot.
      </Typography>
    </Grid>
    <Grid>
      <Typography variant="caption">
        Plots older than 7 days are deleted.
      </Typography>
    </Grid>
  </Grid>
);
