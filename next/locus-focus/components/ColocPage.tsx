"use client";
import { useEffect, useMemo, useState } from "react";
import { useRouter } from "next/navigation";
import Link from "next/link";
import {
  Autocomplete,
  Box,
  Button,
  Checkbox,
  Divider,
  FormControlLabel,
  Grid,
  List,
  ListItem,
  MenuItem,
  Paper,
  SxProps,
  TextField,
  Typography,
  useTheme,
} from "@mui/material";
import {
  Control,
  FieldValues,
  Path,
  useController,
  useForm,
  UseFormRegister,
  Validate,
} from "react-hook-form";
import { ColocFormFields } from "@/lib/ts/types";
import { InfoGrid, LoadingOverlay, Modal, UploadButton } from "@/components";

interface JobStatusResponse {
  error_message?: string;
  error_title?: string;
  redirect_url?: string;
  status: "PENDING" | "RUNNING" | "FAILURE" | "SUCCESS";
  stage_index?: number;
  stage_count?: number;
}
// Validate position, return true if empty string
const validatePosition = (position: string) => {
  if (
    !position ||
    /^(chr)?(X|Y|[1-9]|1[0-9]|2[0-2])[:-]\d+-\d+$/.test(
      position.replaceAll(",", "")
    )
  ) {
    return true;
  } else {
    return "Invalid position format!";
  }
};

// Validate range, assumes either empty string or valid position format,
// in other words, that it's already been validated by validatePosition
const validateRange = (value: string) => {
  if (!value) {
    return true;
  }
  const range = value.replaceAll(",", "").split(":")[1];
  const [p1, p2] = range.split("-");
  if (+p2 - +p1 > 2e6) {
    return "Range too large!";
  }
  if (+p2 - +p1 < 0) {
    return "Range is less than 0!";
  }
  return true;
};

const lDPopsHg19 = [
  { label: "None", value: "" },
  { label: "1000 Genomes 2012 EUR", value: "EUR" },
  { label: "1000 Genomes 2012 AFR", value: "AFR" },
  { label: "1000 Genomes 2012 AMR", value: "AMR" },
  { label: "1000 Genomes 2012 ASN", value: "ASN" },
];

const lDPopsHg38 = [
  { label: "None", value: "" },
  { label: "1000 Genomes 2018 EUR", value: "EUR" },
  { label: "1000 Genomes 2018 AFR", value: "AFR" },
  { label: "1000 Genomes 2018 AMR", value: "AMR" },
  { label: "1000 Genomes 2018 ASN", value: "ASN" },
  { label: "1000 Genomes 2018 SAS", value: "SAS" },
  { label: "1000 Genomes 2018 NFE", value: "NFE" },
];

const submitFieldMapping: Record<string, string> = {
  altCol: "alt-col",
  betaCol: "beta-col",
  chromCol: "chrom-col",
  coloc2check: "coloc2check",
  coordinate: "coordinate",
  GTExTissues: "GTEx-tissues",
  GTExVersion: "GTEx-version",
  htmlFile: "html-file",
  gwasFile: "gwas-file",
  htmlFileCoordinate: "html-file-coordinate",
  LDPopulations: "LD-populations",
  leadSnp: "lead-snp",
  locus: "locus",
  ldFile: "ld-file",
  mafCol: "maf-col",
  markerCheckbox: "markercheckbox",
  multiRegion: "multi-region",
  numcases: "numcases",
  numsamplesCol: "numsamples-col",
  posCol: "pos-col",
  pvalCol: "pval-col",
  refCol: "ref-col",
  regionGenes: "region-genes",
  snpCol: "snp-col",
  separateTestCheckbox: "separate-test-checkbox",
  sessionId: "session-id",
  setbasedP: "setbasedP",
  SSlocus: "SSlocus",
  stderrCol: "stderr-col",
  studytype: "studytype",
};

const coloc2Cols = [
  "betaCol",
  "stderrCol",
  "numsamplesCol",
  "mafCol",
  "studytype",
];

const defaultValues = {
  altCol: "ALT",
  betaCol: "",
  chromCol: "#CHROM",
  coloc2check: false,
  coordinate: "hg38",
  GTExTissues: [] as string[],
  GTExVersion: "V10",
  htmlFile: null,
  htmlFileCoordinate: "gwas",
  gwasFile: null,
  LDPopulations: "",
  leadSnp: "",
  locus: "1:205,500,000-206,000,000",
  ldFile: null,
  mafCol: "",
  markerCheckbox: false,
  multiRegion: "",
  numcases: 0,
  numsamplesCol: "",
  posCol: "POS",
  pvalCol: "P",
  refCol: "REF",
  regionGenes: [],
  snpCol: "SNP",
  separateTestCheckbox: "",
  sessionId: "",
  setbasedP: "",
  SSlocus: "",
  stderrCol: "",
  studytype: "quant",
};

const fetchTissues = async () => {
  const v8 = await fetch(
    `${process.env.NEXT_PUBLIC_BROWSER_API_HOST}/gtex/v8/tissues_list`
  );

  const v8Json = await v8.json();

  const v10 = await fetch(
    `${process.env.NEXT_PUBLIC_BROWSER_API_HOST}/gtex/v10/tissues_list`
  );

  const v10Json = await v10.json();

  const vals = await Promise.all([v8Json, v10Json]);

  return { v8: vals[0], v10: vals[1] };
};

const ColocPage: React.FC = () => {
  const [genesLoading, setGenesLoading] = useState(false);
  const [genes, setGenes] = useState<string[]>([]);
  const [loading, setLoading] = useState(false);
  const [errorModalText, setErrorModalText] = useState("");
  const [runningModalOpen, setRunningModalOpen] = useState<{
    stageIndex?: number;
    stageCount?: number;
  } | null>(null);
  const [pendingModalOpen, setPendingModalOpen] = useState(false);
  const [v8Tissues, setv8Tissues] = useState<string[]>([]);
  const [v10Tissues, setv10Tissues] = useState<string[]>([]);

  const router = useRouter();

  const {
    control,
    formState: { errors, isValid },
    register,
    trigger,
    watch,
  } = useForm<ColocFormFields>({
    mode: "onChange", //validate on change
    defaultValues,
  });

  const { field: gwasFileField } = useController({
    name: "gwasFile",
    control,
    rules: { required: "GWAS file is required!" },
  });
  const { field: ldFileField } = useController({ name: "ldFile", control });
  const { field: htmlFileField } = useController({ name: "htmlFile", control });
  const { field: GTExTissuesField } = useController({
    name: "GTExTissues",
    control,
  });
  const { field: regionGenesField } = useController({
    name: "regionGenes",
    control,
  });

  const formValues = watch();

  const theme = useTheme();

  // run initial validation
  useEffect(() => {
    if (trigger) {
      trigger();
    }
  }, [trigger]);

  useEffect(() => {
    const _fetchTissues = async () => {
      setLoading(true);
      const { v8, v10 } = await fetchTissues();
      setv8Tissues(v8);
      setv10Tissues(v10);
      setLoading(false);
    };

    _fetchTissues();
  }, []);

  //reset tissues and genes when gtex version changes
  useEffect(() => {
    GTExTissuesField.onChange([]);
    regionGenesField.onChange([]);
  }, [formValues.GTExVersion]);

  useEffect(() => {
    if (formValues.locus && !errors.locus) {
      const [chr, range] = formValues.locus.split(":");
      const [start, end] = range.replaceAll(",", "").split("-");
      setGenesLoading(true);
      fetch(
        `${process.env.NEXT_PUBLIC_BROWSER_API_HOST}/genenames/${formValues.coordinate}/${chr}/${start}/${end}`,
        { mode: "cors", headers: { "content-type": "application/json" } }
      )
        .then((r) => r.json())
        .then((r) => setGenes(r as unknown as string[]))
        .finally(() => setGenesLoading(false));
    }
  }, [formValues.locus, formValues.coordinate, errors.locus]);

  const resetModals = () => {
    setLoading(false);
    setPendingModalOpen(false);
    setRunningModalOpen(null);
    setErrorModalText("");
  };

  const handleJobStatus = async (jobStatusURL: string, sessionId: string) => {
    try {
      const statusResponse = await fetch(jobStatusURL);
      const statusData: JobStatusResponse = await statusResponse.json();
      const jobStatus = statusData.status;
      if (jobStatus == "PENDING") {
        setPendingModalOpen(true);
        await new Promise((resolve) =>
          setTimeout(
            () => resolve(handleJobStatus(jobStatusURL, sessionId)),
            10000
          )
        );
      }

      if (jobStatus == "RUNNING") {
        setPendingModalOpen(false);
        const stageIndex = statusData.stage_index
          ? statusData.stage_index - 1
          : undefined;
        const stageCount = statusData.stage_count;
        setRunningModalOpen({ stageIndex, stageCount });
        await new Promise((resolve) =>
          setTimeout(
            () => resolve(handleJobStatus(jobStatusURL, sessionId)),
            500
          )
        );
      }

      if (jobStatus == "FAILURE") {
        setPendingModalOpen(false);
        setRunningModalOpen(null);
        setErrorModalText(
          `Error: ${statusData.error_title}; Details: ${statusData.error_message}`
        );
      }

      if (jobStatus == "SUCCESS") {
        router.push(`colocResults?sessionId=${sessionId}`);
      }
    } catch (error) {
      resetModals();
      setErrorModalText(
        "There was an error fetching the job status. Please try again or contact a system administrator."
      );
      console.error("Error fetching job status:", error);
    }
  };

  return (
    <Grid container direction="column" spacing={3} margin={2}>
      <Grid container>
        <Box width="100%" marginY={6} textAlign="center">
          <Paper elevation={3}>
            <Grid padding={5} spacing={2} container direction="column">
              <Grid>
                <Typography
                  fontWeight={theme.typography.fontWeightLight}
                  variant="h3"
                  color="primary"
                >
                  LocusFocus
                </Typography>
              </Grid>
              <Grid>
                <Typography
                  variant="h6"
                  fontWeight={theme.typography.fontWeightLight}
                >
                  Colocalization Testing Across Datasets
                </Typography>
              </Grid>
              <Grid>
                <Typography
                  variant="h6"
                  fontWeight={theme.typography.fontWeightBold}
                >
                  Version release: 1.6.0 alpha
                </Typography>
              </Grid>
            </Grid>
          </Paper>
        </Box>
      </Grid>
      <Grid container direction="column">
        <Grid>
          <Typography variant="h5">Select a GWAS Coordinate Sytem</Typography>
        </Grid>
        <Grid container spacing={2} direction="row" wrap="nowrap">
          <Grid>
            <HFTextField
              errorText={errors.coordinate?.message}
              hasError={!!errors.coordinate}
              label="Select GWAS Coordinate System"
              name="coordinate"
              defaultValue={defaultValues.coordinate}
              options={[
                { label: "hg19", value: "hg19" },
                { label: "hg38", value: "hg38" },
              ]}
              required={true}
              register={register}
              select
            />
          </Grid>
          <InfoGrid>
            <Typography>
              Select the human coordinate system your GWAS file is stored in.
              Note: LiftOver will be used to convert your GWAS data if necessary
              to match that of your selected GTEx dataset.
            </Typography>
          </InfoGrid>
        </Grid>
      </Grid>
      <Divider />
      <Grid container direction="column">
        <Grid>
          <Typography variant="h5">Upload GWAS Data</Typography>
        </Grid>
        <Grid container alignItems="center" direction="row" spacing={2}>
          <Grid>
            <UploadButton
              onChange={(e) =>
                gwasFileField.onChange((e.currentTarget?.files || [null])[0])
              }
              multiple={false}
              label="Choose file"
              accept=".txt,.tsv"
            />
          </Grid>
          <Grid>
            <Typography variant="caption">
              {formValues.gwasFile ? (
                formValues.gwasFile.name
              ) : (
                <Typography color="error">No file chosen</Typography>
              )}
            </Typography>
          </Grid>
        </Grid>
        <Grid>
          <HFTextField
            register={register}
            name="snpCol"
            required={true}
            label="Marker Column"
            hasError={!!errors.snpCol}
            errorText={errors.snpCol?.message}
          />
        </Grid>
        <Grid>
          <FormControlLabel
            label="Use marker ID column to infer variant position and alleles"
            control={
              <Checkbox
                {...register("markerCheckbox", {
                  deps: ["chromCol", "posCol", "refCol", "altCol"],
                })}
              />
            }
          />
        </Grid>
        {!formValues.markerCheckbox && (
          <Grid container direction="row">
            <Grid>
              <HFTextField
                register={register}
                name="chromCol"
                required={!formValues.markerCheckbox}
                label="Chromomsome Column Name"
                hasError={!!errors.chromCol}
                errorText={errors.chromCol?.message}
              />
            </Grid>
            <Grid>
              <HFTextField
                register={register}
                name="posCol"
                required={!formValues.markerCheckbox}
                label="Position Column Name"
                hasError={!!errors.posCol}
                errorText={errors.posCol?.message}
              />
            </Grid>
            <Grid>
              <HFTextField
                register={register}
                name="refCol"
                required={!formValues.markerCheckbox}
                label="Reference Allele Column Name"
                hasError={!!errors.refCol}
                errorText={errors.refCol?.message}
              />
            </Grid>
            <Grid>
              <HFTextField
                register={register}
                name="altCol"
                required={!formValues.markerCheckbox}
                label="Alternate Allele Column Name"
                hasError={!!errors.altCol}
                errorText={errors.altCol?.message}
              />
            </Grid>
          </Grid>
        )}
        <Grid>
          <FormControlLabel
            label="Add COLOC2 method (additional fields required)"
            control={
              <Checkbox
                {...register("coloc2check", {
                  deps: [
                    "betaCol",
                    "stderrCol",
                    "numsamplesCol",
                    "pvalCol",
                    "mafCol",
                    "studytype",
                  ],
                })}
              />
            }
          />
        </Grid>
        <Grid>
          <HFTextField
            register={register}
            name="pvalCol"
            required={!formValues.pvalCol}
            label="P-value Column Name"
            hasError={!!errors.pvalCol}
            errorText={errors.pvalCol?.message}
          />
        </Grid>
        <Divider />
        {formValues.coloc2check && (
          <Grid container direction="row">
            <Grid>
              <HFTextField
                register={register}
                name="betaCol"
                required={formValues.coloc2check}
                label="Beta Column Name"
                hasError={!!errors.betaCol}
                errorText={errors.betaCol?.message}
              />
            </Grid>
            <Grid>
              <HFTextField
                register={register}
                name="stderrCol"
                required={formValues.coloc2check}
                label="Standard Error Column Name"
                hasError={!!errors.stderrCol}
                errorText={errors.stderrCol?.message}
              />
            </Grid>
            <Grid>
              <HFTextField
                register={register}
                name="numsamplesCol"
                required={formValues.coloc2check}
                label="Number of Samples Column Name"
                hasError={!!errors.numsamplesCol}
                errorText={errors.numsamplesCol?.message}
              />
            </Grid>

            <Grid>
              <HFTextField
                register={register}
                name="mafCol"
                required={formValues.coloc2check}
                label="MAF Column Name"
                hasError={!!errors.mafCol}
                errorText={errors.mafCol?.message}
              />
            </Grid>
            <Grid>
              <HFTextField
                register={register}
                name="studytype"
                required={formValues.coloc2check}
                label="Study Type"
                select
                defaultValue={defaultValues.studytype}
                options={[
                  { label: "Quantative", value: "quant" },
                  { label: "Case-Control", value: "cc" },
                ]}
                hasError={!!errors.studytype}
                errorText={errors.studytype?.message}
              />
            </Grid>
            <Grid>
              <Divider />
            </Grid>
          </Grid>
        )}
        <Grid container direction="row">
          <Grid>
            <HFTextField
              register={register}
              name="locus"
              required
              validate={{
                format: validatePosition,
                range: validateRange,
              }}
              label="Coordinates (max 2Mbp)"
              hasError={!!errors.locus}
              errorText={errors.locus?.message}
              wide
            />
          </Grid>
          <Grid>
            <HFTextField
              hasError={!!errors.leadSnp}
              register={register}
              name="leadSnp"
              placeholder="default: top marker"
              label="Lead Marker name"
              errorText={errors.leadSnp?.message}
            />
          </Grid>
        </Grid>
      </Grid>
      <Grid>
        <Divider />
      </Grid>
      <Grid container direction="column">
        <Grid>
          <Typography variant="h5">
            Select Simple Sum Colocalization Region
          </Typography>
        </Grid>
        <Grid>
          <Typography variant="caption">
            Leave blank for default (+/- 0.1Mbp from the lead SNP)
          </Typography>
        </Grid>
        <Grid container direction="row">
          <Grid>
            <HFTextField
              hasError={!!errors.SSlocus}
              register={register}
              name="SSlocus"
              placeholder="chr:start-end"
              label="Coordinates"
              errorText={errors.SSlocus?.message}
              validate={{
                format: validatePosition,
                range: validateRange,
              }}
              wide
            />
          </Grid>
          <InfoGrid>
            <Typography>
              By convention, colocalization methods default to testing
              colocalization +/- 0.1 Mbp of the top variant in the primary (e.g.
              GWAS) dataset. You can override this behaviour by entering the
              subregion to test colocalization on.
            </Typography>
          </InfoGrid>
        </Grid>
      </Grid>
      <Grid>
        <Divider />
      </Grid>
      <Grid container direction="column">
        <Grid>
          <Typography variant="h5">Linkage Disequilibrium</Typography>
        </Grid>
        <Grid>
          <Typography variant="caption">
            Select an LD population to use from the 1000 Genomes dataset OR
            upload an LD matrix file.
          </Typography>
        </Grid>
        <Grid container direction="row" spacing={10}>
          <Grid>
            <HFTextField
              disabled={!!formValues.ldFile}
              select
              deps={["ldFile"]}
              name="LDPopulations"
              register={register}
              defaultValue={defaultValues.LDPopulations}
              hasError={!!errors.LDPopulations}
              label="LD Population"
              errorText={errors.LDPopulations?.message}
              validate={{
                required: (value) => {
                  return (
                    !!value ||
                    !!formValues.ldFile ||
                    "LD population is required"
                  );
                },
              }}
              options={
                formValues.coordinate === "hg19" ? lDPopsHg19 : lDPopsHg38
              }
              wide
            />
          </Grid>
          <Grid container spacing={1} direction="row" alignItems="center">
            <Grid>
              <UploadButton
                disabled={!!formValues.LDPopulations}
                onChange={(e) =>
                  ldFileField.onChange((e.currentTarget?.files || [null])[0])
                }
                multiple={false}
                label="Choose file"
                accept=".ld"
                tooltipMessage=".ld (optional): LD matrix file generated with PLINK or similar"
              />
            </Grid>
            <Grid>
              <Typography variant="caption">
                {formValues.ldFile ? (
                  formValues.ldFile.name
                ) : (
                  <Typography
                    color={!!formValues.LDPopulations ? undefined : "error"}
                  >
                    No file chosen
                  </Typography>
                )}
              </Typography>
            </Grid>
          </Grid>
        </Grid>
        <InfoGrid>
          <Typography>
            For the most accurate results, the LD (r2) matrix for your primary
            (e.g. GWAS) dataset is recommended to be uploaded as a .ld square
            matrix file. If the .ld file is unavailable, you may choose one of
            the publicly-available 1000 Genomes population datasets.
          </Typography>
          <List>
            <ListItem>EUR: European </ListItem>
            <ListItem>NFE: Non-Finnish European </ListItem>
            <ListItem>AFR: African EAS: East Asian</ListItem>
            <ListItem>SAS: South Asian </ListItem>
            <ListItem>AMR: Ad Mixed American</ListItem>
          </List>
        </InfoGrid>
      </Grid>
      <Grid>
        <Divider />
      </Grid>
      <Grid>
        <Typography variant="h5">Select Secondary Datasets</Typography>
        <Typography variant="caption">
          At least one secondary dataset must be uploaded to perform
          colocalization analysis.
        </Typography>
      </Grid>
      <Grid>
        <Typography variant="h6">Select GTEx eQTL Data to Render</Typography>
        <Typography variant="caption">
          Please allow sufficient time for analysis if selecting many tissues
          and genes. <br /> Analyses may take &gt; 30 minutes when selecting
          many tissue-gene pairs.
        </Typography>
      </Grid>
      <Grid>
        <HFTextField
          select
          name="GTExVersion"
          register={register}
          defaultValue={defaultValues.GTExVersion}
          hasError={!!errors.GTExVersion}
          label="GTEx Version"
          errorText={errors.GTExVersion?.message}
          options={[
            { label: "GTEx V8 (hg38)", value: "V8" },
            { label: "GTEx V10 (hg38)", value: "V10" },
          ]}
          wide
        />
      </Grid>
      <Grid container direction="row">
        <Grid>
          <Multiselect
            control={control}
            required
            name="GTExTissues"
            hasError={!!errors.GTExTissues}
            label={`GTEx (V${formValues.GTExVersion === "V10" ? "10" : "8"}) Tissues`}
            errorText={errors.GTExTissues?.message}
            options={formValues.GTExVersion === "V8" ? v8Tissues : v10Tissues}
            selected={formValues.GTExTissues}
            validate={{
              required: (val) =>
                (!!val && !!val.length) || "GTEx Tissues is required!",
            }}
            wide
          />
        </Grid>
        <Grid>
          <Multiselect
            required
            name="regionGenes"
            disabled={!formValues.locus || genesLoading}
            control={control}
            hasError={!!errors.regionGenes}
            label={`Select genes found in ${formValues.locus}`}
            errorText={errors.regionGenes?.message}
            options={genes}
            selected={formValues.regionGenes}
            validate={{
              required: (val) =>
                (!!val && !!val.length) || "Region Genes is required!",
            }}
            wide
          />
        </Grid>
      </Grid>
      <Grid container direction="column">
        <Grid>
          <Typography variant="h6">Upload Secondary Datasets</Typography>
        </Grid>
        <Grid container spacing={1} direction="row">
          <Grid container direction="column" spacing={2}>
            <Grid container alignItems="center" direction="row">
              <Grid>
                <UploadButton
                  onChange={(e) =>
                    htmlFileField.onChange(
                      (e.currentTarget?.files || [null])[0]
                    )
                  }
                  multiple={false}
                  label="Choose file"
                  accept=".html"
                  tooltipMessage="HTML file generated with LocusFocus"
                />
              </Grid>
              <Grid>
                <Typography variant="caption">
                  {formValues.ldFile ? (
                    formValues.ldFile.name
                  ) : (
                    <Typography>No file chosen</Typography>
                  )}
                </Typography>
              </Grid>
            </Grid>
            <Grid>
              <HFTextField
                defaultValue={defaultValues.htmlFileCoordinate}
                select
                name="htmlFileCoordinate"
                register={register}
                hasError={!!errors.htmlFileCoordinate}
                label="Select Coordinate System"
                errorText={errors.htmlFileCoordinate?.message}
                options={[
                  { label: "Same as GWAS", value: "gwas" },
                  { label: "hg19", value: "hg19" },
                  { label: "hg38", value: "hg38" },
                ]}
                wide
              />
            </Grid>
          </Grid>
          <InfoGrid>
            <Typography>
              For convenience, we are making available eQTL datasets from the
              GTEx project for use as secondary datasets.
            </Typography>
            <Typography>
              LocusFocus, however, can perform more custom colocalization
              analyses utilizing other secondary (e.g. eQTL, mQTL, other
              phenotypes) datasets provided by the user.
            </Typography>
            <Typography>
              Custom secondary datasets may be uploaded after conversion to HTML
              format as described above or in the documentation.
            </Typography>
            <Typography>
              Refer to the documentation on how to generate the HTML file. You
              may use the merge_and_convert_to_html.py script, or
              merge_and_convert_to_html_coloc2.py. You may use provided sample
              datasets as a guide to formatting your files.
            </Typography>
            <List sx={{ listStyleType: "disc", marginLeft: 3 }}>
              <ListItem sx={{ display: "list-item" }}>
                <Typography>
                  Refer to the{" "}
                  <Link
                    target="_blank"
                    href="https://locusfocus.readthedocs.io/en/latest/quick_start.html#formatting-custom-secondary-datasets"
                  >
                    documentation
                  </Link>{" "}
                  on how to generate the HTML file.
                </Typography>
              </ListItem>
              <ListItem sx={{ display: "list-item" }}>
                <Typography>
                  You may use the{" "}
                  <Link
                    target="_blank"
                    href="https://github.com/strug-hub/LocusFocus/blob/master/merge_and_convert_to_html.py"
                  >
                    merge_and_convert_to_html.py
                  </Link>{" "}
                  script, or{" "}
                  <Link
                    target="_blank"
                    href="https://github.com/strug-hub/LocusFocus/blob/master/merge_and_convert_to_html_coloc2.py"
                  >
                    merge_and_convert_to_html_coloc2.py
                  </Link>
                  .
                </Typography>
              </ListItem>
              <ListItem sx={{ display: "list-item" }}>
                <Typography>
                  You may use provided{" "}
                  <Link
                    target="_blank"
                    href="https://github.com/strug-hub/LocusFocus/tree/master/data/sample_datasets"
                  >
                    sample datasets
                  </Link>{" "}
                  as a guide to formatting your files.
                </Typography>
              </ListItem>
            </List>
          </InfoGrid>
        </Grid>
      </Grid>
      <Grid>
        <Divider />
      </Grid>
      <Grid container direction="column">
        <Grid>
          <Typography variant="h5">
            Stage one set-based p-value threshold
          </Typography>
        </Grid>
        <Grid>
          <HFTextField
            register={register}
            hasError={!!errors.setbasedP}
            errorText={errors.setbasedP?.message}
            name="setbasedP"
            label="Threshold"
            placeholder="default: 0.05 / (number of tissues &times; number of genes + additional secondary datasets uploaded)"
            sx={{ width: "100%" }}
          />
        </Grid>
        <Grid>
          <InfoGrid>
            <Typography>
              For the Simple Sum method, a first-stage set-based Bonferroni
              p-value threshold is used for the set of secondary datasets with
              alpha 0.05 (0.05 divided by the number of secondary datasets)
            </Typography>
            <Typography>
              Enter a value if you would like to override the default threshold.
            </Typography>
          </InfoGrid>
        </Grid>
      </Grid>
      <Grid>
        <Divider />
      </Grid>
      <Grid container direction="column">
        <Grid>
          <Button
            disabled={!isValid}
            onClick={async () => {
              //TODO: moveout
              setLoading(true);

              const formData = new FormData();

              Object.entries(formValues)
                // drop coloc2 fields if not selected
                .filter(
                  ([k]) => !formValues.coloc2check && !coloc2Cols.includes(k)
                )
                // map booleans
                .map(([k, v]) => [
                  k,
                  typeof v === "boolean" ? (v === true ? 1 : "") : v,
                ])
                // append to form data
                .forEach(([k, v]) => {
                  // handle arrays
                  if (Array.isArray(v)) {
                    v.forEach((item) => {
                      formData.append(submitFieldMapping[k], item);
                    });
                  } else {
                    formData.append(submitFieldMapping[k], v);
                  }
                });

              try {
                const response = await fetch(
                  `${process.env.NEXT_PUBLIC_BROWSER_API_HOST}`,
                  {
                    method: "POST",
                    body: formData,
                  }
                );

                if (!response.status.toString().startsWith("2")) {
                  throw response;
                }
                const content = await response.json();

                if (!content.queued) {
                  router.push(`/colocResults?sessionId=${content.sessionId}}`);
                } else {
                  handleJobStatus(
                    `${process.env.NEXT_PUBLIC_BROWSER_API_HOST}/job/status/${content.session_id}`,
                    content.session_id
                  );
                }
              } catch (e: unknown) {
                console.error(e);
                if ((e as Response)?.status.toString().startsWith("4")) {
                  const error = await (e as Response).json();
                  const message = error.message;
                  setErrorModalText(message);
                } else {
                  setErrorModalText(
                    "The job failed due to an unexpected error, please try again later."
                  );
                }
              } finally {
                setLoading(false);
              }
            }}
            variant="contained"
          >
            Submit
          </Button>
        </Grid>
        <Grid>
          {Object.entries(errors).map(([k, v]) => (
            <Typography key={k} color="error">
              <>{v?.message}</>
            </Typography>
          ))}
        </Grid>
      </Grid>
      <LoadingOverlay open={loading} />
      <Modal open={!!errorModalText} onClose={() => setErrorModalText("")}>
        <Typography color="error">{errorModalText}</Typography>
      </Modal>
      <Modal
        open={pendingModalOpen}
        sx={{ alignItems: "center" }}
        onClose={() => setPendingModalOpen(false)}
      >
        <Typography variant="h4" color="info">
          Job pending....
        </Typography>
      </Modal>
      {!!runningModalOpen && (
        <Modal
          open={!!runningModalOpen}
          onClose={() => setRunningModalOpen(null)}
        >
          <Typography variant="h4" color="info">
            Job running, stage {runningModalOpen?.stageIndex} of{" "}
            {runningModalOpen?.stageCount}...
          </Typography>
        </Modal>
      )}
    </Grid>
  );
};

export default ColocPage;

interface HFTextFieldProps<F extends FieldValues> {
  defaultValue?: string | number; //needed for select
  deps?: string[];
  disabled?: boolean;
  errorText?: string;
  hasError: boolean;
  name: Path<F>;
  label: string;
  options?: { label: string; value: string | number }[];
  placeholder?: string;
  required?: boolean;
  select?: boolean;
  sx?: SxProps;
  register: UseFormRegister<any>;
  validate?: Validate<any, F> | Record<string, Validate<any, F>> | undefined;
  wide?: boolean;
}

function HFTextField<F extends FieldValues>({
  defaultValue,
  deps,
  disabled,
  errorText,
  hasError,
  name,
  label,
  options,
  placeholder,
  register,
  required,
  select,
  sx,
  validate,
  wide,
}: HFTextFieldProps<F>) {
  return (
    <TextField
      sx={{ width: wide ? 350 : 200, ...sx }}
      defaultValue={defaultValue}
      disabled={disabled}
      error={hasError}
      helperText={errorText}
      label={label}
      placeholder={placeholder}
      required={required}
      select={select}
      title={label}
      {...register(name, {
        required: required ? `${name} is required!` : false,
        validate,
        deps,
      })}
    >
      {options?.map((o) => (
        <MenuItem key={o.value} value={o.value}>
          {o.label}
        </MenuItem>
      ))}
    </TextField>
  );
}

interface MultiselectProps<F extends FieldValues> extends Omit<
  HFTextFieldProps<F>,
  "options" | "register" | "select"
> {
  control: Control<F, any, F>;
  options: string[];
  selected: string[];
}

export function Multiselect<F extends FieldValues>({
  control,
  defaultValue,
  disabled,
  errorText,
  hasError,
  name,
  label,
  options,
  placeholder,
  required,
  selected,
  validate,
  wide,
}: MultiselectProps<F>) {
  const { field } = useController({ name, control, rules: { validate } });

  const allOptions = useMemo(() => {
    return ["Select All"].concat(options);
  }, [options]);

  const handleChange = (value: string) => {
    const optionIsSeleted =
      selected.includes(value) ||
      (value === "Select All" && selected.length === options.length);

    if (optionIsSeleted && value === "Select All") {
      field.onChange([]);
    } else if (!optionIsSeleted && value === "Select All") {
      field.onChange(options);
    } else if (optionIsSeleted) {
      field.onChange(selected.filter((v) => v !== value));
    } else {
      field.onChange(selected.concat(value));
    }
  };

  return (
    <Autocomplete
      sx={{ width: wide ? 350 : 200 }}
      disableCloseOnSelect
      options={allOptions}
      multiple
      onChange={(_, __, reason) => {
        if (reason === "clear") {
          field.onChange([]);
        }
      }}
      renderValue={() => `${selected.length} selected`}
      renderInput={(params) => (
        <TextField
          {...params}
          disabled={disabled}
          defaultValue={defaultValue}
          label={label}
          error={hasError}
          helperText={errorText}
          required={required}
          placeholder={placeholder}
        />
      )}
      renderOption={(props, option) => {
        const { key, ...optionProps } = props;
        return (
          <ListItem key={key} {...optionProps}>
            <Checkbox
              style={{ marginRight: 8 }}
              checked={
                selected.includes(option) ||
                (option === "Select All" &&
                  !!selected.length &&
                  selected.length === options.length)
              }
              onChange={() => handleChange(option)}
            />
            <span onClick={() => handleChange(option)}>{option}</span>
          </ListItem>
        );
      }}
      value={selected}
    />
  );
}
