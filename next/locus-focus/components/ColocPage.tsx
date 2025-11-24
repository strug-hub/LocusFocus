"use client";
import { use, useEffect, useState } from "react";
import Link from "next/link";
import {
  alpha,
  Autocomplete,
  Box,
  Button,
  Checkbox,
  Divider,
  FormControlLabel,
  Grid,
  GridProps,
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
import { UploadButton } from "@/components";

const validatePosition = (position: string) =>
  /^(chr)?(X|Y|[1-9]|1[0-9]|2[0-2])[:-]\d+-\d+$/.test(
    position.replaceAll(",", "")
  );

const validateRange = (value: string) => {
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
  pvalCol: "",
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

const ColocPage: React.FC<{
  _tissuesV8: Promise<string[]>;
  _tissuesV10: Promise<string[]>;
}> = ({ _tissuesV8, _tissuesV10 }) => {
  "use no memo";
  const [genesLoading, setGenesLoading] = useState(false);
  const [genes, setGenes] = useState<string[]>([]);

  const v8Tissues = use(_tissuesV8);
  const v10Tissues = use(_tissuesV10);

  const {
    control,
    formState: { errors },
    register,
    watch,
  } = useForm<ColocFormFields>({
    mode: "onChange", //validate on change
    defaultValues,
  });

  const { field: gwasFileField } = useController({ name: "gwasFile", control });
  const { field: ldFileField } = useController({ name: "ldFile", control });
  const { field: htmlFileField } = useController({ name: "htmlFile", control });

  const formValues = watch();

  const theme = useTheme();

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
              {formValues.gwasFile
                ? formValues.gwasFile.name
                : "No file chosen"}
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
                name="pvalCol"
                required={formValues.coloc2check}
                label="P-value Column Name"
                hasError={!!errors.pvalCol}
                errorText={errors.pvalCol?.message}
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
                format: (value) =>
                  (!!value && validatePosition(value)) ||
                  "Invalid position format!",
                range: (value) => validateRange(value),
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
                format: (value) =>
                  (!!value && validatePosition(value)) ||
                  "Invalid position format!",
                range: (value) => validateRange(value),
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
              name="LDPopulations"
              register={register}
              defaultValue={defaultValues.LDPopulations}
              hasError={!!errors.LDPopulations}
              label="LD Population"
              errorText={errors.LDPopulations?.message}
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
                {formValues.ldFile ? formValues.ldFile.name : "No file chosen"}
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
            name="GTExTissues"
            hasError={!!errors.GTExTissues}
            label={`GTEx (V${formValues.GTExVersion === "V10" ? "10" : "8"}) Tissues`}
            errorText={errors.GTExTissues?.message}
            options={formValues.GTExVersion === "V8" ? v8Tissues : v10Tissues}
            selected={formValues.GTExTissues}
            wide
          />
        </Grid>
        <Grid>
          <Multiselect
            name="regionGenes"
            disabled={!formValues.locus}
            control={control}
            hasError={!!errors.regionGenes}
            label={`Select genes found in ${formValues.locus}`}
            errorText={errors.regionGenes?.message}
            options={genes}
            selected={formValues.regionGenes}
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
                  {formValues.ldFile
                    ? formValues.ldFile.name
                    : "No file chosen"}
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
          <Button variant="contained">Submit</Button>
        </Grid>
      </Grid>
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
      //todo: maybe juse use UseController.... that will give you the default, at least
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

const InfoGrid: React.FC<GridProps> = (props) => (
  <Grid
    {...props}
    size="grow"
    padding={3}
    sx={(theme) => ({
      backgroundColor: alpha(theme.palette.secondary.light, 0.4),
      borderRadius: 2,
      "p + p": {
        marginTop: 2,
      },
    })}
  />
);

interface MultiselectProps<F extends FieldValues>
  extends Omit<HFTextFieldProps<F>, "options" | "register" | "select"> {
  //TODO: clean this up
  control: Control<ColocFormFields, any, ColocFormFields>;
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
  const { field } = useController({ name, control });

  return (
    <Autocomplete
      sx={{ width: wide ? 350 : 200 }}
      disableCloseOnSelect
      options={options}
      multiple
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
          <li key={key} {...optionProps}>
            <Checkbox
              style={{ marginRight: 8 }}
              checked={selected.includes(option)}
              onChange={(_, checked) =>
                field.onChange(
                  checked
                    ? selected.concat(option)
                    : selected.filter((o) => o !== option)
                )
              }
            />
            {option}
          </li>
        );
      }}
      value={selected}
    />
  );
}
