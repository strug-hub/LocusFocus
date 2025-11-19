"use client";

import {
  alpha,
  Box,
  Checkbox,
  Divider,
  FormControlLabel,
  Grid,
  MenuItem,
  Paper,
  TextField,
  Typography,
  useTheme,
} from "@mui/material";
import {
  FieldValues,
  Path,
  useForm,
  UseFormRegister,
  Validate,
} from "react-hook-form";
import { ColocFormFields } from "@/lib/ts/types";

export const validatePosition = (position: string) =>
  /^(chr)?(X|Y|[1-9]|1[0-9]|2[0-2])[:-]\d+-\d+$/.test(
    position.replaceAll(",", "")
  );

const LandingPage: React.FC = () => {
  const {
    register,
    watch,
    formState: { errors },
  } = useForm<ColocFormFields>({
    mode: "onChange", //validate on change
    defaultValues: {
      altCol: "ALT",
      betaCol: "",
      chromCol: "#CHROM",
      coloc2check: false,
      coordinate: "",
      GTExTissues: [],
      GTExVersion: "hg38",
      htmlFile: null,
      htmlFileCoordinate: "",
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
    },
  });

  const formValues = watch();

  const theme = useTheme();

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
              defaultValue="hg38"
              options={[
                { label: "hg19", value: "hg19" },
                { label: "hg38", value: "hg38" },
              ]}
              required={true}
              register={register}
              select
            />
          </Grid>
          <Grid
            size="grow"
            padding={3}
            sx={(theme) => ({
              backgroundColor: alpha(theme.palette.secondary.light, 0.4),
              borderRadius: 2,
            })}
          >
            <Typography>
              Select the human coordinate system your GWAS file is stored in.
              Note: LiftOver will be used to convert your GWAS data if necessary
              to match that of your selected GTEx dataset.
            </Typography>
          </Grid>
        </Grid>
      </Grid>
      <Divider />
      <Grid container direction="column">
        <Grid>
          <Typography variant="h5">Upload GWAS Data</Typography>
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
                label="Reference Allele Column NAme"
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
                defaultValue="quant"
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
                range: (value) => {
                  const range = value.replaceAll(",", "").split(":")[1];
                  const [p1, p2] = range.split("-");
                  if (p2 - p1 > 2e6) {
                    return "Range too large!";
                  }
                  if (p2 - p1 < 0) {
                    return "Range is less than 0!";
                  }
                  return true;
                },
              }}
              label="Coordinates (max 2Mbp)"
              hasError={!!errors.locus}
              errorText={errors.locus?.message}
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
    </Grid>
  );
};

export default LandingPage;

interface HFTextFieldProps<F extends FieldValues> {
  defaultValue?: string | number; //needed for select
  errorText?: string;
  hasError: boolean;
  name: Path<F>;
  label: string;
  options?: { label: string; value: string | number }[];
  placeholder?: string;
  required?: boolean;
  select?: boolean;
  register: UseFormRegister<any>;
  validate?: Validate<any, F> | Record<string, Validate<any, F>> | undefined;
}

function HFTextField<F extends FieldValues>({
  defaultValue,
  errorText,
  hasError,
  name,
  label,
  options,
  placeholder,
  register,
  required,
  select,
  validate,
}: HFTextFieldProps<F>) {
  return (
    <TextField
      sx={{ width: 200 }}
      defaultValue={defaultValue}
      error={hasError}
      helperText={errorText}
      label={label}
      placeholder={placeholder}
      required={true}
      select={select}
      title={label}
      //todo: maybe juse use UseController.... that will give you the default, at least
      {...register(name, {
        required: required ? `${name} is required!` : false,
        validate,
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
