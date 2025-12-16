export interface ColocFormFields extends Record<string, any> {
  altCol: string; // alt-col
  betaCol: string; // beta-col
  chromCol: string; // chrom-col
  coloc2check: boolean;
  coordinate: string;
  GTExTissues: string[]; // GTEx-tissues
  GTExVersion: string; // GTEx-version
  htmlFile: File | null; // html-file
  gwasFile: File | null; // gwas-file
  htmlFileCoordinate: string; // html-file-coordinate
  LDPopulations: string; // LD-populations
  leadSnp: string; // lead-snp
  locus: string;
  ldFile: File | null; // ld-file
  mafCol: string; // maf-col
  markerCheckbox: boolean;
  multiRegion: string; // multi-region
  numcases: number;
  numsamplesCol: string; // numsamples-col
  posCol: string; // pos-col
  pvalCol: string; // pval-col
  refCol: string; // ref-col
  regionGenes: string[]; // region-genes
  snpCol: string; // snp-col
  separateTestCheckbox: string; //separate-test-checkbox
  sessionId: string; //session-id
  setbasedP: string;
  SSlocus: string;
  stderrCol: string; //stderr-col
  studytype: string;
}

export interface Coloc2AnalysisResults {
  SSPvalues_file: SSPvaluesFile;
  coloc2_file: Coloc2File;
  genesfile: GeneRecord[];
  metadata_file: MetadataFile;
  sessionfile: SessionFile;
}

export interface SSPvaluesFile {
  Genes: string[];
  Tissues: string[];
  Secondary_dataset_titles: string[];
  SSPvalues: number[][];
  Num_SNPs_Used_for_SS: number[][];
  Computation_method: string[][];
  SSPvalues_secondary: number[];
  Num_SNPs_Used_for_SS_secondary: number[][];
  Computation_method_secondary: string[][];
  First_stages: string[];
  First_stage_Pvalues: string[];
}

export interface Coloc2File {
  ProbeID: string[];
  PPH4abf: number[];
}

export interface GeneRecord {
  name: string;
  txStart: number;
  txEnd: number;
  exonStarts: number[];
  exonEnds: number[];
}

export interface MetadataFile {
  datetime: string; // ISO timestamp
  files_uploaded: string[];
  session_id: string;
  type: string;
}

// tissue types are keys here as well and have arrays of values
export interface SessionFile {
  chrom: number;
  coordinate: string;
  dataset_title?: string;
  endbp: number;
  first_stages: string[];
  first_stage_Pvalues: string[];
  gene: string;
  gtex_genes: string[];
  gtex_tissues: string[];
  gtex_version: string;
  inferVariant: boolean;
  lead_snp: string;
  ld_populations: string;
  ld_values: number[];
  multiple_tests?: boolean;
  numGTExMatches: number;
  num_SS_snps: number;
  positions: number[];
  pvalues: number[];
  regions?: string[];
  runcoloc2: boolean;
  secondary_dataset_colnames: string[];
  secondary_dataset_titles: string[];
  sessionid: string;
  set_based_p: string;
  snps: string[];
  snps_used_in_test?: string[];
  SS_region: number[];
  startbp: number;
  std_snp_list: string[];
  success: boolean;
  snp_warning: boolean;
  thresh: number;
}
