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
