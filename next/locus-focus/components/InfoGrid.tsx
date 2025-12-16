import { alpha, Grid, GridProps } from "@mui/material";

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

export default InfoGrid;
