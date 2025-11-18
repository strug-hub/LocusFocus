import React from "react";
import Backdrop from "@mui/material/Backdrop";
import CircularProgress from "@mui/material/CircularProgress";
import { Grid, Typography } from "@mui/material";

const LoadingOverlay: React.FC<{ open: boolean; message?: string }> = ({
  message,
  open,
}) => {
  return (
    <Backdrop
      sx={{ color: "#fff", zIndex: (theme) => theme.zIndex.drawer + 1 }}
      open={open}
    >
      <Grid container direction="column" spacing={2} alignItems="center">
        <Grid>
          <CircularProgress color="inherit" />
        </Grid>
        {!!message && (
          <Grid>
            <Typography>{message}</Typography>
          </Grid>
        )}
      </Grid>
    </Backdrop>
  );
};

export default LoadingOverlay;
