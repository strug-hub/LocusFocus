"use client";

import React from "react";
import { Grid, Typography } from "@mui/material";

const Footer: React.FC = () => {
  return (
    <Grid
      container
      justifyContent="center"
      alignItems="center"
      flexGrow={1}
      sx={{
        top: "auto",
        bottom: 0,
        backgroundColor: (theme) => theme.palette.grey["700"],
        height: "8vh",
      }}
    >
      <Grid>
        <Typography
          sx={(theme) => ({
            color: theme.palette.getContrastText(theme.palette.grey["700"]),
          })}
          variant="h6"
        >
          Copyright {new Date().getFullYear()} LocusFocus
        </Typography>
      </Grid>
    </Grid>
  );
};

export default Footer;
