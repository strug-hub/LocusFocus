"use client";

import { Typography } from "@mui/material";
import { useSearchParams } from "next/navigation";

const ColocResultsPage: React.FC = () => {
  const params = useSearchParams();

  return <Typography>Session ID: {params.get("sessionId")}!</Typography>;
};

export default ColocResultsPage;
