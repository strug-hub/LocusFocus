"use client";

import { Suspense } from "react";
import { Coloc2ResultsPage } from "@/components";

const ColocResultsPage: React.FC = () => {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <Coloc2ResultsPage />
    </Suspense>
  );
};

export default ColocResultsPage;
