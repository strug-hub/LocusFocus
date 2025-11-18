"use client";

import React from "react";
import { Box, Container } from "@mui/material";
import { Header } from "@/components";

// interface LocusFocusDataContext {}

// export const LocusFocusContext = createContext<LocusFocusDataContext>({});

interface AppContainerProps {
  children: React.ReactNode;
}

const AppContainer: React.FC<AppContainerProps> = ({ children }) => {
  return (
    //    <LocusFocusContext.Provider value={{}}>
    <Container maxWidth={false} sx={{ minHeight: "92vh" }}>
      <Header />
      <Box flexGrow={1} overflow="auto" padding={2}>
        {children}
      </Box>
    </Container>
    //  </LocusFocusContext.Provider>
  );
};

export default AppContainer;
