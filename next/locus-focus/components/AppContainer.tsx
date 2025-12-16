"use client";

import React from "react";
import { Container } from "@mui/material";
import { Header } from "@/components";

// interface LocusFocusDataContext {}

// export const LocusFocusContext = createContext<LocusFocusDataContext>({});

interface AppContainerProps {
  children: React.ReactNode;
}

const AppContainer: React.FC<AppContainerProps> = ({ children }) => {
  return (
    //    <LocusFocusContext.Provider value={{}}>
    <>
      <Header />
      <Container maxWidth="lg" sx={{ minHeight: "92vh", marginY: 3 }}>
        {children}
      </Container>
    </>
    //  </LocusFocusContext.Provider>
  );
};

export default AppContainer;
