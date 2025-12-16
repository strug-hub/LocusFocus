"use client";

import { alpha, styled } from "@mui/material";
import { DataGrid, DataGridProps, gridClasses } from "@mui/x-data-grid";

const LFDataTable: React.FC<DataGridProps> = (props) => (
  <StripedDataGrid
    initialState={{
      pagination: { paginationModel: { pageSize: 10 } },
    }}
    getRowClassName={(params) =>
      params.indexRelativeToCurrentPage % 2 === 0 ? "even" : "odd"
    }
    density="compact"
    showToolbar
    {...props}
  />
);

export default LFDataTable;

const ODD_OPACITY = 0.2;

const StripedDataGrid = styled(DataGrid)(({ theme }) => ({
  [`& .${gridClasses.row}.even`]: {
    backgroundColor: theme.palette.grey[200],
    "&:hover": {
      backgroundColor: alpha(theme.palette.primary.main, ODD_OPACITY),
      "@media (hover: none)": {
        backgroundColor: "transparent",
      },
    },
    "&.Mui-selected": {
      backgroundColor: alpha(
        theme.palette.primary.main,
        ODD_OPACITY + theme.palette.action.selectedOpacity
      ),
      "&:hover": {
        backgroundColor: alpha(
          theme.palette.primary.main,
          ODD_OPACITY +
            theme.palette.action.selectedOpacity +
            theme.palette.action.hoverOpacity
        ),
        // Reset on touch devices, it doesn't add specificity
        "@media (hover: none)": {
          backgroundColor: alpha(
            theme.palette.primary.main,
            ODD_OPACITY + theme.palette.action.selectedOpacity
          ),
        },
      },
    },
    ...theme.applyStyles("dark", {
      backgroundColor: theme.palette.grey[800],
    }),
  },
}));
