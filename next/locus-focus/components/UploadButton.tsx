import React from "react";
import { Box, Button, Tooltip } from "@mui/material";

interface UploadButtonProps {
  accept: string;
  disabled?: boolean;
  label: string;
  multiple: boolean;
  onChange: (e: React.ChangeEvent<HTMLInputElement>) => any;
  tooltipMessage?: string;
}

const UploadButton: React.FC<UploadButtonProps> = ({
  accept,
  disabled,
  label,
  multiple,
  onChange,
  tooltipMessage,
}) => (
  <Tooltip title={tooltipMessage}>
    <Box>
      <Button
        sx={{ margin: 1 }}
        variant="contained"
        component="label"
        disabled={disabled}
      >
        {label}
        <input
          hidden
          multiple={multiple}
          accept={accept}
          type="file"
          onChange={onChange}
        />
      </Button>
    </Box>
  </Tooltip>
);

export default UploadButton;
