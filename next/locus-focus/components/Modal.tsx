import { Grid, Modal as MuiModal, SxProps } from "@mui/material";

interface ModalProps {
  children: React.ReactElement;
  onClose: () => void;
  open: boolean;
  sx?: SxProps;
}

const Modal: React.FC<ModalProps> = ({ open, onClose, children, sx }) => (
  <MuiModal
    open={open}
    onClose={onClose}
    sx={{
      display: "flex",
      justifyContent: "center",
      marginTop: 5,
      ...sx,
    }}
  >
    <Grid
      sx={{ backgroundColor: "white", borderRadius: 2 }}
      padding={4}
      maxWidth="800px"
      maxHeight="50%"
      container
      justifyContent="center"
      alignItems="center"
    >
      <Grid>{children}</Grid>
    </Grid>
  </MuiModal>
);

export default Modal;
