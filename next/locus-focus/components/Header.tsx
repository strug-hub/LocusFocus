"use client";

import React from "react";
import { AppBar, Grid, Toolbar, Typography } from "@mui/material";
import { NavLinkProps } from "./NavLink";
import { NavLink } from "@/components";

const HeaderNavLink = (props: NavLinkProps) => <NavLink thin {...props} />;

const Header: React.FC = () => {
  return (
    <AppBar position="static" color="primary" sx={{ marginBottom: 3 }}>
      <Toolbar alignItems="center" component={Grid} container>
        <Grid alignItems="center" spacing={2} container size={{ xs: 10 }}>
          <Grid marginRight={4}>
            <HeaderNavLink noDecoration href="/">
              <Typography
                sx={(theme) => ({
                  fontWeight: theme.typography.fontWeightLight,
                })}
                textAlign="center"
                variant="h5"
              >
                LocusFocus
              </Typography>
            </HeaderNavLink>
          </Grid>
          <Grid>
            <HeaderNavLink href="/">Colocalization</HeaderNavLink>
          </Grid>
          <Grid>
            <HeaderNavLink href="/set-based-test">Set-based Test</HeaderNavLink>
          </Grid>
          <Grid>
            <HeaderNavLink
              target="_blank"
              href="https://svatalog.research.sickkids.ca/"
            >
              GWAS SVatalog
            </HeaderNavLink>
          </Grid>
          <Grid>
            <HeaderNavLink
              href="https://locusfocus.readthedocs.io/en/latest/"
              target="_blank"
            >
              Documentation
            </HeaderNavLink>
          </Grid>
          <Grid>
            <HeaderNavLink href="/contact">Contact Us</HeaderNavLink>
          </Grid>
          <Grid>
            <HeaderNavLink href="/citation">Citation</HeaderNavLink>
          </Grid>
        </Grid>
        <Grid
          flexGrow={1}
          size={{ xs: 2 }}
          justifyContent="flex-end"
          container
          spacing={3}
        >
          <Typography sx={(theme) => ({ color: theme.palette.grey[300] })}>
            v1.6.0
          </Typography>
        </Grid>
      </Toolbar>
    </AppBar>
  );
};

export default Header;
