"use client";

import React from "react";
import Link, { LinkProps } from "next/link";
import { Link as MuiLink } from "@mui/material";
import { usePathname } from "next/navigation";

export interface NavLinkProps extends LinkProps {
  children: React.ReactNode;
  noDecoration?: boolean;
  target?: string;
  thin?: boolean;
}

const NavLink: React.FC<NavLinkProps> = ({
  children,
  href,
  noDecoration,
  target,
  thin,
}) => {
  const pathname = usePathname();

  const active = pathname === href;

  return (
    <MuiLink
      sx={({ palette, typography }) => ({
        color: palette.primary.contrastText,
        textDecoration: active && !noDecoration ? "underline" : "inherit",
        fontWeight: thin
          ? typography.fontWeightLight
          : typography.fontWeightRegular,
        "&: hover": {
          textDecoration: !noDecoration ? "underline" : "inherit",
        },
      })}
      component={Link}
      target={target}
      href={href}
    >
      {children}
    </MuiLink>
  );
};

export default NavLink;
