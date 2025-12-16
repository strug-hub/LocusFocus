import { defineConfig, globalIgnores } from "eslint/config";
import nextVitals from "eslint-config-next/core-web-vitals";
import nextTs from "eslint-config-next/typescript";
import hooks from "eslint-plugin-react-hooks";
import importOrder from "eslint-plugin-import";
import prettier from "eslint-config-prettier";

const eslintConfig = defineConfig([
  nextVitals,
  nextTs,
  // Override default ignores of eslint-config-next.
  globalIgnores([
    // Default ignores of eslint-config-next:
    ".next/**",
    "out/**",
    "build/**",
    "next-env.d.ts",
  ]),
  {
    plugins: {
      hooks: [hooks],
      importOrder: [importOrder],
      prettier: [prettier],
    },
    rules: {
      "no-console": [
        "warn",
        {
          allow: ["error"],
        },
      ],
      "dot-notation": "error",
      "react/no-unescaped-entities": "off",
      "react-hooks/rules-of-hooks": "error",
      "react-hooks/exhaustive-deps": "warn",
      "@typescript-eslint/no-explicit-any": "off",
      "@typescript-eslint/no-non-null-assertion": "off",
      "@typescript-eslint/no-this-alias": "off",
      "import/order": [
        "error",
        {
          pathGroups: [
            {
              pattern: "next/dynamic",
              group: "builtin",
              position: "before",
            },
            {
              pattern: "react",
              group: "builtin",
              position: "before",
            },
          ],
          pathGroupsExcludedImportTypes: ["react"],
        },
      ],
      "import/no-useless-path-segments": "error",
    },
  },
]);

export default eslintConfig;
