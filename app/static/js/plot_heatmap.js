function copy(aObject) {
  if (!aObject) {
    return aObject;
  }

  let v;
  let bObject = Array.isArray(aObject) ? [] : {};
  for (const k in aObject) {
    v = aObject[k];
    bObject[k] = typeof v === "object" ? copy(v) : v;
  }

  return bObject;
}

function normalize(value, min, max) {
  return Math.abs((value - min) / (max - min));
}

function plot_heatmap(
  genes,
  tissues,
  SSPvalues,
  SSPvalues_secondary,
  image_width = 1080,
  image_height = 720,
  font_size = 12
) {
  // remember that these are -log10P-values
  // want the different negative p-value statuses to have white/black/grey colors
  // can't really do it that way b/c of the different cases that may occur
  // prior checks ensure that we have at least one positive -log10 SS Pvalue
  let pmax = 0;
  let newSSPvalues = copy(SSPvalues);
  let num_datasets = 0;
  for (i = 0; i < tissues.length; i++) {
    for (j = 0; j < genes.length; j++) {
      newSSPvalues[i][j] = parseFloat(newSSPvalues[i][j]);
      if ([-2, -3].includes(newSSPvalues[i][j])) {
        newSSPvalues[i][j] = -1;
      } else {
        if (newSSPvalues[i][j] > 0) {
          num_datasets += 1;
        }
        if (newSSPvalues[i][j] > pmax) {
          pmax = newSSPvalues[i][j];
        }
      }
    }
  }
  SSPvalues_secondary.forEach((p) => {
    if (p > 0) {
      num_datasets += 1;
    }
  });

  // let colorscale_exception_percentage = 0.05;
  // let new_minp =
  //   -1 *
  //   (colorscale_exception_percentage / (1 - colorscale_exception_percentage)) *
  //   pmax;
  // for (i = 0; i < tissues.length; i++) {
  //   for (j = 0; j < genes.length; j++) {
  //     if (newSSPvalues[i][j] === -1) {
  //       newSSPvalues[i][j] = new_minp;
  //     }
  //   }
  // }

  // Plotly normalizes values before applying color scale.
  // we have to normalize in advance to know what colors should be present

  let suggested_threshold =
    num_datasets > 0 ? -Math.log10(0.05 / num_datasets) : 1.3;

  let color_cutoffs = [
    [-1, "rgb(105,105,105)"], // gray
    [0, "rgb(105,105,105)"],

    // Next trying to follow LD colors

    // dark blue to bright red HSV gradient in 5 steps:
    [0, "rgb(0, 0, 128)"], // dark blue (not significant)
    [suggested_threshold, "rgb(0, 0, 128)"], // dark blue

    [suggested_threshold, "rgb(0, 147, 142)"], // teal
    [suggested_threshold + 1, "rgb(0, 147, 142)"],

    [suggested_threshold + 1, "rgb(10, 166, 0)"], // green
    [suggested_threshold + 2.5, "rgb(10, 166, 0)"],

    [suggested_threshold + 2.5, "rgb(185, 167, 0)"], // yellow/orange
    [suggested_threshold + 5, "rgb(185, 167, 0)"],

    [suggested_threshold + 5, "rgb(204, 0, 24)"], // dark red
  ].filter(([threshold, _]) => threshold >= -1 && threshold <= pmax); // remove colors that are out of bounds of data

  let norm_color_cutoffs = color_cutoffs.map(([threshold, color]) => [
    normalize(threshold, -1, pmax),
    color,
  ]);

  norm_color_cutoffs.push([1, color_cutoffs.slice(-1)[0][1]]);

  color_cutoffs = color_cutoffs.map(
    ([threshold, color]) => Math.round((threshold + Number.EPSILON) * 1e3) / 1e3
  );
  let data = [
    {
      z: newSSPvalues,
      x: genes.map((gene) => `${gene}`),
      y: tissues,
      //colorscale: 'Portland',
      type: "heatmap",
      name: "-log10(Simple Sum P-value)",
      hovertemplate:
        "Gene: %{x}" +
        "<br>Tissue: %{y}<br>" +
        `-log10(SS P-value): %{z: null}`,
      colorbar: {
        title: "-log10(Simple<br>Sum P-value)",
        // dtick0: 0,
        // dtick: 1,
        // autotick: false,
        tickmode: "array",
        tickvals: color_cutoffs,
        ticktext: color_cutoffs,
      },
      colorscale: norm_color_cutoffs,
    },
  ];

  let layout = {
    annotations: [],
    margin: {
      r: 50,
      t: 50,
      b: 125,
      l: 300,
    },
    autosize: false,
    width: image_width,
    height: image_height,
    font: { size: font_size },
    xaxis: {
      //uncomment to set a max range of visible genes and set `dragmode=pan` to pan
      //range: [0, Math.min(genes.length, 10)],
      dtick: 1,
    },
  };

  // Tried to add the SSPvalue numbers, but does not place correctly (they all go into the middle of the plot)
  //   for (var i = 0; i < tissues.length; i++) {
  //     for (var j = 0; j < genes.length; j++) {
  //       var currentValue = SSPvalues[i][j];
  //       if (currentValue != -1) {
  //         var textColor = "white";
  //       } else {
  //         var textColor = "black";
  //       }
  //       var result = {
  //         x: genes[j],
  //         y: tissues[i],
  //         text: SSPvalues[i][j],
  //         showarrow: false,
  //         font: {
  //           color: textColor,
  //         },
  //       };
  //       layout.annotations.push(result);
  //     }
  //   }

  Plotly.newPlot("heatmap", data, layout);
}
