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
  return Math.abs((value - min)/(max - min));
}

function plot_heatmap(
  genes,
  tissues,
  SSPvalues,
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
  for (i = 0; i < tissues.length; i++) {
    for (j = 0; j < genes.length; j++) {
      newSSPvalues[i][j] = parseFloat(newSSPvalues[i][j]);
      if (newSSPvalues[i][j] === -2) {
        newSSPvalues[i][j] = -1;
      } else if (newSSPvalues[i][j] === -3) {
        newSSPvalues[i] = -1;
      } else if (newSSPvalues[i][j] > pmax) {
        pmax = newSSPvalues[i][j];
      }
    }
  }
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

  let color_cutoffs = [
    [-1, "rgb(105,105,105)"], // gray
    [0, "rgb(105,105,105)"],

    // Next trying to follow LD colors

    // dark blue to bright red HSV gradient in 5 steps:
    [0, "rgb(0, 0, 128)"], // dark blue (not significant)
    [1.3, "rgb(0, 0, 128)"], // dark blue

    [1.3, "rgb(0, 147, 142)"], // teal (p < 0.05)
    [5, "rgb(0, 147, 142)"],

    [5, "rgb(10, 166, 0)"], // green (p < 1e-5)
    [6, "rgb(10, 166, 0)"],

    [6, "rgb(185, 167, 0)"], // yellow/orange (p < 1e-6)
    [7.3, "rgb(185, 167, 0)"],

    [7.3, "rgb(204, 0, 24)"], // dark red (p < 5e-8)
  ].filter(([threshold, _]) => (threshold >= -1 && threshold <= pmax)) // remove colors that are out of bounds of data
  .map(([threshold, color]) => [normalize(threshold, -1, pmax), color]);

  color_cutoffs.push([1, color_cutoffs.slice(-1)[0][1]]);

  let data = [
    {
      z: newSSPvalues,
      x: genes.map((gene) => `${gene}`),
      y: tissues,
      //colorscale: 'Portland',
      type: "heatmap",
      name: "-log10(Simple Sum P-value)",
      hovertemplate:
        "Gene: %{x}" + "<br>Tissue: %{y}<br>" + `-log10(SS P-value): %{z}`,
      colorbar: {
        title: "-log10(Simple<br>Sum P-value)",
        dtick0: 0,
        dtick: 1,
        autotick: false,
      },
      colorscale: color_cutoffs,
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
