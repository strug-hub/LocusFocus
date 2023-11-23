function checkLocusInput(regiontext) {
    var build = d3.select("#coordinate").property("value");
    var errorSelect = d3.select("#locusErrorDiv");
    errorSelect.text("");
    if (regiontext !== "") {
        d3.json(`/regionCheck/${build}/${regiontext}`).then(response => {
            var message = response['response'];
            if (message !== "OK") {
                errorSelect.text(message);
            }
            else {
                // var genesButton = d3.select("#load-genes-button");
                // genesButton.property('value', `Load Genes in ${regiontext}`);
                loadGenes(build, regiontext);
            }
        })
    }
}

function checkSSInput(regiontext) {
    var build = d3.select("#coordinate").property("value");
    var errorSelect = d3.select("#locusErrorDiv");
    errorSelect.text("");
    if (regiontext !== "") {
        d3.json(`/regionCheck/${build}/${regiontext}`).then(response => {
            var message = response['response'];
            if (message !== "OK") {
                errorSelect.text(message);
            }
        })
    }
}

function checkNumSamplesInput(numsamples) {
    errordiv = d3.select("#numSamplesError-message");
    errordiv.text("");
    if (Number.isInteger(+numsamples) === false) {
        errordiv.text("Must be integer")
    }
}

/**
 * Validate input format for multiple region inputs.
 * No overlaps, no negative BP positions, and no unusual chromosome values.
 */
async function checkMultipleRegionsInput() {
    let regiontext = $("#multi-region").val();
    let errordiv = $("#multi-region-error");
    errordiv.text("");
    let build = d3.select("#coordinate").property("value");
    let regions = regiontext.trim().split("\n");
    let errors = [];

    for (let i = 0; i < regions.length; i++) {
        let region = regions[i];
        let response = await d3.json(`/regionCheck/${build}/${region}`);
        let message = response['response'];
        if (message !== "OK") {
            errors.push(`Line ${i+1}: ${message}`);
        }    
    }

    errordiv.text(errors.join("\n<br />\n"));
}