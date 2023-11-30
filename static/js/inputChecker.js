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
 * Determine if two regions overlap. 
 * Assume proper format (chrom:start-end).
 * 
 * @param {string} regiontextA 
 * @param {string} regiontextB 
 */
function regionIsOverlapping(regiontextA, regiontextB) {
    let [chromA, startA, endA] = regiontextA.split(":").flatMap((s) => s.split("-"));
    let [chromB, startB, endB] = regiontextB.split(":").flatMap((s) => s.split("-"));
    if (chromA !== chromB) return false;
    
    startA = Number(startA);
    endA = Number(endA);
    startB = Number(startB);
    endB = Number(endB);

    return startA <= endB && startB <= endA;
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
    let regions = regiontext.trim().split("\n").filter(r => r !== "");
    if (regions.length === 0) return;
    let errors = [];

    for (let i = 0; i < regions.length; i++) {
        let region = regions[i];
        if (region === "") continue;
        let response = await d3.json(`/regionCheck/${build}/${region}`);
        let message = response['response'];
        if (message !== "OK") {
            errors.push(`Line ${i+1} ("${region}"): ${message}`);
        }    
    }

    if (errors.length === 0) {
        // overlap check
        for (let i = 0; i < regions.length; i++) {
            for (let j = i+1; j < regions.length; j++) {
                if (regionIsOverlapping(regions[i], regions[j])) {
                    errors.push(`Region on line ${i+1} ("${regions[i]}") overlaps wtih region on line ${j+1} ("${regions[j]}")`);
                }
            }
        }
    }


    errordiv.text(errors.join("\n<br />\n"));
}